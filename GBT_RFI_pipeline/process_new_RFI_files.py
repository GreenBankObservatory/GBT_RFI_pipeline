"""
..module:: process_new_RFI_files.py
    :synopsis: From the beginning of the RFI processing to the end, this script seeks out the unprocessed sdfits RFI files, determines which ones are new or old, 
    gleans the necessary information to process this file in GBTIDL, processes it in GBTIDL, and uploads that information to a mysql database. 
..moduleauthor:: Joy Skipper <jskipper@nrao.edu>
Code Origin: https://github.com/JoySkipper/GBT_RFI_pipeline
"""

import sys
import numpy as np
import os
import csv
import pandas
from rfitrends import GBT_receiver_specs
from rfitrends import RFI_input_for_SQL
#from rfitrends import connection_manager
import rfitrends
import time
import multiprocessing as mp
import subprocess
import argparse
import tempfile
import pathlib
from mysql import connector

class EmptyScans(Exception):
    pass

class BadIDLProcess(Exception):
    pass
class TimeoutError(Exception):
    pass


def determine_new_RFI_files(path_to_current_RFI_files: str,path_to_processed_RFI_files: str):
    """
    :param path_to_current_RFI_files: This is the path to all recent RFI files that have not been pushed to the archive, including those that have been processed into the database and those that haven't
    :param path_to_processed_RFI_files: This is the path that should contain any RFI files you don't want repeated
    :return: Returns RFI_files_to_be_processed, which is all the names of files that still need go through the processing script
    """
    print("beginning to determine files to be processed.")
    processed_RFI_files = []
    RFI_files_to_be_processed = []

    try:
        with open(path_to_processed_RFI_files+'files_not_able_to_be_processed.txt','r') as bad_list_file:
            bad_projects = bad_list_file.read().splitlines()
    except(IOError):
        bad_projects = []
    print("found bad projects, skipping "+str(len(bad_projects))+" bad projects.")

    # This for loop gathers a list of processed RFI files from the directory in which they're contained
    for processed_file in os.listdir(path_to_processed_RFI_files): 
        if processed_file.endswith(".txt") and processed_file != "URLS.txt": # We only care about the .txt files containing actual data of RFI
            processed_RFI_files.append(processed_file)
    print("found files already processed. Found "+str(len(processed_RFI_files))+" already processed files.")
    # This for loop goes through all the RFI files still in sdfits and that haven't been archived, finds those that needs to be processed, and appends them to RFI_files_to_be_processed        
    for current_RFI_file in os.listdir(path_to_current_RFI_files): 
        if current_RFI_file.startswith("TRFI") and not any(current_RFI_file in s for s in processed_RFI_files) and current_RFI_file not in bad_projects:
            RFI_files_to_be_processed.append(current_RFI_file)
    print("found files to be processed. Found "+str(len(RFI_files_to_be_processed))+" files to be processed.") 

    return(RFI_files_to_be_processed)

def determine_all_RFI_files(path_to_current_RFI_files:str):
    """
    :param path_to_current_RFI_files:This is the path to all recent RFI files that have not been pushed to the archive, including those that have been processed into the database and those that haven't
    :param path_to_processed_RFI_files: This is the path to write to all the RFI files
    :return: Returns RFI_files_to_be_processed, which is all the files that need to be processed
    """
    print("User has selected to process all files in the given directory.")
    RFI_files_to_be_processed = []
    # This for loop goes through all the RFI files still in sdfits and that haven't been archived, finds those that needs to be processed, and appends them to RFI_files_to_be_processed        
    for current_RFI_file in os.listdir(path_to_current_RFI_files): 
        if current_RFI_file.startswith("TRFI"):
            RFI_files_to_be_processed.append(current_RFI_file)
    print("Found files to be processed. Found "+str(len(RFI_files_to_be_processed))+" files to be processed.") 
    return(RFI_files_to_be_processed)


def read_header(file_to_be_processed: str, path_to_current_RFI_files: str):
    """
    :param file_to_be_processed: This is the string containing the name of the file that needs to be processed by the IDL processing script
    :param path_to_current_RFI_files: This is the string containing the path to all current RFI files, in which the files to be processed are contained
    :return: header, and line_to_start_reader; the header contains all the header information, and the line_to_start reader is the line in which the header ends, needed for reading in the data later
    """
    #line_to_start reader is counting the line in which the header ends, so that we can later read in the data of the file and skip the header
    line_to_start_reader = 0


    with open(path_to_current_RFI_files+file_to_be_processed+"/"+file_to_be_processed+".raw.vegas/"+file_to_be_processed+".raw.vegas.index") as f:
        header = {}
        file_index = f.readline()
        
        # The equal sign indicates that this is header information, such as filename = myfile.fits
        while(file_index):
            if "=" in file_index: 
                header_entry = file_index.strip().split("=")
                header[header_entry[0]] = header_entry[1].strip()
            # The [] signs indicate that this is a line that can be skipped
            elif "[" in file_index: 
                pass
            # Once we've reached the end of the header, we want to break the line.
            else: 
                break
                
            file_index = f.readline()
            line_to_start_reader += 1

    return(header,line_to_start_reader)

# This dictionary is used to determine the "ymax" of the graphs if one chooses to create them. It depends on the receiver
ymax_determiner = {
    "Rcvr_342":100,
    "Rcvr_450":100,
    "Rcvr_600":100,
    "Rcvr_800":100,
    "RcvrPF_2":100,
    "RcvrPF_1":100,
    "Rcvr1_2":100,
    "Rcvr2_3":10,
    "Rcvr4_6":10,
    "Rcvr8_10":10,
    "Rcvr12_18":10,
    "RcvrArray18_26":10,
    "Rcvr26_40":10,
    "Rcvr40_52":10,
    "Rcvr68_92":10,
    "RcvrArray75_115":10
}

def find_parameters_to_process_file(RFI_files_to_be_processed: list,path_to_current_RFI_files):
    """
    param: RFI_files_to_be_processed: List of all RFI files that need to be processed by the GBTIDL processing script
    param: path_to_current_RFI_files: String containing the path to all current RFI files, in which the files to be processed are contained
    return: data_to_process; which is a list of lists containing each file with the information needed to run the GBTIDL processing script
    """
    data_to_process = []
    
    for file_to_be_processed in RFI_files_to_be_processed:
        try:
            _, line_to_start_reader = read_header(file_to_be_processed, path_to_current_RFI_files)
        except FileNotFoundError:
            print("file not found. Skipping.")
            continue

        # Read the csv-non header portion in
        with open(path_to_current_RFI_files+file_to_be_processed+"/"+file_to_be_processed+".raw.vegas/"+file_to_be_processed+".raw.vegas.index") as f:
            data = pandas.read_csv(f,delimiter='\s+|\t+|\s+\t+|\t+\s+',skiprows=line_to_start_reader,engine='python')
        # If the source is unknown, this is a bad scan. Currently, the processing IDL code we are feeding this into cannot handle files with even one bad scan. 
        # This processing IDL code is in the works to be replaced with one that is more robust. 
        if "Unknown" in data["SOURCE"]:
            print("Unknown source, skipping.")
            continue
        max_scan_number = max(data["SCAN"]) #Scans are 1-indexed
        min_scan_number = min(data["SCAN"]) 
        number_of_feeds = max(data["FDNUM"])+ 1 #FDNUM is zero-indexed
        number_of_IFs = max(data["IFNUM"]) + 1 #IFNUM is zero-indexed

        # Getting the first part of the filename (i.e. TRFI_122001_PF1) and taking the last component (PF1) that actually names the receiver
        receiver_name = file_to_be_processed.split("_")[-1]

        # The receiver name is not consistent in it's naming scheme. This changes the receiver name to the GBT standardized frontend naming scheme
        verified_receiver_name = GBT_receiver_specs.FrontendVerification(receiver_name)
        if verified_receiver_name == 'Unknown':
            print("Unknown Receiver. Skipping.")
            continue
        scanlist = list(range(min_scan_number,max_scan_number))
        
        # Figure out the ymax value using the info from the dictionary above
        ymax = ymax_determiner[verified_receiver_name]


       

        # Creating a dictionary for each file, append to a list of dictionaries 
        oneline_data_to_be_processed = {
            "filename": file_to_be_processed, 
            "frontend": verified_receiver_name,
            "list_of_scans": scanlist,
            "number_of_feeds":number_of_feeds,
            "number_of_IFs":number_of_IFs,
            "ymax":ymax
        }

        # Appending oneline data to be processed to create a list of lists for the file info
        data_to_process.append(oneline_data_to_be_processed)

        print(str(file_to_be_processed)+" parameter data gleaned")
    return(data_to_process)
    

def analyze_file(file_to_process,output_directory):
    """
    param: file_to_process:: if the data has passed all checks up to this point, it is a dictionary containing metadata needed to process the RFI file.
    """
    if file_to_process['list_of_scans'] == []:
        raise(EmptyScans)
    # The parameters for running the process are different if the receiver is ka (26_40) so it needs to be called separately
    # Unfortunately IDL does not have an easy way to access module data, so we have to actually navigate to the appropriate directory
    path = str(pathlib.Path(__file__).parent.absolute())+'/'
    temp_path = tempfile.gettempdir()+'/'
    temp_file = open(temp_path+"temp_file.pro","w+")
    temp_file.write('.compile '+path+'scalUtils_wilsonedit.pro\n')
    temp_file.write('.compile '+path+'rfiDisplay_wilsonedit.pro\n')
    temp_file.write('.compile '+path+'process_file\n')
    temp_file.write('offline,  \''+str(file_to_process["filename"])+'\' \n')
    temp_file.write('process_file, '+str(file_to_process["list_of_scans"])+', fdnum='+str(file_to_process["number_of_feeds"])+'-1, ymax='+str(file_to_process["ymax"])+', ifmax = '+str(file_to_process["number_of_IFs"])+'-1, nzoom = 0, output_file=\''+output_directory+'\', temp_path=\''+temp_path+'\', /blnkChans, /makefile')
    if file_to_process['frontend'] == 'Rcvr26_40':
        temp_file.write(', /ka\n')
    else:
        temp_file.write('\n')
    temp_file.write('exit')
    temp_file.close()
    # Create a subprocess that calls the idl script that can process the file. 
    process = subprocess.Popen(['gbtidl', temp_path+'temp_file.pro'])
    # Wait 5 minutes (300 seconds) and if the file does not finish processing, kill it. It has entered an infinite loop
    # Since this process uses gettp, an official GBTIDL module, which does not return errors when it enters an infinite loop, this is the best we can do.
    try:
        print('Running in process', process.pid)
        process.wait(timeout=300)
        # Prints the status of whether or not the subprocess occured to a file. Because of the difficulties of output communications between the IDL subprocess
        # And the python process, it was best to communicate this status through a file
        subprocess_success = open(temp_path+"stat.txt","r").read().strip('\n')
        if subprocess_success == "bad_data":
            raise(BadIDLProcess)
    except subprocess.TimeoutExpired:
        print('Timed out - killing', process.pid)
        process.kill()
        raise(TimeoutError)

    # After we're done with getting the status, go ahead and remove the stat and temporary idl file
    def remove_file(temp_path,filename):
        if os.path.exists(temp_path+filename):
            os.remove(temp_path+filename)
    remove_file(temp_path,'stat.txt')
    remove_file(temp_path,'temp_file.pro')

def add_slash_if_needed(path):
    # Tired of forgetting to put '/' at the end of my path and having it break everything
    if path[-1] != '/':
        path += '/'
    return(path)

def main():
    parser = argparse.ArgumentParser(description="Processes new RFI files from the Green Bank Telescope and prints them as .txt files to the current directory")
    parser.add_argument("current_path",help="The path to the current RFI files, of which some will be the new files waiting to be processed")
    parser.add_argument("-skipalreadyprocessed",help="a flag to determine if you want to reprocess files that have already been processed or no. If this is selected, you must give the path to files that you know have already been processed.",type=str)
    parser.add_argument('-output_directory',help='The directory to which you want the data to be written',type=str)
    parser.add_argument("--upload_to_database",help="a flag to determine if you want to upload the txt files to a given database",action="store_true")
    parser.add_argument("-host_name",help="The host name of the machine containing the database to which you want to upload (required if you have selected -upload_to_database) ",type=str)
    parser.add_argument("-database_name",help="The name of the database to which you want to upload (required if you have selected -upload_to_database)",type=str)
    parser.add_argument("-main_table",help="The string name of the table to which you'd like to upload your clean RFI data (required if you have selected -upload_to_database)",type=str)
    parser.add_argument("-bad_table",help="The string name of the table to which you'd like to upload your flagged or bad RFI data (required if you have selected -upload_to_database)",type=str)
    
    args = parser.parse_args()
    path_to_current_RFI_files = args.current_path
    
    # If an output directory is provided, use that. Otherwise, output to the current directory. 
    if not args.output_directory:
        output_directory = "./"
    else: 
        output_directory = args.output_directory
        output_directory = add_slash_if_needed(output_directory)


    if args.upload_to_database: 
        if args.host_name is None or args.database_name is None or args.main_table is None or args.bad_table is None:
            parser.error("--upload_to_database requires -IP_address, -database_name, -main_table, and -bad_table.")  
        host_name = args.host_name
        database = args.database_name
        #connection_manager = rfitrends.connection_manager.connection_manager(host_name,database)
        username= os.environ['RFI_DB_USERNAME']
        password = os.environ['RFI_DB_PASSWORD']
        connector.connect(user=username, password=password,
                    host=host,
                    database=database)     
   
    if args.skipalreadyprocessed:
        skipalreadyprocessed = args.skipalreadyprocessed
        path_to_processed_RFI_files = add_slash_if_needed(skipalreadyprocessed)
        # If you specified 'output_directory' as the path to processed RFI files, then we need to check that you also specified the output_directory flag
        if path_to_processed_RFI_files == 'output_directory/':
            if args.output_directory is None: 
                # If you didn't, raise an error
                parser.error("Specifying -skipalreadyprocessed with 'output_directory' requires you to specify the -output_directory flag.")
            # If you did, then we can set the path to the output directory given
            else:
                path_to_processed_RFI_files = output_directory
        print(path_to_current_RFI_files)
        print(path_to_processed_RFI_files)
        # Regardless of where you get the path to processed RFI files, you need to then get the data to be processed from each file
        RFI_files_to_be_processed = determine_new_RFI_files(path_to_current_RFI_files,path_to_processed_RFI_files)
        # Get the data to be processed from each file
    else:
        RFI_files_to_be_processed = determine_all_RFI_files(path_to_current_RFI_files)
        # Get all data from the directory given of files to be processed
    data_to_be_processed = find_parameters_to_process_file(RFI_files_to_be_processed,path_to_current_RFI_files)
    # Go through each file and process it, and tallying the number of problem files as well
    problem_tally = 0
    
    with open(path_to_processed_RFI_files+'files_not_able_to_be_processed.txt','a+') as bad_list_file:
        for file_to_process in data_to_be_processed:
            print("processing file: "+str(file_to_process['filename']))
            try:
                analyze_file(file_to_process,output_directory)
            except(EmptyScans):
                problem_tally += 1
                bad_list_file.write(file_to_process['filename']+'\n')
                print("File had no scans. Skipping.")
                continue
            except(BadIDLProcess):
                problem_tally += 1
                bad_list_file.write(file_to_process['filename']+'\n')
                print("Was not able to IDL reduce file. Skipping.")
                continue
            except(TimeoutError):
                problem_tally += 1 
                bad_list_file.write(file_to_process['filename']+'\n')
                print("File processing timed out. Skipping.")

            print("file "+str(file_to_process['filename'])+" processed.")
        
    print("All new files processed and loaded as .txt files")
    # Let the user know how many bad files there were, if any:
    if problem_tally > 0:
        print(str(problem_tally)+" file out of "+str(len(data_to_be_processed))+" failed to process due to bad data.")
    else: 
        print("all files processed successfully")
    if args.upload_to_database: 
        # Find which file to be processed
        main_table = args.main_table
        bad_table = args.bad_table
        print("Uploading .txt files to database")
        filepaths_to_process = RFI_input_for_SQL.gather_filepaths_to_process(path_to_processed_RFI_files)
        RFI_input_for_SQL.upload_files(filepaths_to_process,connection_manager,main_table,bad_table)
        print("All files uploaded to database")
    
if __name__ == '__main__': 
    import ptvsd 
    main()