; 
;..module:: rfiDisplay_wilsonedit.pro
    ;:synopsis:: Processes RFI sdfits data and exports it to a .txt file. Initially written by Ron Maddalena and Jim Braats, it has been edited by Joy Skipper to accomodate error handling
;..moduleauthor:: Joy Skipper <jskipper@nrao.edu>, Ron Maddalena <rmaddale@nrao.edu>, Jim Braatz <jbraatz@nrao.edu>
;Code Origin: https://github.com/JoySkipper/GBT_RFI_pipeline

; In version 1.0 of GBT_RFI_Pipeline, the hope is to either go through this code, try to understand it, and document it, or to redo it in a well-documented fashion.

; makegif.pro by Jim Braatz
; 
; 


pro makegif,filename,notrim=notrim,reverse=reverse
    common gbtplot_common,mystate,xarray
    if n_elements(filename) eq 0 then filename = 'mygif.gif'
    print,'Making a GIF image in file ',filename
    !g.frozen = 1
    print, "frozen 3 rfiDisplay"
    widget_control,mystate.main,/show
    reshow
    if n_elements(reverse) ne 0 then begin
        tmp = !g.background
        !g.frozen = 1
        print, "frozen 4 rfiDisplay"
        !g.background = !g.foreground
        !g.foreground = tmp
        reshow
        end
        wait, 1
        spawn,'xwd -name "GBTIDL Plotter" -out temp.xwd'
        spawn,'convert temp.xwd gif:temp.gif'
        spawn,'rm temp.xwd',/sh
        if n_elements(notrim) ne 1 then begin
        jnk = query_gif('temp.gif',gif_info)
        x = gif_info.dimensions[0]-11
        y = gif_info.dimensions[1]-47
        s = "gifclip -q -i 8 52 "+strtrim(string(x),2)+" "+strtrim(string(y),2)+ $
        " temp.gif > "+strtrim(filename,2)
        spawn,s
        spawn,'rm temp.gif',/sh
        end else begin
        spawn,'mv temp.gif '+strtrim(filename,2)
    end
    if n_elements(reverse) ne 0 then begin
        tmp = !g.background
        !g.frozen = 1
        print, "frozen 5 rfiDisplay"
        !g.background = !g.foreground
        !g.foreground = tmp
        reshow
    end
end

pro displayRFI 
    hanning
    data0 = getdata()
    !g.frozen = 1
    print, "frozen 6 rfiDisplay"
    setdata, abs(convol(0.5 * data0 * !g.s[0].tsys/mean(data0), digital_filter(0.1, 1, 200, 64)))
    !g.s[0].units = "Jy"
end

function rfiFilename
    projID = strtrim(!g.s[0].projid,2)
    src = strtrim(!g.s[0].source,2)
    scan = string(!g.s[0].scan_number, format='("_s",i4.4)')
    feed = string(!g.s[0].feed, format='("_f",i3.3)')
    pol = strtrim(!g.s[0].polarization,2)
    polType = "_Linr"
    if (pol eq "LL" or pol eq "RR" or pol eq "RL" or pol eq "LR") then begin
        polType = "_Circ"
    endif
    az = string(floor(!g.s[0].azimuth+0.5), format='("_az",i3.3)')
    el = string(floor(!g.s[0].elevation+0.5), format='("_el",i3.3)')
    filename = projID + "_" + src + scan + feed + polType + az + el
    return, filename
end

pro getRFIScan, s1, plnum=plnum, ifnum=ifnum, fdnum=fdnum, instance=instance, nbox=nbox, tau=tau, ap_eff=ap_eff, $
        fltrParms=fltrParms, intnum=intnum, blnkWdth=blnkWdth, blnkChans=blnkChans, blnkFreqs=blnkFreqs, gain=gain

    if (n_elements(ifnum) eq 0) then ifnum = 0
    if (n_elements(plnum) eq 0) then plnum = 0
    if (n_elements(fdnum) eq 0) then fdnum = 0
    if (n_elements(instance) eq 0) then instance = 0
    if (n_elements(nbox) eq 0) then nbox=64
    if (n_elements(blnkWdth) eq 0) then blnkWdth = 0
    if (n_elements(fltrParms) eq 0) then begin
        fltrParms=[0.1,200,64]
    end

    fltrCoeff = digital_filter(fltrParms[0], 1, fltrParms[1], fltrParms[2])

    if (n_elements(gain) eq 0 or gain[0] eq 0) then begin
        !g.frozen = 1
    print, "frozen 7 rfiDisplay"
        gettp,s1,plnum=plnum,ifnum=ifnum,fdnum=fdnum,instance=instance,cal_state=1,intnum=intnum, /quiet & hanning

        nchans = n_elements(getdata(0))
        if keyword_set(blnkFreqs) then begin
            blankFreqs, blnkWdth
        end
        if keyword_set(blnkChans) then begin
            blankChans, blnkWdth
        end
        on = getdata()

        ; on[5000]=10.

        !g.frozen = 1
        print, "frozen 8 rfiDisplay"
        gettp,s1,plnum=plnum,ifnum=ifnum,fdnum=fdnum,instance=instance,cal_state=0,intnum=intnum, /quiet & hanning 

        if keyword_set(blnkFreqs) then begin
            blankFreqs, blnkWdth
        end
        if keyword_set(blnkChans) then begin
            blankChans, blnkWdth
        end
        off = getdata()
        ; off[5000]=10.
        
        setdata, 0.5*!g.s[0].mean_tcal*abs(convol( (on+off)/(doboxcar1d(on,nbox,/nan,/edge_truncate)-doboxcar1d(off,nbox,/nan,/edge_truncate))-1., fltrCoeff))
 
    endif else begin
        !g.frozen = 1
        print, "frozen 9 rfiDisplay"
        gettp,s1,plnum=plnum,ifnum=ifnum,fdnum=fdnum,instance=instance,intnum=intnum, /quiet & hanning

        if keyword_set(blnkFreqs) then begin
            blankFreqs, blnkWdth 
        end
        if keyword_set(blnkChans) then begin
            blankChans, blnkWdth
        end
        setdata, gain*abs(convol( getdata(), fltrCoeff))

    end
    nchans=n_elements(getdata())

    if n_elements(tau) eq 1 then begin
        if tau eq -1 then tau = getForecastedTau(dateToMJD(!g.s[0].timestamp),!g.s[0].center_frequency/1.e6)
    endif
    atten = getTau(!g.s[0].center_frequency/1.e6, coeffs=tau)

    eff = getApEff(!g.s[0].elevation, !g.s[0].center_frequency/1.e6, coeffs=ap_eff)
    Ta2Flux, tau=atten, ap_eff=eff

    replace, 0, fltrParms[2]+0.03*nchans-1, /blank
    replace, nchans-fltrParms[2]-0.005*nchans, nchans-1, /blank
end

pro getRFITP, s1, plnum=plnum, ifnum=ifnum, fdnum=fdnum, instance=instance, fltrParms=fltrParms, intnum=intnum, blnkWdth=blnkWdth, blnkChans=blnkChans, blnkFreqs=blnkFreqs

    if (n_elements(ifnum) eq 0) then ifnum = 0
    if (n_elements(plnum) eq 0) then plnum = 0
    if (n_elements(fdnum) eq 0) then fdnum = 0
    if (n_elements(instance) eq 0) then instance = 0
    if (n_elements(blnkWdth) eq 0) then blnkWdth = 0
    if (n_elements(fltrParms) eq 0) then begin
        fltrParms=[0.1,200,64]
    end

    !g.frozen = 1
    print, "frozen 10 rfiDisplay"
    gettp,s1,plnum=plnum,ifnum=ifnum,fdnum=fdnum,instance=instance,intnum=intnum, /quiet & hanning
    nchans=n_elements(getdata())

    if keyword_set(blnkFreqs) then begin
        blankFreqs, blnkWdth
    end
    if keyword_set(blnkChans) then begin
        blankChans, blnkWdth
    end
    on = getdata()

    setdata, abs(convol( on, digital_filter(fltrParms[0], 1, fltrParms[1], fltrParms[2])))
    replace, 0, fltrParms[2]+0.03*nchans-1, /blank
    replace, nchans-fltrParms[2]-0.005*nchans, nchans-1, /blank

end

pro blankFreqs, width

    if (width le 0) then begin
        return
    end

    openr, lun, "/users/rmaddale/mypros/blankFreqs.dat", /get_lun
    nchans=n_elements(getdata())
    while (eof(lun) ne 1) do begin
        on_ioerror, bad_rec
        rcvr=string(20)
        readf, lun, freq, rcvr
        rcvr=strtrim(rcvr,2)
        bchan = xtochan(freq)-width/2
        echan=bchan+width
        if (rcvr eq '*' || rcvr eq !g.s[0].frontend) and (echan ge 0) and (bchan lt nchans) then begin
            replace, bchan, echan, /blank
        end
        bad_rec:
    endwhile
    free_lun, lun

end

pro blankChans, width

    if (width le 0) then begin
        return
    end

    openr, lun, "/users/rmaddale/mypros/blankChans.dat", /get_lun
    nchans=n_elements(getdata())
    while (eof(lun) ne 1) do begin
        on_ioerror, bad_rec
        ; readf, lun, bw, n, c, rcvr, format='(f5, i7, i7, a20)'
        rcvr=string(20)
        readf, lun, bw, n, c, rcvr
        rcvr=strtrim(rcvr,2)
        n=long(n)
        c=long(c)
        if ((rcvr eq '*' || rcvr eq !g.s[0].frontend) and abs(bw - !g.s[0].bandwidth/1.e6) lt 0.1 and n eq nchans) then begin
           bchan=c-width/2
           echan=bchan+width
           if (echan ge 0) and (bchan lt nchans) then begin
               print, 'Blanking Channels:', bchan, echan
               replace, bchan, echan, /blank
           end
        end
        bad_rec:
    endwhile
    free_lun, lun

end


pro flagFreqs
    openr, lun, "/users/rmaddale/mypros/flagFreqs.dat", /get_lun
    iLines=0
    while (eof(lun) ne 1) do begin
        on_ioerror, bad_rec
        rcvr=string(20)
        readf, lun, freq, rcvr
        rcvr=strtrim(rcvr,2)
        if (rcvr eq '*' || rcvr eq !g.s[0].frontend) then begin
            vline, freq, label=strtrim(string(freq,format='(f10.6)'),2), ylabel=0.95-iLines*0.05, /noshow, /ynorm
            iLines=iLines+1
            if (iLines GT 5) then begin
                iLines=0
            end
        end
        bad_rec:
    endwhile
    free_lun, lun
    !g.frozen = 1
    print, "frozen 11 rfiDisplay"
    reshow
end

pro zoomGifs, filename, zfactor=zfactor, nzoom=nzoom

    if (n_elements(nzoom) eq 0) then nzoom = 0
    if (n_elements(zfactor) eq 0) then zfactor = 10

    makegif, filename + "_z00.gif", reverse=1
    yr = getyrange()
    for iz=1, nzoom do begin
        z = string(iz,format='("_z",i2.2)')
        sety, yr[0]/(zfactor^iz), yr[1]/(zfactor^iz)
        makegif, filename + z + ".gif", reverse=1
    end
end

pro writeDC, filename, buffer=buffer, freqMin=freqMin, freqMax=freqMax, output_file=output_file

    if (n_elements(buffer) eq 0) then buffer = 0

    nchans=n_elements(getdata(buffer))

    ; Open the TXT file that will contain the data
    print,'Writing ASCII data to ', filename + '.txt, using writeDC'
    openw, lun, output_file+filename + '.txt', /GET_LUN

    ; Add a header to the text file
    printf, lun, "################ HEADER #################"
    printf, lun, "# projid: ", !g.s[buffer].projid
    printf, lun, "# date: ", !g.s[buffer].date
    printf, lun, "# utc (hrs): ", !g.s[buffer].utc/3600.
    printf, lun, "# mjd: ", !g.s[buffer].mjd
    printf, lun, "# lst (hrs): ", !g.s[buffer].lst/3600.
    printf, lun, "# scan_number: ", !g.s[buffer].scan_number
    printf, lun, "# frontend: ", !g.s[buffer].frontend
    printf, lun, "# feed: ", !g.s[buffer].feed
    printf, lun, "# polarization: ", !g.s[buffer].polarization
    printf, lun, "# backend: ", !g.s[buffer].backend
    printf, lun, "# if number: ", !g.s[buffer].if_number
    printf, lun, "# exposure (sec): ", !g.s[buffer].exposure
    printf, lun, "# tsys (K): ", !g.s[buffer].tsys
    printf, lun, "# frequency_type: ", !g.s[buffer].frequency_type
    printf, lun, "# frequency_resolution (MHz): ", !g.s[buffer].frequency_resolution/1.e6
    printf, lun, "# source: ", !g.s[buffer].source
    printf, lun, "# azimuth (deg): ", !g.s[buffer].azimuth
    printf, lun, "# elevation (deg): ", !g.s[buffer].elevation
    printf, lun, "# units: ", !g.s[buffer].units
    printf, lun, "################   Data  ################"
    printf, lun, "# Channel Frequency(GHz)  Intensity(Jy)"

    setxunit, 'GHz'
    setframe, 'TOPO'
    bdrop, 0
    edrop, 0

    ; Retrieve X and Y arrays.  Sort data according to frequency.
    !g.frozen = 1
    print, "frozen 12 rfiDisplay"
    unzoom
    show, buffer
    x = getxarray()
    y = getyarray()
    if (n_elements(freqMin) eq 0) then freqMin=min(x)
    if (n_elements(freqMax) eq 0) then freqMax=max(x)
    for ic=long(0), long(nchans-1) do begin
        if (x[ic] GT freqMin AND x[ic] LT freqMax) then begin
            printf, lun, ic+1, x[ic], y[ic], format='(i10, 2d15.6)'
        endif
    end
    free_lun, lun
end

function calseq,scan,tcold=tcold,ifnum=ifnum,plnum=plnum,fdnum=fdnum
;
;;Computes gain for W-band (modified from version provided by D. Frayer).
;
;;Inputs:
;;scan = auto calseq scan
;;tcold = effective temperature of cold load (e.g., 50K)
;;ifnum = IFnum of spectral window
;;plnum = pol-number
;;fdnum = beam-number
;;  
    if (n_elements(tcold) eq 0) then tcold = 54.

    !g.frozen = 1
    print, "frozen 13 rfiDisplay"
    gettp,scan,plnum=plnum,fdnum=fdnum,ifnum=ifnum,quiet=1,wcalpos='Observing'
    vsky=getdata(0)
    twarm=!g.s[0].twarm
    gettp,scan,plnum=plnum,fdnum=fdnum,ifnum=ifnum,quiet=1,wcalpos='Cold1'
    vcold1=getdata(0)
    gettp,scan,plnum=plnum,fdnum=fdnum,ifnum=ifnum,quiet=1,wcalpos='Cold2'
    vcold2=getdata(0)

    if (!g.s[0].feed eq 1) then begin
        gain=(twarm-tcold)/(vcold2-vcold1)
    endif else begin 
        gain=(twarm-tcold)/(vcold1-vcold2)
    end

    return, gain
end


FUNCTION rfiScans_Mod, scanList, ifmax=ifmax, fdnum=fdnum, intnum=intnum, nzoom=nzoom, ymax=ymax, $
                zfactor=zfactor, instance=instance, makegifs=makegifs, makefile=makefile, ka=ka, $
                nbox=nbox, tau=tau, ap_eff=ap_eff, fltrParms=fltrParms, blnkWdth=blnkWdth, $
                blnkChans=blnkChans, blnkFreqs=blnkFreqs, flagFreqs=flagFreqs, colors=colors, pols=pols, $
                calseqList=calseqList, output_file=output_file

; Necessary parameters: 
;   scanList = A scan number or a list of scan numbers.  Examples: 5 or [1,3,5]

; Data filtering:
;   Parameters with default values:
;       fdnum = feed to process.  Default=0.  Ignored if '/KA' is specified.
;       intnum = integration number.  Default=all integrations.
;       instance = which instance of a redundant scan number is desired.  Default is the 1st (0) instance.
;       pols = 0 or 1, the polarization you want to use.  If not set (the default), will average both 
;            polarizations.
;   Flags: 
;       /ka = a flag that specifies whether to process as if this were Ka-band, a single-pol, dual feed 
;       receiver.

; Specifying plot, GIF parameters 
;   Parameters with default values:
;       nzoom = number of zoom levels.  Default = 0 (1 zoom level).  Used only if /makegif is specified.
;       zfactor = factor between each zoom level.  Default = 10.  Used only if /makegif is specified.
;       ymax = maximum y-axis for 1st display.  Default is auto range.
;   Flags: 
;       /colors = a flag that specifies that each spectral window should be colored differently.
;       /makegifs = a flag that specifies that GIF files are to be produced

; Calibration specifics:
;   Parameters with default values:
;       calseqList = A list of scans to use for W-band chopper-wheel calibration.  The length of the list  
;              must be the same as the length of scanList.  The first scan in scanList will be calibrated  
;              using the first scan in calseqList.  If not specified (the default), the routine will  
;              assume a noise diode exists and calibrate using the noise diode on & off phases.
;       nbox = number of channels to smooth CalOn and CalOff.  Default=64.  Ignored if calseqList is 
;              specified.
;       tau = atmospheric opacity (see getTau in scalUtils.pro for details).  Default = 0 (no correction).  
;               If -1, then the algorithm will use the value from the archive of weather forecasts (see   
;               getForecastedTau).  Only practical for frequencies above about 5 GHz.
;       ap_eff = aperture efficiency (see getApEff in scalUtils.pro for details).  Default = standard 
;               model

; Algorithm details:
;   Parameters with default values:
;       fltrParms = three coefficient to pass to IDL routine digital_filter to create a high-pass filter.  
;             Defaults = [0.1,200,64], which seems to be a reasonable compromise in almost all cases
;       blnkWdth = width of blanking window in channels for any frequencies or channels listed in the 
;            blankFreqs.dat or blankChans.dat files.  If <= zero, the default, no channels/frequencies
;            will be blanked.
;       ifmax = maximum possible number of spectral windows to process.  Probably never needs to be 
;            specified since the default (=7, for 8 windows), is sufficient for all known cases.
;   Flags: 
;       /makefile = a flag that specifies an ASCII file is to be produced
;       /flagFreqs = a flag that specifies whether frequencies in the flagFreqs.dat file are to be 
;            'flagged'
;       /blnkChans = a flag that specifies the channels listed in the blankChans.dat file will be blanked.  
;            You must also specify a value for blankWdth that is greater than zero.
;       /blnkFreqs = a flag that specifies the channels listed in the blankFreqs.dat file will be blanked. 
;            You must also specify a value for blankWdth that is greater than zero.
;       Note that you can specify  /flagFreqs, /blnkChans, and /blnkFreqs at the same time.
    
    CATCH, Error_status
    if (Error_status NE 0) and !Error_state.msg eq "Attempt to subscript SORTORDER with IC is out of range." THEN BEGIN
        print, 'Error Index: ', Error_status
        print, 'Error Message: ', !Error_state.msg
        RETURN, "bad_data"
        
        
    endif else if (Error_status NE 0) then begin
        print, 'Error Index: ', Error_status
        print, 'Error Message: ', !Error_state.msg
        RETURN, "bad_data"
        
        
    endif
    

    if (n_elements(ifmax) eq 0) then ifmax = 7
    if (n_elements(fdnum) eq 0) then fdnum = 0
    if (n_elements(nzoom) eq 0) then nzoom = 0
    if (n_elements(zfactor) eq 0) then zfactor = 10
    if (n_elements(instance) eq 0) then instance = 0
    if (n_elements(blnkWdth) eq 0) then blnkWdth = 0
    if (n_elements(pols) eq 0) then pols = -1

    nscans = n_elements(scanList)
    nArrays = nscans*(ifmax+1)

    ; Create arrays to hold sorted X and Y arrays, etc.
    xsort = make_array(nArrays, 262144, /float)
    ysort = make_array(nArrays, 262144, /float)
    xmin = make_array(nArrays, /float)
    xmax = make_array(nArrays, /float)
    freqmin = make_array(nArrays, /float)
    freqmax = make_array(nArrays, /float)
    c0 = make_array(nArrays, /float)
    c1 = make_array(nArrays, /float)
    ifSort = make_array(nArrays, /long)
    ifs = make_array(nArrays, /long)
    cntrFreq = make_array(nArrays, /float)

    ; Create DC's to hold the individual plots
    dc_Arr=replicate({spectrum_struct}, nArrays)

    ; FLIP holds the sense of the X Frequqncy axis vs the x channel axis
    flip=0

    ; NPLOTS = number of plots that are to be generated
    nplots=0

    ; Get the first scan so as to determine a filename, number of channels
    if keyword_set(ka) then begin
        !g.frozen = 1
        print, "frozen 14 rfiDisplay"
        gettp, scanList(0), plnum=0, fdnum=1, instance=instance, /quiet
    endif else begin
        !g.frozen = 1
        print, "frozen 15 rfiDisplay"
        gettp, scanList(0), plnum=0, fdnum=fdnum, instance=instance, /quiet
    end

    setxunit, 'GHz'
    setframe, 'TOPO'
    bdrop, 0
    edrop, 0
    clear & freeze
    
    filename = rfiFilename()
    nchans=n_elements(getdata())

    for idScan=0, n_elements(scanList)-1 do begin

        s1 = scanList(idScan)
        info = scan_info(s1)   
        
        print, "Scan:", s1, "   Bandwith(MHz):", info.bandwidths[instance]/1e6, "   Channels:", info[instance].n_channels[0]

        ifs[nplots] = info[instance].n_ifs - 1

        ; calScan = calseqList(idScan)
        if (n_elements(calseqList) ne 0) then begin
            calScan = calseqList[idScan]
        end

        ; No longer can we assume the observations will be in increasing frequency.
        ; So, creating a mapping between ifNum and a sorted array of center frequencies.
        for ifTest = 0, ifs[nplots] do begin
            if keyword_set(ka) then begin
                !g.frozen = 1
                print, "frozen 16 rfiDisplay"
                gettp, s1, plnum=0, fdnum=1, ifnum=ifTest, intnum=0, instance=instance, /quiet
            endif else begin
                !g.frozen = 1
                print, "frozen 17 rfiDisplay"
                gettp, s1, plnum=0, fdnum=fdnum, ifnum=ifTest, intnum=0, instance=instance, /quiet
            end
            cntrFreq(ifTest) = chantox(0)
        end
        ifSort = sort(cntrFreq[0:ifs[nplots]])

        for ifTest = 0, ifs[nplots] do begin

            ifn = ifSort(ifTest)

            print, "Processing Scan:", s1, " IFNum:", ifn

            ; Average the polarizations
            sclear
            gain = 0
            case pols of
                0 : begin
                    if keyword_set(ka) then fdnum=1
                    if (n_elements(calseqList) ne 0) then gain = calseq(calScan, plnum=0, ifnum=ifn, fdnum=fdnum)

                    getRFIScan, s1, plnum=0, ifnum=ifn, fdnum=fdnum, intnum=intnum, instance=instance, nbox=nbox, $
                            tau=tau, ap_eff=ap_eff, fltrParms=fltrParms, blnkWdth=blnkWdth, blnkChans=blnkChans, blnkFreqs=blnkFreqs, gain=gain & accum
                    end
                1 : begin
                    if keyword_set(ka) then fdnum=0
                    if (n_elements(calseqList) ne 0) then gain = calseq(calScan, plnum=1, ifnum=ifn, fdnum=fdnum)

                    getRFIScan, s1, plnum=1, ifnum=ifn, fdnum=fdnum, intnum=intnum, instance=instance, nbox=nbox, $
                            tau=tau, ap_eff=ap_eff, fltrParms=fltrParms, blnkWdth=blnkWdth, blnkChans=blnkChans, blnkFreqs=blnkFreqs, gain=gain & accum
                    end
                else: begin
                    if keyword_set(ka) then fdnum=1
                    if (n_elements(calseqList) ne 0) then gain = calseq(calScan, plnum=0, ifnum=ifn, fdnum=fdnum)
                    ; print, gain

                    getRFIScan, s1, plnum=0, ifnum=ifn, fdnum=fdnum, intnum=intnum, instance=instance, nbox=nbox, $
                            tau=tau, ap_eff=ap_eff, fltrParms=fltrParms, blnkWdth=blnkWdth, blnkChans=blnkChans, blnkFreqs=blnkFreqs, gain=gain & accum

                    if keyword_set(ka) then fdnum=0
                    if (n_elements(calseqList) ne 0) then gain = calseq(calScan, plnum=1, ifnum=ifn, fdnum=fdnum)

                    ; gain2 = doboxcar1d(gain[floor(nchans*0.2):floor(nchans-nchans*0.2)], 100, /nan, /edge_truncate)
                    ; mn = mean(gain2, /nan)
                    ; print, max(gain2, /nan), min(gain2, /nan), mn, median(gain), stddev(gain2, /nan)
                    ; print, max(gain2, /nan)/mn - 1, min(gain2, /nan)/mn - 1, mn, median(gain)/mn - 1, stddev(gain2, /nan)/mn

                    getRFIScan, s1, plnum=1, ifnum=ifn, fdnum=fdnum, intnum=intnum, instance=instance, nbox=nbox, $
                            tau=tau, ap_eff=ap_eff, fltrParms=fltrParms, blnkWdth=blnkWdth, blnkChans=blnkChans, blnkFreqs=blnkFreqs, gain=gain & accum
                    end
            endcase
            ave, /quiet

            ; show

            ; Copy DC for later plotting
            tmp = dc_Arr[nplots]
            data_copy, !g.s[0], tmp
            dc_Arr[nplots] = tmp

            ; Must now worry about the 1500 MHZ mode of VEGAS only having a 1250 MHz usable range
            chanFirst = 0
            chanLast = nchans
            if !g.s[0].bandwidth GT 1499.e6 then begin
                if !g.s[0].frequency_interval LT 0 then begin
                    chanLast = nchans+250e6/!g.s[0].frequency_interval
                endif else begin
                    chanFirst = 250e6/!g.s[0].frequency_interval
                end
            end
            ; print, chanFirst, chanLast, nchans, chantox(0), chantox(chanFirst), chantox(chanLast), chanTox(nchans-1)

            ; Determine frequency limits for each plot
            f0=chantox(chanFirst)
            fnchans=chantox(chanLast-1)
            xmin[nplots] = f0
            xmax[nplots] = fnchans
            if f0 GT fnchans then begin
                flip=1
                xmin[nplots] = fnchans
                xmax[nplots] = f0
            end
        nplots=nplots+1
        end
    end

    nplots=nplots-1
    print, "Minimum window frequencies:", xmin[0:nplots]
    print, "Maximum window frequencies:", xmax[0:nplots]

    ; Determine order in which to print out the IF's
    printOrder = sort(xmin[0:nplots])

    ; Determine cutoff frequencies for overlapping areas.  
    ; Cutoffs are at the 30%/70% of the plot overlaps so as to reduce the possibility of 
    ; edge effects.  Use FLIP to determine how the 30/70 is to be taken
    freqMin[0] = xmin[printOrder[0]]
    freqMax[nplots] = xmax[printOrder[nplots]]
    for ifn = 0, nplots-1 do begin
        if flip EQ 1 then begin
            freqMax[ifn] = (5*xmax[printOrder[ifn]] + 5*xmin[printOrder[ifn+1]])/10.
            freqMin[ifn+1] = freqMax[ifn]
        endif else begin
            freqMax[ifn] = (5*xmax[printOrder[ifn]] + 5*xmin[printOrder[ifn+1]])/10.
            freqMin[ifn+1] = freqMax[ifn]
        end
    end

    ; Find the chanel numbers corresponding to each plots max/min frequency
    freex
    for ifn = 0, nplots do begin
        c0(ifn)=xtochan(freqMin(ifn),dc=dc_Arr[printOrder[ifn]])
        c1(ifn)=xtochan(freqMax(ifn),dc=dc_Arr[printOrder[ifn]])
        if c0(ifn) GT c1(ifn) then begin
            c0(ifn)=xtochan(freqMax(ifn),dc=dc_Arr[printOrder[ifn]])
            c1(ifn)=xtochan(freqMin(ifn),dc=dc_Arr[printOrder[ifn]])
        end 
    end

    ; Blank out the data that is beyond the overlaps.   Then, plot up the various DC's 
    !g.frozen = 1
    print, "frozen 18 rfiDisplay"
    set_data_container, dc_Arr[0]
    replace, 0, c0(0), /blank
    replace, c1(0), nchans-1, /blank
    show

    clrs = [!white, !red, !green, !blue, !magenta, !orange, !cyan, !yellow]
    iclr = 0
    for ifn = 0, nplots do begin
        !g.frozen = 1
        print, "frozen 19 rfiDisplay"
        set_data_container, dc_Arr[ifn]
        replace, 0, c0(ifn), /blank
        replace, c1(ifn), nchans-1, /blank
        oshow, color=clrs[iclr]
        if keyword_set(colors) then begin
            iclr=iclr+1
            if (iclr eq n_elements(clrs)) then begin
                iclr = 0
            endif
        endif
    end
    unzoom

    setxunit, 'GHz'
    setframe, 'TOPO'
    bdrop, 0
    edrop, 0

    if (n_elements(ymax) ne 0) then begin
        sety, -ymax/10., ymax
    endif
    setx, min(xmin[0:nplots]), max(xmax[0:nplots])

    ; Draw vertical lines
    if keyword_set(flagFreqs) then begin
        flagFreqs
    end

    ; Create GIF files.  Perform muiltiple zooms
    if keyword_set(makegifs) then begin
        zoomGifs, filename, zfactor=zfactor, nzoom=nzoom
     end

    ; Do to some vagaries with GBTIDL and DC's and plots, we must reprocess 
    ; evreything to generate the ASCII file
    if keyword_set(makefile) then begin 

        ; Open the TXT file that will contain the data
        print,'Writing ASCII data to ', output_file+filename + '.txt, using rfiscans_mod fxn'
        openw, lun, output_file+filename + '.txt', /GET_LUN

        ; Add a header to the text file
        printf, lun, "################ HEADER #################"
        printf, lun, "# projid: ", !g.s[0].projid
        printf, lun, "# date: ", !g.s[0].date
        printf, lun, "# utc (hrs): ", !g.s[0].utc/3600.
        printf, lun, "# mjd: ", !g.s[0].mjd
        printf, lun, "# lst (hrs): ", !g.s[0].lst/3600.
        printf, lun, "# scan_numbers: ", scanList
        printf, lun, "# frontend: ", !g.s[0].frontend
        printf, lun, "# feed: ", !g.s[0].feed
        printf, lun, "# polarization: ", !g.s[0].polarization
        printf, lun, "# backend: ", !g.s[0].backend
        printf, lun, "# exposure (sec): ", !g.s[0].exposure
        printf, lun, "# tsys (K): ", !g.s[0].tsys
        printf, lun, "# frequency_type: ", !g.s[0].frequency_type
        printf, lun, "# frequency_resolution (MHz): ", !g.s[0].frequency_resolution/1.e6
        printf, lun, "# source: ", !g.s[0].source
        printf, lun, "# azimuth (deg): ", !g.s[0].azimuth
        printf, lun, "# elevation (deg): ", !g.s[0].elevation
        printf, lun, "# units: ", !g.s[0].units
        printf, lun, "################   Data  ################"
        printf, lun, "# Window   Channel Frequency(GHz)  Intensity(Jy)"

        setxunit, 'GHz'
        setframe, 'TOPO'
        bdrop, 0
        edrop, 0

        for ifn = 0, nplots do begin
            set_data_container, dc_Arr[ifn]

            ; Retrieve X and Y arrays.  Sort data according to frequency.
            !g.frozen = 1
            print, "frozen 19 rfiDisplay"
            unzoom
            show
            setx, min(xmin), max(xmax)

            x = getxarray()
            y = getyarray()
            sortOrder = sort(x)
            ; print, x[0], x[nchans-1]

            
            for ic=long(0), long(nchans-1) do begin
                xsort[ifn,ic] = x[sortOrder[ic]]
                ysort[ifn,ic] = y[sortOrder[ic]]    
            end
        end

        ; print, printOrder, n_elements(xsort)

        ; Print out the data values.  Use the sorted order of IF's, remove overlapping
        ; frequency channels
        for ifn = 0, nplots do begin
            ifOrdered = printOrder[ifn]
            ; print, ifn, ifOrdered, freqMin(ifn), freqMax(ifn), xsort[ifOrdered,0], xsort[ifOrdered,nchans-1]
            for ic=long(0), long(nchans-1) do begin
                xvalue = xsort[ifOrdered,ic]
                if (xvalue GT freqMin[ifn] AND xvalue LT freqMax[ifn]) then begin
                    printf, lun, ifn+1, ic+1, xsort[ifOrdered,ic], ysort[ifOrdered,ic], format='(2i10, 2d15.6)'
                endif
            end
        end
        free_lun, lun

     end

    data_free, dc_Arr
    emptystack
    clear
    sclear
    print, 'Data Successfully written.'
    return, "good_data"
end

