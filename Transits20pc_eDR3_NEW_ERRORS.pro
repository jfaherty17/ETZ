pro transits_errors,nums

;nums is the total number of iterations used to determine uncertainty range in ETZ (and rETZ) entry and exit time

path='/Users/jfaherty/Dropbox/For-Lisa/'
;;;-------------------------------------
;;;SETTING UP FOR PLOTTING
;;;-------------------------------------
set_plot,'ps'
!p.font=0
device,helvetica=1
device,isolatin1=1
@/Users/jfaherty/Dropbox/CPAPIR/symbols_ps_kc
@/Users/jfaherty/Dropbox/CPAPIR/colors_kc


;;;-------------------------------------
;;;READING IN ALL PARAMETERS
;;;-------------------------------------
readcol,path+'Errors2.txt',F='(A,D,D,D,D,F,F,F,F,F,F)',designation,RA,RA_ERR,DEC,DEC_Err,parallax,parallaxe,pmra,pmrae,pmdec,pmdece
ra_err=ra_err/(3600.*1000.) & dec_err=dec_err/(3600.*1000.)  ;;CONVERT RA AND DEC ERROR FROM DEG TO MAS
distance=double(1000.)/parallax
distanceh=double(1000.)/(parallax+parallaxe)
distancel=double(1000.)/(parallax-parallaxe)
d_err=(ABS(distancel-distanceh)/2.)*0.0  ;;COMPUTE THE AVERAGE DISTANCE ERROR
pmra=pmra/1000. & pmrae=pmrae/1000. & pmdec=pmdec/1000. & pmdece=pmdece/1000.
vrad=ra*0.0 & vrade=ra*0.0

;;;-------------------------------------
;;;INITIAL XYZ AND UVW CALCULATIONS WITH ASSOCIATED UNCERTAINTIES
;;;-------------------------------------
xyz_errors, ra, dec, distance, ra_err, dec_err, d_err, X, Y, Z, EX, EY, EZ, FAST=fast
uvw_errors,ra, dec, pmra*1000. ,pmdec*1000., vrad,distance,ra_err, dec_err,pmrae*1000.,pmdece*1000.,vrade,d_err,U_J, V_J, W_J, E_U_J, E_V_J, E_W_J,/PRECISE

;;;-------------------------------------
;;;BEGIN THE PROPOGATION OF TIME OVER THE GIVEN INTERVAL
;;;-------------------------------------
USEYEAR=5000 ;;THE CHOSEN INTERVAL OF TIME FORWARD AND BACKWARD TO SEARCH FOR THE ETZ ENTRY AND EXIT
index=findgen(USEYEAR)
xline=findgen(USEYEAR)
xline2=-xline
xline=[xline,xline2]
so=sort(xline)
xline=xline[so] ;;A CREATED TIMELINE TO BE USED IN THE CALCULATION

;;;-------------------------------------
;;;OPEN FILE FOR ERROR PROPOGATION
;;;-------------------------------------
openw,lun,path+'Error-Final.txt',/get_lun

;;;-------------------------------------
;;;LOOP OVER ALL STARS, PROPOGATING WITH UVW INFLUENCED CHANGES IN XYZ AND DETERMINE ETZ ENTRY AND EXIT
;;;-------------------------------------

for j=0,n_elements(ra)-1 do begin
        Xn_J=X[j]+(U_J[j]*xline*60*60*24*365.*3.24078e-14) ;;CONVERTING KM/S TO PC OVER 1 YEAR
        Yn_J=Y[j]+(V_J[j]*xline*60*60*24*365.*3.24078e-14) ;;CONVERTING KM/S TO PC OVER 1 YEAR
        Zn_J=Z[j]+(W_J[j]*xline*60*60*24*365.*3.24078e-14) ;;CONVERTING KM/S TO PC OVER 1 YEAR
        xyz_errors_inverse, Xn_J, Yn_J, Zn_J, EX, EY, EZ, ra_n_J, dec_n_J, dist_n_J, era_n_J, edec_n_J, edist_n_J, FAST=fast
        euler,RA_N_J,DEC_N_J,elong,elat,3
        
    
        si_etz=where(elat ge -0.264 and elat lt 0.264,c) ;;ETZ ZONE DEFINED AS +/-0.264 DEGREES IN THE ECLIPTIC
            if c gt 0 then begin
                si_etz2=where(elat ge -0.132 and elat lt 0.132,c2) ;;rETZ ZONE DEFINED AS +/-0.132 DEGREES IN THE ECLIPTIC
                if c2 gt 0 then begin
                    etzboth='yes' & etz2S=min(si_etz2) & etz2E=max(si_etz2)
                endif
                if c2 lt 0 then begin
                    etzboth='no' & etz2S=0 & etz2E=0
                endif
                if c2 eq 0 then begin
                    if si_etz2 ge 0 then begin
                        etzboth='yes' & etz2S=min(si_etz2) & etz2E=max(si_etz2)
                    endif
                endif
                if c2 eq 0 then begin
                    if si_etz2 lt 0 then begin
                        etzboth='no' & etz2S=0 & etz2E=0
                    endif
                endif
           
;;;-------------------------------------
;;;PLOT THE TWO ZONES FOR VISUAL INSPECTION
;;;-------------------------------------
            device,filename=path+'Error3/Obj_'+string(designation[j])+'.eps',/encapsulated,color=1 ;;HAVE TO MAKE THIS DIRECTORY TO STORE PLOTS
            plot,elong,elat,xstyle=1,ystyle=1,xtitle='Ecliptic Long',ytitle='Ecliptic Latitude',title=Designation[j],yrange=[-0.35,0.35],xrange=[min(elong)-.01,max(elong)+.01]
            oplot,[0,500],[-0.264,-0.264]
            oplot,[0,500],[0.264,0.264]
            oplot,[0,500],[-0.132,-0.132],linestyle=2
            oplot,[0,500],[0.132,0.132],linestyle=2
            xyouts,elong[max(si_etz)],elat[max(si_etz)],xline[max(si_etz)]
            xyouts,elong[min(si_etz)],elat[min(si_etz)],xline[min(si_etz)]
        endif

;;;-------------------------------------
;;;SAVE THE INITIAL VALUES FOR VALIDATION IN LATER STEPS
;;;-------------------------------------
        OFFICIAL264s=xline[min(si_etz)]
        OFFICIAL264e=xline[max(si_etz)]
        OFFICIAL132s=xline[etz2S]
        OFFICIAL132e=xline[etz2E]

;;;-------------------------------------
;;;FOR RANGE OF ETZ ZONES WE WILL LOOP OVER UNCERTAINTIES IN POSITION OVER TIME AND PLACE VALUES IN NEW ARRAYS
;;;-------------------------------------
;;NUMS IS THE DECIDED NUMBER OF LOOPS TO RUN OVER FOR UNCERTAINTY CHECK
        start262_B=fltarr(nums)
        end262_B=fltarr(nums)
        start132_B=fltarr(nums)
        end132_B=fltarr(nums)
            for k=0,nums-1 do begin
                U_add=(randomn(seed,1)*E_U_J[j]) & U_add=U_add[0] & X_add=(randomn(seed,1)*EX[j]) & X_add=X_add[0] & Xn_J=(X[j]+X_add)+((U_J[j]+U_add)*xline*60*60*24*365.*3.24078e-14)
                V_add=(randomn(seed,1)*E_V_J[j]) & V_add=V_add[0] & Y_add=(randomn(seed,1)*EY[j]) & Y_add=Y_add[0] & Yn_J=(Y[j]+Y_add)+((V_J[j]+V_add)*xline*60*60*24*365.*3.24078e-14)
                W_add=(randomn(seed,1)*E_W_J[j]) & W_add=W_add[0] & Z_add=(randomn(seed,1)*EZ[j]) & Z_add=Z_add[0] & Zn_J=(Z[j]+Z_add)+((W_J[j]+W_add)*xline*60*60*24*365.*3.24078e-14)
                xyz_errors_inverse, Xn_J, Yn_J, Zn_J, EX, EY, EZ, ra_n_J, dec_n_J, dist_n_J, era_n_J, edec_n_J, edist_n_J, FAST=fast
                euler,RA_N_J,DEC_N_J,elong,elat,3

                si_etz=where(elat ge -0.264 and elat lt 0.264,c)
                    if c gt 0 then begin
                        si_etz2=where(elat ge -0.132 and elat lt 0.132,c2)
                    if c2 gt 0 then begin
                        etzboth='yes' & etz2S=min(si_etz2) & etz2E=max(si_etz2)
                    endif
                    if c2 lt 0 then begin
                        etzboth='no' & etz2S=0 & etz2E=0
                    endif
                    if c2 eq 0 then begin
                    if si_etz2 ge 0 then begin
                        etzboth='yes' & etz2S=min(si_etz2) & etz2E=max(si_etz2)
                    endif
                    endif
                    if c2 eq 0 then begin
                    if si_etz2 lt 0 then begin
                        etzboth='no' & etz2S=0 & etz2E=0
                    endif
                    endif

            oplot,elong,elat,color=blue ;;OPLOT THE RANGE OF VALUES FROM PROPOGATED UNCERTAINTIES
            start262_B[k]=xline[min(si_etz)] & end262_B[k]=xline[max(si_etz)] & start132_B[k]=xline[etz2S] & end132_B[k]=xline[etz2E] ;;PRESERVE NEW VALUES FOR STDEV CHECK
            endif
    endfor
                        device,/close

;;;-------------------------------------
;;;PRINT VALUES TO SCREEN AND TO FILE
;;;-------------------------------------

print,j,designation[j],ra[j],dec[j],distance[j],official264s,official264e,official132s,official132e,median(start262_B),median(end262_B),median(start132_B),median(end132_B),stdev(start262_B),stdev(end262_B),stdev(start132_B),stdev(end132_B),nums,F='(I10,A40,F20.5,F20.5,F20.5,I10,I10,I10,I10,I10,I10,I10,I10,F10.2,F10.2,F10.2,F10.2,I10)

printf,lun,j,designation[j],ra[j],dec[j],distance[j],official264s,official264e,official132s,official132e,median(start262_B),median(end262_B),median(start132_B),median(end132_B),stdev(start262_B),stdev(end262_B),stdev(start132_B),stdev(end132_B),nums,F='(I10,A40,F20.5,F20.5,F20.5,I10,I10,I10,I10,I10,I10,I10,I10,F10.2,F10.2,F10.2,F10.2,I10)

endfor
free_lun,lun
;;;-----------------------------------------------------------------
;;;-----------------------------------------------------------------
;;;-----------------------------------------------------------------


end


