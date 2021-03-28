pro transits,USEYEAR

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
;;;READING IN ALL PARAMETERS USING THE GCNS AS INPUT
;;;-------------------------------------

data=mrdfits(path+'GCNS_cat.fits',1)
VRAD=double(data.RADIAL_VELOCITY)
parallax=double(data.parallax)
parallaxe=double(data.parallax_error)
pmra=double(data.pmra)
pmrae=double(data.pmra_error)
pmdec=double(data.pmdec)
pmdece=double(data.pmdec_error)
RA=double(data.RA)
DEC=double(data.DEC)
designation=data.source_ID
RA_Err=double(data.RA_ERROR)
DEC_Err=double(data.DEC_ERROR)
parallax_error=data.PARALLAX_ERROR

;;;-------------------------------------
;;;WE WILL IGNORE RV AS MOST DON'T HAVE IT AND WE ARE LOOKING FOR PROJECTION POSITION ON SKY.  CALCULATE DISTANCE, DISTANCE ERROR AND SET PM VALUES
;;;-------------------------------------
VRAD=parallax*0.0 & VRADE=RV*0.0
distance=double(1000.)/parallax
distanceh=double(1000.)/(parallax+parallax_error)
distancel=double(1000.)/(parallax-parallax_error)
d_err=ABS(distancel-distanceh)/2.
pmra=pmra/1000. & pmrae=pmrae/1000. & pmdec=pmdec/1000. & pmdece=pmdece/1000.

;;;-------------------------------------
;;;SET THE EXOPLANET HOST STAR LIST
;;;-------------------------------------
readcol,path+'Exoplanets.txt',F='(A,F,F)',name_E,ra_Exo,dec_Exo
srcor,ra_Exo,dec_Exo,ra,dec,0.6,ind1,ind2,option=1
VRAD_E=VRAD[ind2]
VRADE_E=VRADE[ind2]
Parallax_E=parallax[ind2]
parallaxe_E=parallaxe[ind2]
distance_E=distance[ind2]
d_err_E=d_err[ind2]
PMRA_E=pmra[ind2]
PMRAe_E=pmrae[ind2]
PMDEC_E=pmdec[ind2]
PMDECe_E=pmdece[ind2]
RA_E=RA[ind2]
DEC_E=DEC[ind2]
RA_ERR_E=RA_ERR[ind2]
DEC_ERR_E=DEC_ERR[ind2]
name_E_n=name_E[ind1]

;;;-------------------------------------
;;;-------------------------------------
;;;-------------------------------------
;;;BEGIN WITH EXOPLANET HOST STAR LIST
;;;-------------------------------------
;;;-------------------------------------
;;;-------------------------------------


;;;-------------------------------------
;;;BEGIN THE PROPOGATION OF TIME OVER THE GIVEN USEYEAR INTERVAL
;;;-------------------------------------
;;USEYEAR IS THE CHOSEN INTERVAL OF TIME FORWARD AND BACKWARD TO SEARCH FOR THE ETZ ENTRY AND EXIT
index=findgen(USEYEAR)
xline=findgen(USEYEAR)
xline2=-xline
xline=[xline,xline2]
so=sort(xline)
xline=xline[so]  ;;A CREATED TIMELINE TO BE USED IN THE CALCULATION

;;;-------------------------------------
;;;SET INITIAL VALUES AND CALCULATE XYZUVW AT CURRENT EPOCH
;;;-------------------------------------
etzboth='no' & etz2S=0 & etz2E=0
xyz_errors, ra_E, dec_E, distance_E, ra_err_E, dec_err_E, d_err_E, X, Y, Z, EX, EY, EZ, FAST=fast
uvw_errors,ra_E, dec_E, pmra_E*1000. ,pmdec_E*1000., vrad_E,distance_E,ra_err_E, dec_err_E,pmrae_E*1000.,pmdece_E*1000.,vrade_E,d_err_E,U_J, V_J, W_J, E_U_J, E_V_J, E_W_J,/PRECISE

;;;-------------------------------------
;;;OPEN FILE FOR DATA RETENTION
;;;-------------------------------------
openw,lun,path+'NEW-TEST-Exoplanet-'+string(USEYEAR)+'.txt',/get_lun


;;;-------------------------------------
;;;LOOP OVER THE GIVEN TIME INTERVAL AND DETERMINE ETZ ENTRY AND EXIT
;;;-------------------------------------
for j=0,n_elements(ra_E)-1 do begin
        Xn_J=X+((U_J[j])*xline*60*60*24*365.*3.24078e-14)
        Yn_J=Y+(V_J[j]*xline*60*60*24*365.*3.24078e-14)
        Zn_J=Z+(W_J[j]*xline*60*60*24*365.*3.24078e-14)
        xyz_errors_inverse, Xn_J, Yn_J, Zn_J, EX_J, EY_J, EZ_J, ra_n_J, dec_n_J, dist_n_J, era_n_J, edec_n_J, edist_n_J, FAST=fast
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
;;;PRINT VALUES TO SCREEN AND TO FILE
;;;-------------------------------------

        printf,lun,F='(I10,A40,F20.5,F20.5,F20.5,I10,I10,I10,A5,I10,I10,I10)',j,name_E_n[j],ra_E[j],dec_E[j],distance_E[j],n_elements(si_etz),xline[min(si_etz)],xline[max(si_etz)],etzboth,n_elements(si_etz2),xline[etz2S],xline[etz2E]
           
        print,F='(I10,A40,F20.5,F20.5,F20.5,I10,I10,I10,A5,I10,I10,I10)',j,name_E_n[j],ra_E[j],dec_E[j],distance_E[j],n_elements(si_etz),xline[min(si_etz)],xline[max(si_etz)],etzboth,n_elements(si_etz2),xline[etz2S],xline[etz2E]

;;;-------------------------------------
;;;PLOT THE TWO ZONES FOR VISUAL INSPECTION
;;;-------------------------------------
       device,filename=path+'TEST-Figures-10000-Exoplanets/Obj_'+string(name_E_n[j])+'.eps',/encapsulated,color=1 ;;HAVE TO MAKE THIS DIRECTORY TO STORE PLOTS
           plot,elong,elat,yrange=[-0.35,0.35],xrange=[min(elong)-.01,max(elong)+.01],xstyle=1,ystyle=1,xtitle='Ecliptic Long',ytitle='Ecliptic Latitude',title=name_E_n[j]
           oplot,[0,500],[-0.264,-0.264]
           oplot,[0,500],[0.264,0.264]
           oplot,[0,500],[-0.132,-0.132],linestyle=2
           oplot,[0,500],[0.132,0.132],linestyle=2
           xyouts,elong[max(si_etz)],elat[max(si_etz)],xline[max(si_etz)]
           xyouts,elong[min(si_etz)],elat[min(si_etz)],xline[min(si_etz)]
           device,/close
           
            endif
endfor
free_lun,lun
        
;;;-----------------------------------------------------------------
;;;-----------------------------------------------------------------
;;;-----------------------------------------------------------------























;;;-------------------------------------
;;;-------------------------------------
;;;-------------------------------------
;;;USING THE GCNS FILE AS INPUT
;;;-------------------------------------
;;;-------------------------------------
;;;-------------------------------------


;;;-------------------------------------
;;;BEGIN THE PROPOGATION OF TIME OVER THE GIVEN USEYEAR INTERVAL
;;;-------------------------------------
;; USEYEAR IS THE CHOSEN INTERVAL OF TIME FORWARD AND BACKWARD TO SEARCH FOR THE ETZ ENTRY AND EXIT
index=findgen(USEYEAR)
xline=findgen(USEYEAR)
xline2=-xline
xline=[xline,xline2]
so=sort(xline)
xline=xline[so]  ;;A CREATED TIMELINE TO BE USED IN THE CALCULATION

;;;-------------------------------------
;;;SET INITIAL VALUES AND CALCULATE XYZUVW AT CURRENT EPOCH
;;;-------------------------------------
etzboth='no' & etz2S=0 & etz2E=0
xyz_errors, ra, dec, distance, ra_err, dec_err, d_err, X, Y, Z, EX, EY, EZ, FAST=fast
uvw_errors,ra, dec, pmra*1000. ,pmdec*1000., vrad,distance,ra_err, dec_err,pmrae*1000.,pmdece*1000.,vrade,d_err,U_J, V_J, W_J, E_U_J, E_V_J, E_W_J,/PRECISE

;;;-------------------------------------
;;;OPEN FILE FOR DATA RETENTION
;;;-------------------------------------
openw,lun,path+'NEW-TEST-100pc-'+string(USEYEAR)+'.txt',/get_lun


;;;-------------------------------------
;;;LOOP OVER THE GIVEN TIME INTERVAL AND DETERMINE ETZ ENTRY AND EXIT
;;;-------------------------------------
for j=0,n_elements(ra)-1 do begin
        Xn_J=X+((U_J[j])*xline*60*60*24*365.*3.24078e-14)
        Yn_J=Y+(V_J[j]*xline*60*60*24*365.*3.24078e-14)
        Zn_J=Z+(W_J[j]*xline*60*60*24*365.*3.24078e-14)
        xyz_errors_inverse, Xn_J, Yn_J, Zn_J, EX_J, EY_J, EZ_J, ra_n_J, dec_n_J, dist_n_J, era_n_J, edec_n_J, edist_n_J, FAST=fast
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
;;;PRINT VALUES TO SCREEN AND TO FILE
;;;-------------------------------------
        printf,lun,F='(I10,A40,F20.5,F20.5,F20.5,I10,I10,I10,A5,I10,I10,I10)',j,designation[j],ra[j],dec[j],distance[j],n_elements(si_etz),xline[min(si_etz)],xline[max(si_etz)],etzboth,n_elements(si_etz2),xline[etz2S],xline[etz2E]
        print,F='(I10,A40,F20.5,F20.5,F20.5,I10,I10,I10,A5,I10,I10,I10)',j,designation[j],ra[j],dec[j],distance[j],n_elements(si_etz),xline[min(si_etz)],xline[max(si_etz)],etzboth,n_elements(si_etz2),xline[etz2S],xline[etz2E]

;;;-------------------------------------
;;;PLOT THE TWO ZONES FOR VISUAL INSPECTION
;;;-------------------------------------
            device,filename=path+'TEST-Figures-10000-100pc/Obj_'+string(designation[j])+'.eps',/encapsulated,color=1 ;;HAVE TO MAKE THIS DIRECTORY TO STORE PLOTS
            plot,elong,elat,yrange=[-0.35,0.35],xrange=[min(elong)-.01,max(elong)+.01],xstyle=1,ystyle=1,xtitle='Ecliptic Long',ytitle='Ecliptic Latitude',title=Designation[j]
            oplot,[0,500],[-0.264,-0.264]
            oplot,[0,500],[0.264,0.264]
            oplot,[0,500],[-0.132,-0.132],linestyle=2
            oplot,[0,500],[0.132,0.132],linestyle=2
            xyouts,elong[max(si_etz)],elat[max(si_etz)],xline[max(si_etz)]
            xyouts,elong[min(si_etz)],elat[min(si_etz)],xline[min(si_etz)]
            device,/close
        endif
endfor
free_lun,lun
        
;;;-----------------------------------------------------------------
;;;-----------------------------------------------------------------
;;;-----------------------------------------------------------------



end


