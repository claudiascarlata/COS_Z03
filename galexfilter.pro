function galexfilter, wave,flux,filter,weff,stmag=stmag,$
                      fnu=fnu,pivot=pivot,abmag=abmag,z=z
;+
; NAME
;   GALEXFILTER
; PURPOSE:
;   Convolve a spectrum through a GALEX filter to obtain a flux. If the
;   input spectrum is in units of erg/s/cm^2/Angstrom flambda(lambda), then
;   returns flambda(filter) with options for effective wavelength, ST
;   magnitude, fnu(filter), pivot wavelength, and AB magnitude. May also
;   blueshift or redshift the filter. 
;
; Note: GALEX filter curves are very coarse and irregular. Therefore they are
;       interpolated to 1 angstrom resolution wavelength grid.
;
;      Also modify fpath (path to filter curves) to place in personal library
;      location. see below
;
; CALLING SEQUENCE:
;
;   result = GALEXFILTER( wave, flux, [ filter, weff, ... ])  
;
; INPUTS:
;   WAVE - wavelength vector, should overlap the wavelength range of the
;          GALEX filter
;   FLUX - flux vector, same number of elements as W. 
; OPTIONAL INPUT:
;   FILTER - 1 character string giving the name of GALEX filter to be used.
;            N = nuv filter             
;            F = fuv filter  
;   z      - positive(negative) redshift(blueshift) bandpass
;
; OUTPUT:
;   RESULT - scalar giving intensity in same units as FLUX convolved through 
;          specified filter  
; OPTIONAL OUTPUT:
;   WEFF  - scalar giving the effective wavelength of the flux distribution
;          through the specified filter
;   STmag - STmag = -2.5alog10(flambda_filter) - 21.1
;   Fnu   - erg/s/cm^2/Hz through filter 
;   PIVOT - pivot wavelength
;   ABmag - ABmag = -2.5alog10(fnu_filter) - 48.6
;
; FILES USED:
;   GALEXFILTER reads the files galex.fuv.dat and galex.nuv.dat
;              stored in the text table format.
;
; EXAMPLE:
;
;   NUV = UITFILTER(W,F,'N') ;Relative intensity in the near UV
;   FUV = UITFILTER(W,F,'F') ;Relative intensity in the far UV.
;
; REVISION HISTORY:
;   Written    M. Seibert 5/14/2002
;   Revised    M. Seibert 5/17/2002    
;              Based on UITFILTER 
;-
 On_error,2
 if N_params() LT 2 then begin
  print,'Sytax - result = GALEXFILTER( wave, flux, [ filter, weff, stmag=stmag ... ])'
  return, -1
 endif

 ;specify directory of filter curves (galex.nuv.dat anf galexfuv.dat)
 fpath='./' ; assumes same directory 

 name = ['F','N' ]
 peakwl = [ 1475, 2200 ]
 bandps = [  275, 800  ]  ;these are just an estimate

 if N_params() LT 3 then readfilter = 1b else readfilter = 0b

 GETFILTER: if readfilter then begin
    filter = ''
    print,'The following GALEX filters are available'
    print,'Filter    Peak WL (A)   Bandpass (A)'
    for i = 0,1 do print,f ='(A4,6X,I6,8X,I6)',name(i),peakwl(i),bandps(i)
    print,' '        
    read,'Enter 1 character ID of desired filter: ',filter
 endif

 remchar,filter,' '
 if ( filter ne 'F') and ( filter ne 'N' ) then begin
      message,'ERROR - Unknown filter '+strupcase(filter), /CON
      readfilter = 1b
      goto, GETFILTER 
 endif


 case strupcase(filter) of 
 'N':  begin                      ;NUV Filter (N)
       readcol, fpath+'galex.nuv.dat',w1,frel1,format='f,f',skip=3,/silent
       w=findgen(1500)+1500.
       linterp,w1,frel1,w,frel
       end
 'F':  begin                      ;FUV Filter (F)
       readcol, fpath+'galex.fuv.dat',w1,frel1,format='f,f',skip=3,/silent
       w=findgen(1000)+1240.
       linterp,w1,frel1,w,frel
       end
ELSE: begin
      message,'ERROR - Unknown Filter '+strupcase(filter),/CON
      readfilter = 1b
      goto, GETFILTER   
      end
 endcase

IF keyword_set(z) THEN BEGIN 
 IF z Ge 0. THEN w=w*(1+z)
 IF z lt 0. THEN w=w/(1-z)
ENDIF


 linterp, wave, flux, w, f  ;Linearly interpolate onto new wavelength grid
 feff = f*frel              ;Effective flux through filter
 
 flambda_filter=tsum(w,w*feff)/tsum(w,w*frel);Integrate flux*(rel. intensity)
 weff = tsum(w,w*feff)/tsum(w,feff)          ;effective wavelength
 pivot2=tsum(w,w*frel)/tsum(w,frel/w)        ;pivot wavelength squared
 pivot=sqrt(pivot2)
 fnu_filter = flambda_filter*pivot2/2.99e18  ;convert to erg/s/cm2/Hz
 stmag = -2.5*alog10(flambda_filter) -21.1
 abmag = -2.41 - 2.5*alog10(tsum(w,w*feff)/tsum(w,frel/w)) ;more precise
 
 return,flambda_filter ;Integrate flux*(relative intensity)

 end
