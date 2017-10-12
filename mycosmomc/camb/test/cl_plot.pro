 loadct, 39

 readcol,'my_lcdm_w7_lensedCls.dat',l, cl_lens                                                               

 readcol,'../../../../needlets/indata/wmap7/wmap_binned_tt_spectrum_7yr_v4p1.txt',wl,wcl,werr,format='l,x,x,f,f'

 readcol,'my_lcdm_w7_scalCls.dat',ml,mcl                                                                     

 readcol,'my_tgh_w7_scalCls.dat',l,cl_tgh

 window, 0, tit='!17Angular Power Spectrum'
 plot, l, mcl, chars=1.5, xtit='!17l', ytit='!17l(l+1)C!dl!n/2/!pi !7l!17K!u2'
 oplot, l, cl_tgh, col=245, line=2
 oplot, wl, wcl, psym=1;, col=70
 errplot, wl, wcl-werr, wcl+werr;, col=70

end
