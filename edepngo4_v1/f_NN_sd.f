    	subroutine f_NN_sd(is1,is2,is1p,is2p,
     &                           s,t,f_NN_onre,f_NN_onim,ins,ini)
*******************************************************************************
* parameterization of forward NN scattering amplitude based on 
* the combination of SAID and diffractive amplitudes   
* Up to plab  (1.3) 3.68 it works as SAID after it is diffrective
*******************************************************************************
*                    spin polarized amplitudes in GeV-2
*             Normalization is: Im(f)=sigma_tot (so, 4*pi/pcm seats in f)
*input
* s                  --  (p_1^mu+p_2^mu)^2
* t                  --  (p_1^mu-p_3^mu)^2 
*
* is1,is2,is1p,is2p  -- initial and final spins of scattered nucleons 1-up -1 down
*
* ins                -- NN flag 21=pn->pn, 22 = pp->pp  
*output 
*f_NN_onre           -- real part of amplitude
*f_NN_onim           -- imaginary part of amplitude
*          when plab < 1.4 said amplitudes are involved
*          when plab > 1.4 old diffraction model is involved
*ini                 -- initializes Saclay's(Bystricky's) amplitudes if ini=0 
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        dimension B(10),BB(10) 
        common/par/pi,pm,pmp,pmn,tm,eb  
        common/Bystr/ ibys
        ibys = 1
***************************************************************************
*ibys=1 one needs to call for Bystricky_dat(0.01,0.,'negative number',BB) *
*       to read data from the files  B_ampl_pp.dat, B_ampl_pn.dat         *   
*       and B_ampl_pp_nocoul.dat previously obtained from the SAID code   *
*                                                                         *
*ibys=0 one needs to activate CALL NNSOLA(tlab,acm,IR,TR,TI,H,TTL)        *
*       subroutine (not included with this file) in subroutine Bystricky  * 
*       to obntain amplitudes  directly from the SAID's NNSOLA code       *
***************************************************************************

        if (ini.eq.0)then
        call Bystricky_dat(0.01,0.,-1,BB)
        pi  = acos(-1.0)       
        pm  = 0.938279         
        return
        endif        
        f_NN_onre = 0.0
        f_NN_onim = 0.0

        ep = (s - 2.0*pm**2)/(2.0*pm)
        if(ep.le.pm)return

        p = sqrt(ep**2-pm**2)
        pcm = sqrt(s-4*pm**2)/2.0 ! in GeV
*        argum = 1.0+ t/(2.0*pcm**2) ! (2.*t+s-4.*pm**2)/(s-4.*pm**2)      
*        if(argum.gt.1..and.argum.lt.1.01) argum =  0.999
*        if(argum.lt.-1.0.and.argum.gt.-1.01) argum = -0.999
*        if(argum**2.gt.1.0)return
*        acm_rad = acos(argum)
*        acm  = (acm_rad/pi)*180.
*        print *,ins,p,ini,acm,argum

        sn_thkes2 = -t/4.0/pcm**2
        if(sn_thkes2.gt.1.0)return
        acm_rad = 2.0*asin(sqrt(sn_thkes2))
        acm = acm_rad*180.0/pi

***********************************************************************
      IF(abs(ins).eq.21.or.abs(ins).gt.100)then        !pn scattering
***********************************************************************
*        if (p.le.3.68) then
        if (p.le.1.4) then
          fMre = 0.0; fMim = 0.0
   
          call f_M(s,t,is1,is2,is1p,is2p,
     &                      ins,fMre,fMim,B)! Look in file Bystricky.f
*************************************************************************************
*           output fMre and fMim are in sqrt(mbn)
*           Normalization is dsg=(Tr(fMre**2)+Tr(fMim**2))/4
*                 (dsg is dsigma/dOmega)
*          Now changing to normalization Im(f)=sigma_tot
*                     (Im(f) = Tr(fMim)/4)
*************************************************************************************        
          fMre = fMre/sqrt(0.389385)  ! sqrt(mbn) to GeV-1
          fMim = fMim/sqrt(0.389385)
          
          f_NN_onre = 4.*pi/pcm*fMre  ! in GeV-2
          f_NN_onim = 4.*pi/pcm*fMim  ! in GeV-2
*          return
*        if(is1.eq.is1p.and.is2.eq.is2p)print *,p,f_NN_onre/f_NN_onim,acm
        else 
          f_NN_onre =0.0;  f_NN_onim = 0.0
          if(is1.eq.is1p.and.is2.eq.is2p) then   
            sigma_tot = sigma_pn_tot_d(s) / 0.389385   ! mbn to Gev^{-2}
            an   = a_pn_d(s)
            expo = expo_pn_d(s,t)
            f_NN_onre = sigma_tot*an*expo ! in GeV-2
            f_NN_onim = sigma_tot*expo    ! in GeV-2
          endif
! If not diagonal then automaticly f_NN_onre =0.0;  f_NN_onim = 0.0
*          return

        endif

********************************************************************
      ELSEIF(ins.eq.22.or.ins.eq.-22)THEN                   ! pp scattering 
********************************************************************
        if (p.le.3.68) then
          fMre = 0.0; fMim = 0.0
          call f_M(s,t,is1,is2,is1p,is2p,
     &                      ins,fMre,fMim,B)! Look in file Bystricky.f

*             output fMre and fMim are in sqrt(mbn)
*             Normalization is dsg=(Tr(fMre**2)+Tr(fMim**2))/4
*                     (dsg is dsigma/dOmega)
*             Now changing to normalization Im(f)=sigma_tot
*                     (Im(f) = Tr(fMim)/4)
        
          fMre = fMre/sqrt(0.389385)  ! sqrt(mbn) to GeV-1
          fMim = fMim/sqrt(0.389385)

          f_NN_onre = 4.*pi/pcm*fMre  ! in GeV-2
          f_NN_onim = 4.*pi/pcm*fMim  ! in GeV-2
          return
        else 
          f_NN_onre =0.0;  f_NN_onim = 0.0
          if(is1.eq.is1p.and.is2.eq.is2p) then   
            sigma_tot = sigma_pp_tot(s) / 0.389385   ! mbn to Gev^{-2}
            ap   = a_pp(s)
            expo = expo_pp(s,t)
            f_NN_onre = sigma_tot*ap*expo ! in GeV-2
            f_NN_onim = sigma_tot*expo    ! in GeV-2
          endif
          return

        endif
      
      ENDIF

        return
        end


*************************************************************************
* exponential for pp scattering amplitude
*************************************************************************
        function expo_pp(s,t)
        common/par/pi,pm,pmp,pmn,tm,eb  
        pcm2 = (s-4.*pm**2)/4. 
        sn_thkes2 = -t/4./pcm2
        if(sn_thkes2.gt.1.0)then
        expo_pp = 0.0
        return
        endif
  	expo_pp =  exp(bp(s)*t/2.0)
   	return
	end

**************************************************************************
* exponential for pn scattering amplitude
* _d means that this is the duplicate of fpn.f
**************************************************************************
        function expo_pn_d(s,t)
        common/par/pi,pm,pmp,pmn,tm,eb  
        pcm2 = (s-4.*pm**2)/4. 
        sn_thkes2 = -t/4./pcm2
        if(sn_thkes2.gt.1.0)then
*         write(6,*) "ERROR: function expo_pp: sn_thkes2 is gt 1 "
        expo_pn_d = 0.0
        return
        endif
  	expo_pn_d =  exp(bn_d(s)*t/2.0)
   	return
	end

******************************************************************
* Parameterizations for pp and pn scattering
*
* 18-Dec-03
* FIU, Miami
******************************************************************


       function sigma_pp_tot(s)
******************************************************************
*      pp total cross section fit in [mb] (sigma_pp)
******************************************************************
       common/par/pi,pm,pmp,pmn,tm,eb  
       dimension B(10)
       ep = (s - 2.0*pm**2)/(2.0*pm)
       p = sqrt(ep**2-pm**2)
       pcm = pm*p/sqrt(s)              ! in GeV
       pcm = pcm/sqrt(0.389385)        ! in 1/sqrt(mb)
       if(p.lt.3.68) then
        call Bystricky_dat(p,0.0,0,B)!using SAID's analysis
        stot = 4.*pi/pcm *  ( B(6)+B(7) )/2.
       elseif ( p.le.5.87) then
        x1 =  3.68
        y1 =  41.39135
        x2 =  3.9
        y2 =  41.19135
        x3 =  5.870000                       
        y3 =  40.42657 
        stot =  cubic_spline_d(p,x1,x2,x3,y1,y2,y3)
       elseif ( p.le.12.) then
        x1 =  5.87
        y1 =  40.42657
        x2 =  9.
        y2 =  40.19135
        x3 =  12.00000                      
        y3 =  40.01630
        stot =  cubic_spline_d(p,x1,x2,x3,y1,y2,y3)
       else
        stot = 48.0 + 0.522*alog(p)**2 - 4.51*alog(p)
       endif
       sigma_pp_tot = stot
       return
       end



       function sigma_pn_tot_d(s)
***************************************************************************
*      pn total cross section fit in [mb]  (sigma_pn)
* _d means that this is the duplicate of fpn.f
***************************************************************************
       common/par/pi,pm,pmp,pmn,tm,eb  
       dimension B(10)
       ep = (s - 2.0*pm**2)/(2.0*pm)
       p = sqrt(ep**2-pm**2)
       pcm = pm*p/sqrt(s)             ! in GeV
       pcm = pcm/sqrt(0.389385)       ! in 1/sqrt(mb)
       if(p.lt.1.4)then
        call Bystricky_dat(p,0.0,1,B)!using SAID's analysis
        stot = 4.*pi/pcm *  ( B(6)+B(7) )/2.
       elseif(p.le.3.1) then
        x1 =  1.4
        y1 =  38.42728
        x2 =  2.07800
        y2 =  41.905
        x3 =  3.1                     
        y3 =  43.13159
        stot =  cubic_spline_d(p,x1,x2,x3,y1,y2,y3)
       else
        stot = 47.3 + 0.513*alog(p)**2 - 4.27*alog(p)
       endif
       sigma_pn_tot_d = stot
       return
       end



       function a_pp(s)
****************************************************************************
*      pp Re/Im amplidtudes ratio  parameter (alpha_pp)
****************************************************************************
       common/par/pi,pm,pmp,pmn,tm,eb  
       dimension B(10)
       ep = (s - 2.0*pm**2)/(2.0*pm)
       p = sqrt(ep**2-pm**2)
       a_pp = 0.0
       if(p.lt.3.68) then
        call Bystricky_dat(p,0.0,0,B)!using SAID's analysis
        a_pp = ( B(1)+B(2) )/( B(6)+B(7) )
       return
       elseif(p.ge.3.68.and.p.lt.5.) then
        x1 =  3.68
        y1 =  -0.3316863
        x2 =  4.34
        y2 =  -0.352568200
        x3 =  5.0
        y3 =  -0.3734501
        a_pp =  cubic_spline_d(p,x1,x2,x3,y1,y2,y3)
       return
       elseif(p.ge.5..and.p.lt.101.) then
        a_pp=-0.43377 +0.13158E-01*p -0.22820E-03*p**2 
     &       +0.19088E-05*p**3 -0.59601E-08*p**4
       return
       endif
       end



       function a_pn_d(s)
****************************************************************************
*      pn Re/Im amplidtudes ratio  parameter (alpha_pn)
* _d means that this is the duplicate of fpn.f
****************************************************************************
       common/par/pi,pm,pmp,pmn,tm,eb  
       dimension B(10)
       ep = (s - 2.0*pm**2)/(2.0*pm)
       p = sqrt(ep**2-pm**2)
       a_pn_d = 0.0
       if(p.lt.1.4) then
        call Bystricky_dat(p,0.0,1,B)!using SAID's analysis
       a_pn_d = ( B(1)+B(2))/( B(6)+B(7))
       return
       elseif(p.ge.1.4.and.p.lt.3.) then
        x1 =  1.4000
        y1 =  -0.2845016
        x2 =  1.9
        y2 =  -0.43915290
        x3 =  3.0
        y3 =  -0.4938042
        a_pn_d =  cubic_spline_d(p,x1,x2,x3,y1,y2,y3)
       return
       elseif(p.ge.3..and.p.lt.101.) then
        a_pn_d=-0.56207 +0.24223E-01*p -0.50362E-03*p**2 
     &       +0.48408E-05*p**3 -0.17331E-07*p**4
       endif
       return
       end



       function bp(s)
***************************************************************************
*      pp Slope parameter (Bpp)
***************************************************************************
       common/par/pi,pm,pmp,pmn,tm,eb  
       ep = (s - 2.0*pm**2)/(2.0*pm)
       p = sqrt(ep**2-pm**2)
       if(p.lt.0.926)then
        x1=  0.       
        y1=  0.              
        x2=  0.745
        y2=  0.6
        x3=  0.926  
        y3=  0.39
        bp = cubic_spline_d(p,x1,x2,x3,y1,y2,y3)
       return
       elseif(p.lt.1.2)then
        x1=  0.9300001       
        y1=  0.39           
        x2=  1.15
        y2=  2.0
        x3=  1.2
        y3=  3.177635
        bp = cubic_spline_d(p,x1,x2,x3,y1,y2,y3)
       return
       elseif(p.lt.1.7)then
        bp =  55.861-125.66*p+95.365*p**2-22.695*p**3 
       return
       elseif(p.le.3.0)then
        x1 =  1.7
        y1 =  6.343315
        x2 =  2.4
        y2 =  7.5
        x3 =  3.0
        y3 =  7.903550
        bp = cubic_spline_d(p,x1,x2,x3,y1,y2,y3)
       return
       else
        bp  = 8.22  + 1.10*log(p/4.)
       endif
       return
       end



       function bn_d(s)
***************************************************************************
*      pn Slope parameter (Bpn)
* _d means that this is the duplicate of fpn.f
***************************************************************************
       common/par/pi,pm,pmp,pmn,tm,eb  
       dimension B1(10),B2(10)
       ep = (s - 2.0*pm**2)/(2.0*pm)
       p = sqrt(ep**2-pm**2)
       if(p.lt.0.8)then
        x1=  0.       
        y1=  0.              
        x2=  0.545
        y2=  1.0 
        x3=  0.8000000       
        y3=  2.232239 
        bn_d =  cubic_spline_d(p,x1,x2,x3,y1,y2,y3)
       return  
       elseif(p.lt.1.74)then  
        bn_d = -0.41579 + 1.8333*p + 2.6484*p**2  -1.0031*p**3              
       return  
       elseif(p.le.3.0)then
        x1 =  1.740000   
        y1 =  5.508093 
        x2 =  2.4
        y2 =  7.202996 
        x3 =  3.0
        y3 =  7.903550 
        bn_d =  cubic_spline_d(p,x1,x2,x3,y1,y2,y3)
       return  
       else  
        bn_d  = 8.22  + 1.10*alog(p/4.)
       endif  
       return 
       end

       function cubic_spline_d(p,x1,x2,x3,y1,y2,y3)
****************************************************************************
*      Approximates function with the 3-nodes cubic spline
****************************************************************************
       dyx2 = (y3-y2)/(x3-x2)
       dyx1 = (y2-y1)/(x2-x1) 
       gamma2 = 6.*( dyx2 - dyx1 )
       yyy2 = gamma2/2./(x3-x1)
       yy1  = dyx1 - yyy2*(x2-x1)/6.
       yy2  = dyx2 - yyy2*(x3-x2)/3.
       C1   = y1+yy1*(p-x1)+yyy2*(p-x1)**3/(x2-x1)/6. 
       C2   = y2+yy2*(p-x2)+yyy2*(p-x2)**2/2.
     &                     -yyy2*(p-x2)**3/(x3-x2)/6.
       if(p.le.x2) cubic_spline_d = C1
       if(p.gt.x2) cubic_spline_d = C2
       return 
       end





