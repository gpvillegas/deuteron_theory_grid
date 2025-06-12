    	subroutine f_NN_on(s,t,f_NN_onre,f_NN_onim,ins)
*         t         = -q**2
        sigma_tot = sigma_pn_tot(s) / 0.389385   ! Gev^{-2}
        an        = a_pn(s)
        expo      = expo_pn(s,t)

        f_NN_onre = sigma_tot*an*expo
        f_NN_onim = sigma_tot*expo
        return
        end
       
        subroutine f_NN_off(s,t,f_NN_offre,f_NN_offim,ins)
*         t          = -q**2
        sigma_tot  = sigma_pn_tot(s) / 0.389385   ! Gev^{-2}
        an         = a_pn(s)
        expo       = expo_pn(s,t)

        f_NN_offre = sigma_tot*an*expo
        f_NN_offim = sigma_tot*expo
*        print *,an,s
        return
        end



       function sigma_pn_tot(s)
***************************************************************************
*      Tabra
*      pn total cross section fit in [mb]  (sigma_pn)
***************************************************************************
       common/par/pi,pm,pmp,pmn,tm,eb ! pm -> mass of the proton
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
        stot =  cubic_spline(p,x1,x2,x3,y1,y2,y3)
       else
        stot = 47.3 + 0.513*alog(p)**2 - 4.27*alog(p)
       endif
*       sigma_pn_tot = 40.0
       sigma_pn_tot = stot
       return
       end

        function expo_pn(s,t)   
        common/par/pi,pm,pmp,pmn,tm,eb
        expo_pn = 0.0
        pcm2 = (s-4.*pm**2)/4. 
        sn_thkes2 = -t/4./pcm2
        if(sn_thkes2.gt.1.0)then
*        write(6,*) "ERROR: function expo_pp: sn_thkes2 is gt 1 "
        expo_pn = 0.0
        return
        endif
 	expo_pn =  exp(bn(s)*t/2.0)
   	return
	end

       function bn(s)
***************************************************************************
*      tabra
*      pn Slope parameter (Bpn)
***************************************************************************
       common/par/pi,pm,pmp,pmn,tm,eb ! pm -> mass of the proton
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
        bno =  cubic_spline(p,x1,x2,x3,y1,y2,y3)
       elseif(p.lt.1.74)then  
        bno = -0.41579 + 1.8333*p + 2.6484*p**2  -1.0031*p**3              
       elseif(p.le.3.0)then
        x1 =  1.740000   
        y1 =  5.508093 
        x2 =  2.4
        y2 =  7.202996 
        x3 =  3.0
        y3 =  7.903550 
        bno =  cubic_spline(p,x1,x2,x3,y1,y2,y3)
       else  
        bno  = 8.22  + 1.10*alog(p/4.)
       endif  
       bn = bno
       return 
       end

       function a_pn(s)
****************************************************************************
*      tabra
*      pn Re/Im amplidtudes ratio  parameter (alpha_pn)
****************************************************************************
       common/par/pi,pm,pmp,pmn,tm,eb  ! pm -> mass of the proton
       dimension B(10)
       ep = (s - 2.0*pm**2)/(2.0*pm)
       p = sqrt(ep**2-pm**2)
       a_pn = 0.0
       if(p.lt.1.4) then
        call Bystricky_dat(p,0.0,1,B)!using SAID's analysis
        a_pn = ( B(1)+B(2))/( B(6)+B(7))
       return
       elseif(p.ge.1.4.and.p.lt.3.) then
        x1 =  1.4000
        y1 =  -0.2845016
        x2 =  1.9
        y2 =  -0.43915290
        x3 =  3.0
        y3 =  -0.4938042
        a_pn =  cubic_spline(p,x1,x2,x3,y1,y2,y3)
       return
       elseif(p.ge.3..and.p.lt.101.) then
        a_pn=-0.56207 +0.24223E-01*p -0.50362E-03*p**2 
     &       +0.48408E-05*p**3 -0.17331E-07*p**4
       endif
       return
       end

        

       function cubic_spline(p,x1,x2,x3,y1,y2,y3)
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
       if(p.le.x2) cubic_spline = C1
       if(p.gt.x2) cubic_spline = C2
       return 
       end
