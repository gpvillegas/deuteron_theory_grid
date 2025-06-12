*      program test
*      common/parms/pi,pm,dm
*      dimension q(3),p2(3),p1(3),s1(3),s2(3)
*      real, dimension(2,2) :: jN_re,jN_im
*      complex, dimension(0:3,2,2) :: jN,jN_off
*      pi = acos(-1.0)
*      pm = 0.938279
*      dm = 1.875
*      nuc = 0 ! initialization
*      call  j_N_matrix_s_onof(nuc,q2,q,p2,s2,p1,s1,jN,jN_off,e_on,e_off,noff)

      
*      nuc = 1
*      q2 = 1.0
*      q(1)  = 0.0; q(2)  = 0.0;  q(3) = 1.5 

*      p1(1) = 0.2; p1(2) = 0.1; p1(3) = 0.2
*      theta1 = 0.0 !pi/2.0
*      s1(1) = sin(theta1)
*      s1(2) = 0.0
*      s1(3) = cos(theta1) !pi/2.0

*      p2(1) = 0.2; p2(2) = 0.1; p2(3) = 0.2+q(3)
*      theta2 = 0.0 !pi/2.0
*      s2(1) = sin(theta2)
*      s2(2) = 0.0
*      s2(3) = cos(theta2) !pi/2.0

*      p1m = sqrt(p1(1)**2 + p1(2)**2 + p1(3)**2)
*      e_on  = sqrt(pm**2 + p1m**2)
*      e_off = dm - e_on 

*      e_on  = 0.0
*      e_off = 0.0
      
*      mu = 3
*      call  j_N_matrix_s_onof(nuc,q2,q,p2,s2,p1,s1,jN,jN_off,e_on,e_off,noff)

*      print *,nuc,q2,q,p2,s2,p1,s1

*      jN_re = real(jN(mu,1:2,1:2))
*      jN_im = imag(jN(mu,1:2,1:2))
*      print *,jN_re
*      print *,jN_im
*      print *,jN(3,1:2,1:2)
*      print *,jN_off(3,1:2,1:2)
*      end


      subroutine j_N_matrix_s_onof(nuc,q2,q,p2,s2,p1,s1,jN,jN_off,
     &                                                 e_on,e_off,noff)
********************************************************************************************
* Calculates electromagnetic current of nucleon
* for initial nucleon polarization s1 and final s2
* up - means polarization along s1 or s2
* down - polarization opposite to s1 and s2
*  <s2,p2|J^\mu|p1,s1>  = \bar u(p_2,s_2)\Gamma^\mu u(p_1,s_1) - jN_re/im(s2,s1), where
*  1,1 = u2,u1, 1,2=u2,d1,  2,1=d2,u1, 2,2=d2,d1
*
* and polarization axises for initial lepton is defined by s1 
*     and for final lepton s2
* here \Gamma_\mu = \gamma_mu F1 + i/(2M) F2 q^\nu \sigma_{\mu\nu}
*
* nuc - (0) initializes the \gamma and \sigma matrices
* nuc - (1) - proton (-1) neutron
* q2 - Q2
* q  - q(3)  three momentum of the virtual photon
* p1 - p1(3) three momentum of the initial nucleon
* s1 - s1(3) polarization of initial nucleon
* p2 - p2(3) three momentum of the final nucleon 
* s2 - s2(3) polarization of final nucleon
* jN(0:3,2f,2i) - complex out - returns re and im part of the current
*               0:3 are mu components,  2f - final   nucleon's spin 1-up 2-down
*                                       2i - initial nucleon's spin 1-up 2-down
*up and down means parallel and aniparallel to s1 (for initial) and s2 (for final)
*! common/parms/pi,pm,dm should be provided 
*
* this version also calculates off-shell current like
* ubar Gamma gamma_0 u (e_off-e_on)/2pm
*
* ics - defines the case for form-factor parameterization
*  
*
* Misak Sargsian 
* March 2007,FIU, Miami
* March-April 2008,FIU, Miami
* error in Kelly's parameterization is corrected
* June 2010, Miami
*
****************************************************************************************
      implicit none
      complex, dimension(0:3,4,4):: gamma,bgamma
      complex, dimension(0:3,0:3,4,4) :: sigma
      complex, dimension(4):: u1,uc1,u2,uc2,u2bar,ul0,ulx,uly,ulz,u1bar
      complex, dimension(4)::g0u1 
      complex, dimension(0:3,4):: ucgamma,ucmat,u2bar_bgamma,gammau
      real,    dimension(3) :: q,p2,p1
      real,    dimension(3)   :: s2,s1
      complex, dimension(0:3,2,2):: jN,jN_off

      complex c,cm,j0,jx,jy,jz
      real pi,pm,dm
      real p1x,p1y,p1z,e1,p2x,p2y,p2z,e2,qx,qy,qz,q0,q2,e_on,e_off,f_off 
      integer nuc,is1,isp1,is2,isp2

      real norm_test
      common/gammamatices_s/gamma,sigma
      common/parms/pi,pm,dm
      common/formfcase/ics
      integer ics,noff

      jN = cmplx(0.0,0.0)
      jN_off = cmplx(0.0,0.0)
      if(nuc.eq.0)then

*      ics = 1  ! defines the parameterization of nucleon form-factors
               ! 1 - SLAC, 2 - Kelly, 3-Bodek, BB, Arrington
               ! 4 - Kelly with different Gen ( =0)
               ! 5 - Kelly with different Gen (2Gen_Galster)
      call gammamatricess  ! initializes gamma and sigma matrices

      return
      endif

      p1x = p1(1)
      p1y = p1(2)
      p1z = p1(3)
      e1  = sqrt(pm**2 + p1x**2 + p1y**2 + p1z**2)

      p2x = p2(1)
      p2y = p2(2)
      p2z = p2(3)
      e2  = sqrt(pm**2 + p2x**2 + p2y**2 + p2z**2)

      
      qx = q(1)
      qy = q(2)
      qz = q(3)
      q0 = sqrt(qx**2+qy**2+qz**2-q2)



      do is1 = 1,2
                     isp1 =   1
         if(is1.eq.2)isp1 =  -1
      do is2 = 1,2
                     isp2 =  1
         if(is2.eq.2)isp2 = -1
      
*******************************************************
*      Spinor of initial nucleon
*******************************************************
      call  us_s(pm,p1x,p1y,p1z,e1,s1,isp1,u1,uc1)
*      print *,"us",u1
      call vec_mat_prods_inv(gamma,u1,gammau)  
      g0u1 = gammau(0,1:4)
*      write(24,*)u1
*      write(24,*)"*************"
*      write(24,*)g0u1
*      norm_test  = u1bar(1)*u1(1) + u1bar(2)*u1(2) + 
*     &             u1bar(3)*u1(3) + u1bar(4)*u1(4)
*      print *,"norm",norm_test
*      call  u_s(pm,p1x,p1y,p1z,e1,isp1,u1,uc1)
*      print *,u1
*******************************************************
*      Spinor of final nucleon
*******************************************************
      call  us_s(pm,p2x,p2y,p2z,e2,s2,isp2,u2,uc2)
*      print *,"s",s
*      call  u_s(pm,p2x,p2y,p2z,e2,isp2,u2,uc2)
*      a =   matmul(gamma(1,1:4,1:4),gamma(1,1:4,1:4))
*      b =   matmul(gamma(1,1:4,1:4),gamma(1,1:4,1:4))
*      call vec_vec_sprod(uc1,u1,c)
********************************************************
*     Calculating u2bar
********************************************************
      call vec_mat_prods(uc2,gamma,ucgamma)  
      u2bar = ucgamma(0,1:4)
*      print *,u2bar
******************************************************************
*      Calculating nucleon vertex
******************************************************************
      call nucleon_vertex(nuc,q0,qx,qy,qz,bgamma)

*      print *, "bgamma",bgamma(0,1:4,1:4)
*****************************************************************
*      Calculating \bar u2 \Gamma part
*****************************************************************
      call vec_mat_prods(u2bar,bgamma,u2bar_bgamma)  
      
      ul0 = u2bar_bgamma(0,1:4)
      ulx = u2bar_bgamma(1,1:4)
      uly = u2bar_bgamma(2,1:4)
      ulz = u2bar_bgamma(3,1:4)

*      call vec_vec_sprod(ul0,u2,c)
*      print *, "c",c,2.*e2
*      call vec_vec_sprod(uc2,u1,c)
*      print *, "c",c

****************************************************************
*      Completing the calculation of current components
****************************************************************
      call vec_vec_sprods(ul0,u1,j0)
      call vec_vec_sprods(ulx,u1,jx)
      call vec_vec_sprods(uly,u1,jy)
      call vec_vec_sprods(ulz,u1,jz)

      jN(0,is2,is1) = j0
      jN(1,is2,is1) = jx
      jN(2,is2,is1) = jy
      jN(3,is2,is1) = jz


****************************************************************
*      Completing the calculation of off shell current components
****************************************************************
      call vec_vec_sprods(ul0,g0u1,j0)
      call vec_vec_sprods(ulx,g0u1,jx)
      call vec_vec_sprods(uly,g0u1,jy)
      call vec_vec_sprods(ulz,g0u1,jz)

      f_off = (e_off-e_on)/(2.0*pm)

      jN_off(0,is2,is1) = j0*f_off
      jN_off(1,is2,is1) = jx*f_off
      jN_off(2,is2,is1) = jy*f_off
      jN_off(3,is2,is1) = jz*f_off

*      print *, "j0",j0
      enddo
      enddo

      return
      end

      subroutine vec_vec_sprods(a,b,c)
      complex, dimension(4) :: a, b
      complex c
      c = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)+a(4)*b(4)
      return
      end

      subroutine vec_mat_prods(a,b,pr)
      complex, dimension(4) :: a
      complex, dimension(0:3,4,4) :: b
      complex, dimension(0:3,4):: pr
      pr = cmplx(0.0,0.0)
      do mu = 0,3
      do j = 1,4
      do i = 1,4
      pr(mu,j)  = pr(mu,j) + a(i)*b(mu,i,j)
      enddo
      enddo
      enddo
*      print *,"*",a
*      print *,"***", pr(0,1:4)
      return
      end

      subroutine vec_mat_prods_inv(bm,a,pr)
      complex, dimension(4) :: a
      complex, dimension(0:3,4,4) :: bm
      complex, dimension(0:3,4):: pr
      pr = cmplx(0.0,0.0)
      do mu = 0,3
      do j = 1,4
      do i = 1,4
      pr(mu,j)  = pr(mu,j) + a(i)*bm(mu,j,i)
      enddo
      enddo
      enddo
*      print *,"*",a
*      print *,"***", pr(0,1:4)
      return
      end




      subroutine nucleon_vertex(nuc,q0,qx,qy,qz,bgamma)
***************************************************************************
*  This subroutine calculates the elastic electromagnetc vertex of nucleon
*   \Gamma_\mu = \gamma_mu F1 + i/(2M) F2 q^\nu \sigma_{\mu\nu}
* 
*   nuc = 1 proton, -1 neutron
*   q0,qx,qy,qz - energy and momentum components of virtual photon
*   bgamma(0:3,4,4)  \gamma_m  4x4 matrix for \mu = 0-3
*
*
* Misak Sargsian
* March 2007, FIU, Miami
******************************************************************************
      implicit none
      common/parms/pi,pm,dm
      common/gammamatices_s/gamma,sigma
      complex, dimension(0:3,4,4):: gamma,bgamma,sigmaq
      complex, dimension(0:3,0:3,4,4) :: sigma
      complex fdir_1,fdir_2
      real pi,pm,dm
      real q0,qx,qy,qz,qv2,q2,F1s,F2s
      integer nuc,mu,j,i
      qv2 = qx**2 + qy**2 + qz**2
      q2  = qv2 - q0**2
      do mu = 0,3
      do j = 1,4
      do i = 1,4
      sigmaq(mu,i,j) = q0*sigma(mu,0,i,j) - qx*sigma(mu,1,i,j)
     &                -qy*sigma(mu,2,i,j) - qz*sigma(mu,3,i,j)

      enddo
      enddo
      enddo
     
      fdir_1 = cmplx(F1s(nuc,q2), 0.0)
      fdir_2 = cmplx(0.0       , 1.0/(2.0*pm)*F2s(nuc,q2))

      bgamma =fdir_1*gamma + fdir_2*sigmaq
*      print *,"sigma",fdir_2,q2,pm
      return
      end





     
**********************************************************
*  Dirac Nucleon Form Factors for Nucleon
**********************************************************
      function F1s(nuc,q2)
*      change
      if(nuc.eq. 1)F1s = F1ps(q2)
      if(nuc.eq.-1)F1s = F1ns(q2)
      return
      end

      function F2s(nuc,q2)
*       change
      if(nuc.eq. 1)F2s = F2ps(q2)
      if(nuc.eq.-1)F2s = F2ns(q2)
      return
      end


      function F1ps(Q2)                                                  
      common/parms/pi,pm,dm
      common/formfcase/ics
      call proton_formfactors(Q2,GEp,GMp,ics)
      tau = Q2/4./pm**2             
      F1ps = (tau*GMp+GEp)/(1.+tau)                                        
      return                                                            
      end                                                            
C     ----------------------------------------                          
      function F2ps(Q2)                                                  
      common/parms/pi,pm,dm
      common/formfcase/ics
      call proton_formfactors(Q2,GEp,GMp,ics)
      tau = Q2/4./pm**2                                                   
      F2ps = (GMp-GEp)/(1.+tau)    
      return                                             
      end                                                            
C     ---------------------------------------------                     
      function F1ns(Q2)                                                  
      common/parms/pi,pm,dm
      common/formfcase/ics
      call neutron_formfactors(Q2,GEn,GMn,ics)
      tau = Q2/4./pm**2                                         
      F1ns = (tau*GMn+GEn)/(1.+tau)                                        
      return                                                            
      end                                                            
C     ---------------------------------------------                     
      function F2ns(Q2)                                                  
      common/parms/pi,pm,dm
      common/formfcase/ics
      call neutron_formfactors(Q2,GEn,GMn,ics)
      tau = Q2/4./pm**2                                                   
      F2ns = (GMn-GEn)/(1.+tau)                                            
      return                                                            
      end                                                            
C     -------------------------------------                             
*******************************************************
*   Charge and Magnetic Form Factor Parameterizations
*   Standard parameterization (one I used before)
*******************************************************
      subroutine proton_formfactors(q2,GEp,GMp,ics)
      common/parms/pi,pm,dm
      tau = Q2/4./pm**2                           

      if(ics.eq.1)then
************************************************************
*  SLAC Parameterization
************************************************************
      GEp = Gs(Q2)                                                         
      GMp = Gs(Q2)*2.79

      elseif(ics.eq.2.or.ics.ge.4)then
************************************************************
*    J.J. Kelly's Parameterization
*    From Phys. Rev. C70, 068202, 2004
************************************************************
      GEp = (1.-0.24*tau)/(1.+10.98*tau+12.82*tau**2+21.97*tau**3)
      GMp = 2.79*(1.+0.12*tau)/(1.+10.97*tau+18.86*tau**2+6.55*tau**3)

      elseif(ics.eq.3)then
************************************************************
*    Bradford, Bodek, Budd, Arrington
*    Hep-ex/0602017
************************************************************
      GEp = (1.-0.0578*tau)/(1.+11.1*tau+13.6*tau**2+33.0*tau**3)
      GMp = 2.79*(1.+0.15*tau)/(1.+11.1*tau+19.6*tau**2+7.54*tau**3)
***************************************************************************      
      endif

      return
      end

      subroutine neutron_formfactors(q2,GEn,GMn,ics)
      common/parms/pi,pm,dm
      tau = Q2/4./pm**2    

      if(ics.eq.1)then
************************************************************
*  SLAC Parameterization
************************************************************
      GEn =  1.91*Gs(Q2)*tau/(1.+5.6*tau)                                   
      GMn = -1.91*Gs(Q2) 

      elseif(ics.eq.2.or.ics.ge.4)then
************************************************************
*    J.J. Kelly's Parameterization
*    From Phys. Rev. C70, 068202, 2004
************************************************************
      GEn0=  1.7*tau/(1.+3.3*tau)*Gs(Q2)                                  
      GMn = -1.91*(1.+2.33*tau)/(1.+14.72*tau+24.20*tau**2+84.1*tau**3)

      GEn = Gen0
      if(ics.eq.4)GEn = 0.0
      if(ics.eq.5)GEn = 2.0*GEn0

      elseif(ics.eq.3)then   
************************************************************
*    Bradford, Bodek, Budd, Arrington
*    Hep-ex/0602017
************************************************************
      GEn = (1.25*tau+1.3*tau**2)/(1.-9.86*tau+305.*tau**2-758.0*tau**3
     &                                                    +802.0*tau**4)

      GMn = -1.91*(1.+1.81*tau)/(1.+14.1*tau+20.7*tau**2+68.7*tau**3)
      endif

      return
      end



      function Gs(Q2)                                                    
      Gs=1./(1.+Q2/0.71)**2                                               
      return                                                            
      end                                                            





      subroutine us_s(m,px,py,pz,e,s,isp,us,usc)
      implicit none
      complex, dimension(4)   :: us, usc
      real,    dimension(3)   :: s
      real,    dimension(0:3) :: s4
      real m,px,py,pz,e, s_dot_p
      complex p_pl,p_mn
      complex s_0,s_pl,s_mn,s_z
      integer isp
      real eplm,norm

      p_pl = cmplx(px,py)
      p_mn = cmplx(px,-py)
*      e = sqrt(m**2 + px**2 + py**2 + pz**2)
      eplm = e+m

      s_dot_p = s(1)*px + s(2)*py + s(3)*pz

      s4(0) = s_dot_p/m
      s4(1) = s(1) + s_dot_p/(m*(e+m)) * px
      s4(2) = s(2) + s_dot_p/(m*(e+m)) * py
      s4(3) = s(3) + s_dot_p/(m*(e+m)) * pz

      s_0   = cmplx(s4(0),0.0)
      s_pl  = cmplx(s4(1),s4(2))
      s_mn  = cmplx(s4(1),-s4(2))
      s_z   = cmplx(s4(3),0.0)
      norm = sqrt(2.0/(1.0+s(3)))*sqrt(eplm)/2.0
      if(isp.eq.1)then
      us(1)=norm*(1.0  + s_z - s_0*pz/eplm)
      us(2)=norm*(s_pl - s_0*p_pl/eplm)
      us(3)=norm*(s_0 + (1.0-s_z)*pz/eplm-s_mn*p_pl/eplm)
      us(4)=norm*(-s_pl*pz/eplm + (1.0 + s_z)*p_pl/eplm)

      usc(1)=norm*(1.0  + s_z - s_0*pz/eplm)
      usc(2)=norm*(s_mn - s_0*p_mn/eplm)
      usc(3)=norm*(s_0 +(1.0-s_z)*pz/eplm-s_pl*p_mn/eplm)
      usc(4)=norm*(-s_mn*pz/eplm + (1.0+s_z)*p_mn/eplm)


      elseif(isp.eq.-1.0)then
      us(1)=norm*(-s_mn + s_0*p_mn/eplm)
      us(2)=norm*(1.0 + s_z - s_0*pz/eplm)
      us(3)=norm*((1.0 +s_z)*p_mn/eplm  - s_mn*pz/eplm)
      us(4)=norm*(-s_0 + s_pl*p_mn/eplm - (1.0-s_z)*pz/eplm)

      usc(1)=norm*(-s_pl + s_0*p_pl/eplm)
      usc(2)=norm*(1.0 + s_z - s_0*pz/eplm)
      usc(3)=norm*((1.0 +s_z)*p_pl/eplm - s_pl*pz/eplm)
      usc(4)=norm*(-s_0 +s_mn*p_pl/eplm - (1.0-s_z)*pz/eplm)
      endif

      return
      end



      subroutine u_s(m,px,py,pz,e,isp,u,uc)
***************************************************************
*   This subroutine calculates Dirac Spinor
*   m - mass
*   px,py,pz,e - three components of momenta and energy
*   isp - 1 spin - up  -1 spin down
*
*   u   - spinor
*   uc  - complex conjugate spinor
*
*
*
* Misak Sargsian
* March 2007, FIU, Miami
*
***************************************************************
      implicit none
      complex, dimension(4):: u,uc
      real m,px,py,pz,e,fct
      integer isp
      fct = sqrt(e+m)
      if(isp.eq.1)then
      u(1) = fct*cmplx(1.0,0.0)
      u(2) =     cmplx(0.0,0.0)
      u(3) = fct/(e+m)*cmplx(pz,0)
      u(4) = fct/(e+m)*cmplx(px,py)

      uc(1) = fct*cmplx(1.0,0.0)
      uc(2) =     cmplx(0.0,0.0)
      uc(3) = fct/(e+m)*cmplx(pz,0)
      uc(4) = fct/(e+m)*cmplx(px,-py)

      elseif(isp.eq.-1)then
      u(1) =     cmplx(0.0,0.0)
      u(2) = fct*cmplx(1.0,0.0)
      u(3) = fct/(e+m)*cmplx(px,-py)
      u(4) = fct/(e+m)*cmplx(-pz,0.0)

      uc(1) =     cmplx(0.0,0.0)
      uc(2) = fct*cmplx(1.0,0.0)
      uc(3) = fct/(e+m)*cmplx(px,py)
      uc(4) = fct/(e+m)*cmplx(-pz,0.0)
      endif
      return
      end


      


      subroutine e_vector(e_vec)
      implicit none
      complex, dimension(-1:1,3) :: e_vec
      complex im,one,zero
      zero = cmplx(0.0,0.0)
      one  = cmplx(1.0,0.0)      
      im   = cmplx(0.0,1.0)
      e_vec(-1,1:3) =  1.0/sqrt(2.0)*(/one, -im  ,zero/)
      e_vec(0,1:3)  =  1.0/sqrt(2.0)*(/zero, zero, one/)
      e_vec(1,1:3)  = -1.0/sqrt(2.0)*(/one,  im  ,zero/)
      return
      end




      subroutine clebsh_gordan(i1,mi1,i2,mi2,i,m,coeff)
*****************************************************************
*      Calculates the Clebsh-Gordan Coefficients
*      spins and projections defined as multiplied by 2
*      for example -1/2 is -1,  -3/2 is -3
*
*      Misak Sargsian
*      March 2007, FIU, Miami
*     
******************************************************************
      implicit none 
      integer i1,mi1,i2,mi2,i,m
      real coeff
      coeff = 0.0
      if(i1.eq.2.and.i2.eq.1.and.i.eq.3)then
      coeff = 0.0
      if(mi1.eq. 2.and.mi2.eq. 1.and.m.eq. 3)coeff = 1.0
      if(mi1.eq. 2.and.mi2.eq.-1.and.m.eq. 1)coeff = sqrt(1.0/3.0)
      if(mi1.eq. 0.and.mi2.eq. 1.and.m.eq. 1)coeff = sqrt(2.0/3.0)
      if(mi1.eq. 0.and.mi2.eq.-1.and.m.eq.-1)coeff = sqrt(2.0/3.0)
      if(mi1.eq.-2.and.mi2.eq. 1.and.m.eq.-1)coeff = sqrt(1.0/3.0)
      if(mi1.eq.-2.and.mi2.eq.-1.and.m.eq.-3)coeff = 1.0
      endif
      return
      end
      

      subroutine gammamatricess
*********************************************************************************
*     This subroutine calculates Dirac gamma and sigma 
*     matrices.
*
*    sigma_{\mu,\nu} = {1\over 2}[\gamma_\mu\gamma_nu - \gamma_\nu\gamma_\mu}
*    gamma(0:3,4,4) 0:3 - are the values of \mu
*    (4,4) components are such
*     a(4,4)
*     a11 a12 a13 a14
*     a21 a22 a23 a24
*     a31 a32 a33 a34
*     a41 a42 a43 a44
*
*   March 2007, FIU, Miami
*********************************************************************************
      complex, dimension(0:3,4,4) :: gamma
      complex, dimension(0:3,0:3,4,4) :: sigma
      common/gammamatices_s/gamma, sigma
      complex fct
***********************************************
*  gamma_0
***********************************************
      gamma(0,1,1) = cmplx(1.,0.) 
      gamma(0,1,2) = cmplx(0.,0.)
      gamma(0,1,3) = cmplx(0.,0.)
      gamma(0,1,4) = cmplx(0.,0.)

      gamma(0,2,1) = cmplx(0.,0.) 
      gamma(0,2,2) = cmplx(1.,0.)
      gamma(0,2,3) = cmplx(0.,0.)
      gamma(0,2,4) = cmplx(0.,0.)

      gamma(0,3,1) = cmplx(0.,0.) 
      gamma(0,3,2) = cmplx(0.,0.)
      gamma(0,3,3) = cmplx(-1.,0.)
      gamma(0,3,4) = cmplx(0.,0.)

      gamma(0,4,1) = cmplx(0.,0.) 
      gamma(0,4,2) = cmplx(0.,0.)
      gamma(0,4,3) = cmplx(0.,0.)
      gamma(0,4,4) = cmplx(-1.,0.)
***********************************************
*  gamma_x
***********************************************
      gamma(1,1,1) = cmplx(0.,0.) 
      gamma(1,1,2) = cmplx(0.,0.)
      gamma(1,1,3) = cmplx(0.,0.)
      gamma(1,1,4) = cmplx(1.,0.)

      gamma(1,2,1) = cmplx(0.,0.) 
      gamma(1,2,2) = cmplx(0.,0.)
      gamma(1,2,3) = cmplx(1.,0.)
      gamma(1,2,4) = cmplx(0.,0.)

      gamma(1,3,1) = cmplx(0.,0.) 
      gamma(1,3,2) = cmplx(-1.,0.)
      gamma(1,3,3) = cmplx(0.,0.)
      gamma(1,3,4) = cmplx(0.,0.)

      gamma(1,4,1) = cmplx(-1.,0.) 
      gamma(1,4,2) = cmplx(0.,0.)
      gamma(1,4,3) = cmplx(0.,0.)
      gamma(1,4,4) = cmplx(0.,0.)
***********************************************
*  gamma_y
***********************************************
      gamma(2,1,1) = cmplx(0.,0.) 
      gamma(2,1,2) = cmplx(0.,0.)
      gamma(2,1,3) = cmplx(0.,0.)
      gamma(2,1,4) = cmplx(0.,-1.)

      gamma(2,2,1) = cmplx(0.,0.) 
      gamma(2,2,2) = cmplx(0.,0.)
      gamma(2,2,3) = cmplx(0.,1.)
      gamma(2,2,4) = cmplx(0.,0.)

      gamma(2,3,1) = cmplx(0.,0.) 
      gamma(2,3,2) = cmplx(0.,1.)
      gamma(2,3,3) = cmplx(0.,0.)
      gamma(2,3,4) = cmplx(0.,0.)

      gamma(2,4,1) = cmplx(0.,-1.) 
      gamma(2,4,2) = cmplx(0.,0.)
      gamma(2,4,3) = cmplx(0.,0.)
      gamma(2,4,4) = cmplx(0.,0.)

***********************************************
*  gamma_y
***********************************************
      gamma(3,1,1) = cmplx(0.,0.) 
      gamma(3,1,2) = cmplx(0.,0.)
      gamma(3,1,3) = cmplx(1.,0.)
      gamma(3,1,4) = cmplx(0.,0.)

      gamma(3,2,1) = cmplx(0.,0.) 
      gamma(3,2,2) = cmplx(0.,0.)
      gamma(3,2,3) = cmplx(0.,0.)
      gamma(3,2,4) = cmplx(-1.,0.)

      gamma(3,3,1) = cmplx(-1.,0.) 
      gamma(3,3,2) = cmplx(0.,0.)
      gamma(3,3,3) = cmplx(0.,0.)
      gamma(3,3,4) = cmplx(0.,0.)

      gamma(3,4,1) = cmplx(0.,0.) 
      gamma(3,4,2) = cmplx(1.,0.)
      gamma(3,4,3) = cmplx(0.,0.)
      gamma(3,4,4) = cmplx(0.,0.)

      fct = cmplx(0.0, 0.5)
      do mu = 0,3
      do nu = 0,3
      sigma(mu,nu,1:4,1:4) = fct*
     &  (matmul(gamma(mu,1:4,1:4),gamma(nu,1:4,1:4)) - 
     &   matmul(gamma(nu,1:4,1:4),gamma(mu,1:4,1:4)) )
      enddo
      enddo

      return
      end
