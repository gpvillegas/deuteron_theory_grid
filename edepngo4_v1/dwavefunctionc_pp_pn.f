      subroutine dfunction(kj,kt,ksp,ksn,p,theta,phi,wf)
****************************************************************************************
*       Subroutine calculates the deuteron wave function by spin and isospin 
*       components for given deuteron spin projection
*       kj  = 1,0,-1 - deuteron spin projection
*       kt  = 1,-1 - struck out nucleon isospin 1-proton -1 neutron
*       ksp = 1,-1 - spin projection of proton
*       ksn = 1,-1 - spin projection of neutron
*       p, theta, phi - momentum(GeV/c), polar (rad) and azimuthal(rad) angles of 
*                       relative p-n momentum
*       wf - complex - wave function for given spin 
*      20-July-2003  
*      Created
*      16-Feb-2005
*      Modified: complex structure is added
*      30-Nov-2006
*      common/pp_or_pn/ipp
*      add to switch off d wave in wd, wd1 and wd2 to imitate the pp state SRC
*      20-Jan-2007
*      we added common/isd/isd, which makes sure that 
*      wave function changes the sign when proton replaced by neutron 
*      only for the case of the direct mechanism
* 
*      27-March-2010
*      error in the sign of kj=-1 wave function is corrected for up dn + un dp component 
*      of the wave function
*      
*      Miami
*
*****************************************************************************************
      complex wf
      common/isd/isd


                  fis =  1.0 ! struck proton
*      if(kt.eq.-1.and.isd.eq.1)fis = -1.0 ! struck neutron
      if(kt.eq.-1)fis = -1.0 ! struck neutron
                  

      if(kj.eq.1)then ! deuteron spin projection 1
          if(ksp.eq.1.and.ksn.eq.1)then ! spins up
      wf_re = uu(p) + wd(p)/sqrt(8.0)*(3.0*cos(theta)**2-1.0)
      wf_im = 0.0
      elseif(ksp.eq.1.and.ksn.eq.-1)then ! proton up neutron down
      wf_re = wd(p)/sqrt(8.0)*3.0*cos(theta)*sin(theta)*cos(phi)   
      wf_im = wd(p)/sqrt(8.0)*3.0*cos(theta)*sin(theta)*sin(phi)   
      elseif(ksp.eq.-1.and.ksn.eq.1)then ! proton down neutron up
      wf_re = wd(p)/sqrt(8.0)*3.0*cos(theta)*sin(theta)*cos(phi)   
      wf_im = wd(p)/sqrt(8.0)*3.0*cos(theta)*sin(theta)*sin(phi)   
      elseif(ksp.eq.-1.and.ksn.eq.-1)then ! proton down neutron down
      wf_re = wd(p)/sqrt(8.0)*3.0*sin(theta)**2*cos(2.0*phi)
      wf_im = wd(p)/sqrt(8.0)*3.0*sin(theta)**2*sin(2.0*phi)
      endif

      elseif(kj.eq.0)then ! deuteron spin projection 0
          if(ksp.eq.1.and.ksn.eq.1)then ! spins up
      wf_re =  wd(p)*3.0/2.0*cos(theta)*sin(theta)*cos(phi)
      wf_im = -wd(p)*3.0/2.0*cos(theta)*sin(theta)*sin(phi)
      elseif(ksp.eq.1.and.ksn.eq.-1)then ! proton up neutron down
      wf_re = uu(p)/sqrt(2.0) - wd(p)/2.0 * (3.0*cos(theta)**2-1.0)
      wf_im = 0.0
      elseif(ksp.eq.-1.and.ksn.eq.1)then ! proton down neutron up
      wf_re = uu(p)/sqrt(2.0) - wd(p)/2.0 * (3.0*cos(theta)**2-1.0)
      wf_im = 0.0
      elseif(ksp.eq.-1.and.ksn.eq.-1)then ! proton down neutron down
      wf_re = -wd(p)*3.0/2.0*cos(theta)*sin(theta)*cos(phi)
      wf_im = -wd(p)*3.0/2.0*cos(theta)*sin(theta)*sin(phi)
      endif

      elseif(kj.eq.-1)then ! deuteron spin projection -1
          if(ksp.eq.1.and.ksn.eq.1)then ! spins up
      wf_re =  wd(p)/sqrt(8.0)*3.0*sin(theta)**2*cos(2.0*phi)
      wf_im = -wd(p)/sqrt(8.0)*3.0*sin(theta)**2*sin(2.0*phi)
      elseif(ksp.eq.1.and.ksn.eq.-1)then ! proton up neutron down
      wf_re = -wd(p)/sqrt(8.0)*3.0*cos(theta)*sin(theta)*cos(phi)   
      wf_im =  wd(p)/sqrt(8.0)*3.0*cos(theta)*sin(theta)*sin(phi)        
      elseif(ksp.eq.-1.and.ksn.eq.1)then ! proton down neutron up
      wf_re = -wd(p)/sqrt(8.0)*3.0*cos(theta)*sin(theta)*cos(phi)   
      wf_im =  wd(p)/sqrt(8.0)*3.0*cos(theta)*sin(theta)*sin(phi)   
      elseif(ksp.eq.-1.and.ksn.eq.-1)then ! proton down neutron down
      wf_re = uu(p) + wd(p)/sqrt(8.0)*(3.0*cos(theta)**2-1.0)
      wf_im = 0.0
      endif
      endif
      wf_rea = fis*wf_re 
      wf_ima = fis*wf_im
      wf = cmplx(wf_rea,wf_ima)
      return
      end

      subroutine dfunctionc(kj,kt,ksp,ksn,p,theta,phi,wf_rea,wf_ima)
****************************************************************************************
*       Subroutine calculates the deuteron wave function by spin and isospin 
*       components for given deuteron spin projection
*       kj  = 1,0,-1 - deuteron spin projection
*       kt  = 1,-1 - struck out nucleon isospin 1-proton -1 neutron
*       ksp = 1,-1 - spin projection of proton
*       ksn = 1,-1 - spin projection of neutron
*       p, theta, phi - momentum(GeV/c), polar (rad) and azimuthal(rad) angles of 
*                       relative p-n momentum
*       wf_rea - real part of the wave function
*       wf_ima - imaginary part of the wave function
*
*      1-December-2003  
*      Miami
*
*****************************************************************************************
                  fis =  1.0 ! struck proton
      if(kt.eq.-1)fis = -1.0 ! struck neutron
      call zeta(kj,ksp,ksn,ze_re,ze_im)
      call stensor(kj,ksp,ksn,theta,phi,s_re,s_im)
      s8 = sqrt(8.0)

      wf_re = uu(p)*ze_re + wd(p)/s8 * s_re 
      wf_im =             + wd(p)/s8 * s_im 
      wf_rea = fis*wf_re 
      wf_ima = fis*wf_im
      return
      end


      subroutine dfunction_a(kj,kt,ksp,ksn,p,theta,phi,awf)
****************************************************************************************
*       Subroutine calculates the a-deuteron wave function by spin and isospin 
*       components for given deuteron spin projection
*       kj  = 1,0,-1 - deuteron spin projection
*       kt  = 1,-1 - knocked out nucleon isospin 1-proton -1 neutron
*       ksp = 1,-1 - spin projection of proton
*       ksn = 1,-1 - spin projection of neutron
*       p, theta, phi - momentum(GeV/c), polar (rad) and azimuthal(rad) angles of 
*                       relative p-n momentum
*       awf_rea - real part of the wave function
*       awf_ima - imaginary part of the wave function
*
*      20-July-2003  
*      Miami
*
*****************************************************************************************
      complex awf
      pz = p*cos(theta) 
      pmin = 0.005
                  fis =  1.0 ! knocked-out proton
      if(kt.eq.-1)fis = -1.0 ! knocked-out neutron

      call zeta(kj,ksp,ksn,ze_re,ze_im)

      call stensor(kj,ksp,ksn,theta,phi,s_re,s_im)
      theta0 = acos(-1.0)/2.0
      call stensor(kj,ksp,ksn,theta0,phi,s0_re,s0_im)
                           tn2 = 0.0
      if(cos(theta).ne.0.0)tn2 = (sin(theta)/cos(theta))**2
      
      pt = abs(p*sin(theta))
      s8 = sqrt(8.0)

      awf_re= uu1(p,pt)*ze_re + 
     &        wd1(p,pt)/s8*s_re + tn2*wd2(pt)/s8*(s_re-s0_re) 
      awf_im= wd1(p,pt)/s8*s_im + tn2*wd2(pt)/s8*(s_im-s0_im) 

      awf_rea = pz*(fis*awf_re)
      awf_ima = pz*(fis*awf_im)
      awf = cmplx(awf_rea,awf_ima)
      return
      end


      subroutine zeta(kj,ksp,ksn,ze_re,ze_im)
****************************************************************************************
*       Subroutine calculates the spin wave function of the s component of the deuteron 
*       for given deuteron spin projection
*       kj  = 1,0,-1 - deuteron spin projection
*       ksp = 1,-1 - spin projection of proton
*       ksn = 1,-1 - spin projection of neutron
*       ze_re - real part of the tensor 
*       ze_im - imaginary part of the tensor wave 
*
*      1-December-2003  
*      FIU
*      Miami
*
*****************************************************************************************
      ze_re = 0.0
      ze_im = 0.0
      if(kj.eq.1)then      ! deuteron spin projection 1
      if(ksp.eq.1.and.ksn.eq.1)ze_re = 1             ! spins up
      elseif(kj.eq.0)then  ! deuteron spin projection 0
      if(ksp.eq.1.and.ksn.eq.-1)ze_re = 1.0/sqrt(2.0)! proton up neutron down
      if(ksp.eq.-1.and.ksn.eq.1)ze_re = 1.0/sqrt(2.0)! proton down neutron up
      elseif(kj.eq.-1)then ! deuteron spin projection -1
      if(ksp.eq.-1.and.ksn.eq.-1)ze_re = 1.0         ! spins down
      endif
      return
      end



      subroutine stensor(kj,ksp,ksn,theta,phi,s_re,s_im)
****************************************************************************************
*       Subroutine calculates the tensor S by spin and isospin 
*       components for given deuteron spin projection
*       S = [3 (\sigma_p p)(\sigma_n p)/p^2 - \sigma_p\sigma_n]\zeta
*       kj  = 1,0,-1 - deuteron spin projection
*       kt  = 1,-1 - struck out nucleon isospin 1-proton -1 neutron
*       ksp = 1,-1 - spin projection of proton
*       ksn = 1,-1 - spin projection of neutron
*       p, theta, phi - momentum(GeV/c), polar (rad) and azimuthal(rad) angles of 
*                       relative p-n momentum
*       s_re - real part of the tensor 
*       s_im - imaginary part of the tensor wave 
*
*      1-December-2003  
*      FIU
*      Miami
*
*****************************************************************************************

      if(kj.eq.1)then ! deuteron spin projection 1
          if(ksp.eq.1.and.ksn.eq.1)then ! spins up
      s_re = 3.0*cos(theta)**2-1.0
      s_im = 0.0
      elseif(ksp.eq.1.and.ksn.eq.-1)then ! proton up neutron down
      s_re = 3.0*cos(theta)*sin(theta)*cos(phi)   
      s_im = 3.0*cos(theta)*sin(theta)*sin(phi)   
      elseif(ksp.eq.-1.and.ksn.eq.1)then ! proton down neutron up
      s_re = 3.0*cos(theta)*sin(theta)*cos(phi)   
      s_im = 3.0*cos(theta)*sin(theta)*sin(phi)   
      elseif(ksp.eq.-1.and.ksn.eq.-1)then ! proton down neutron down
      s_re = 3.0*sin(theta)**2*cos(2.0*phi)
      s_im = 3.0*sin(theta)**2*sin(2.0*phi)
      endif

      elseif(kj.eq.0)then ! deuteron spin projection 0
          if(ksp.eq.1.and.ksn.eq.1)then ! spins up
      s_re =  sqrt(2.0)*3.0*cos(theta)*sin(theta)*cos(phi)
      s_im = -sqrt(2.0)*3.0*cos(theta)*sin(theta)*sin(phi)
      elseif(ksp.eq.1.and.ksn.eq.-1)then ! proton up neutron down
      s_re = -sqrt(2.0)*(3.0*cos(theta)**2-1.0)
      s_im =  0.0
      elseif(ksp.eq.-1.and.ksn.eq.1)then ! proton down neutron up
      s_re = -sqrt(2.0)*(3.0*cos(theta)**2-1.0)
      s_im = 0.0
      elseif(ksp.eq.-1.and.ksn.eq.-1)then ! proton down neutron down
      s_re = -sqrt(2.0)*3.0*cos(theta)*sin(theta)*cos(phi)
      s_im = -sqrt(2.0)*3.0*cos(theta)*sin(theta)*sin(phi)
      endif

      elseif(kj.eq.-1)then ! deuteron spin projection -1
          if(ksp.eq.1.and.ksn.eq.1)then ! spins up
      s_re =  3.0*sin(theta)**2*cos(2.0*phi)
      s_im = -3.0*sin(theta)**2*sin(2.0*phi)
      elseif(ksp.eq.1.and.ksn.eq.-1)then ! proton up neutron down
      s_re = -3.0*cos(theta)*sin(theta)*cos(phi)   
      s_im =  3.0*cos(theta)*sin(theta)*sin(phi)        
      elseif(ksp.eq.-1.and.ksn.eq.1)then ! proton down neutron up
      s_re = -3.0*cos(theta)*sin(theta)*cos(phi)   
      s_im =  3.0*cos(theta)*sin(theta)*sin(phi)   
      elseif(ksp.eq.-1.and.ksn.eq.-1)then ! proton down neutron down
      s_re = 3.0*cos(theta)**2-1.0
      s_im = 0.0
      endif
      endif
      return
      end



      function uu(p)
      common/which_wave/iw
      if(iw.eq.1)uu = U(p/0.197328)/sqrt(0.197328**3)
      if(iw.eq.2)uu = U_V18(p/0.197328)/sqrt(0.197328**3)
      if(iw.eq.3)uu = uub(p/0.197328)/sqrt(0.197328**3)
      if(iw.eq.4)uu = U_V18sb(p/0.197328)/sqrt(0.197328**3)
      return
      end
      function uu1(p,pt)
      common/which_wave/iw
      x  = p /0.197328
      xt = pt/0.197328
      if(iw.eq.1)uu1 = U_A(x,xt)/sqrt(0.197328**3)/0.197328
      if(iw.eq.2)uu1 = U_A_V18(x,xt)/sqrt(0.197328**3)/0.197328
      if(iw.eq.3)uu1 = uub_A(x,xt)/sqrt(0.197328**3)/0.197328
      if(iw.eq.4)uu1 = U_A_V18sb(x,xt)/sqrt(0.197328**3)/0.197328
*      uu1 = U_A(x,xt)/0.197328**(5.0/2.0)
      return
      end

      function wd(p)
      common/which_wave/iw
      common/pp_or_pn/ipp
      if(iw.eq.1)wd = W(p/0.197328)/sqrt(0.197328**3)
      if(iw.eq.2)wd = W_V18(p/0.197328)/sqrt(0.197328**3)
      if(iw.eq.3)wd = wwb(p/0.197328)/sqrt(0.197328**3)
      if(iw.eq.4)wd = W_V18sb(p/0.197328)/sqrt(0.197328**3)
      if(ipp.eq.1)wd = 0.0
      return
      end
      function wd1(p,pt)
      common/which_wave/iw
      common/pp_or_pn/ipp
      x  = p /0.197328
      xt = pt/0.197328
      if(iw.eq.1)wd1 = W_A(x,xt)/sqrt(0.197328**3)/0.197328
      if(iw.eq.2)wd1 = W_A_V18(x,xt)/sqrt(0.197328**3)/0.197328
      if(iw.eq.3)wd1 = wwb_A(x,xt)/sqrt(0.197328**3)/0.197328
      if(iw.eq.4)wd1 = W_A_V18sb(x,xt)/sqrt(0.197328**3)/0.197328
      if(ipp.eq.1)wd1 = 0.0
      return
      end

      function wd2(pt)
      common/which_wave/iw
      common/pp_or_pn/ipp
      xt = pt/0.197328
      if(iw.eq.1)wd2 = W_AA(xt)/sqrt(0.197328**3)/0.197328
      if(iw.eq.2)wd2 = W_AA_V18(xt)/sqrt(0.197328**3)/0.197328
      if(iw.eq.3)wd2 = wwb_AA(xt)/sqrt(0.197328**3)/0.197328
      if(iw.eq.4)wd2 = W_AA_V18sb(xt)/sqrt(0.197328**3)/0.197328
      if(ipp.eq.1)wd2 = 0.0
      return
      end



***************************************************************************************
* Paris wave function
***************************************************************************************
      FUNCTION fd_Paris(X,i)                                                  
C ************************************************                      
C *  DEUTRON WAVE FUNCTION WITH PARIS POTENTIAL  *                      
C ************************************************                      
      COMMON/PARIS/C(13),D(13),BM(13)                                   
      if(i.eq.1)then
      C(1)=0.88688076                                                   
      C(2)=-0.34717093                                                  
      C(3)=-3.050238                                                    
      C(4)=56.207766                                                    
      C(5)=-749.57334                                                   
      C(6)=5336.5279                                                    
      C(7)=-22706.863                                                   A1507090
      C(8)=60434.4690                                                   A1507100
      C(9)=-102920.58                                                   A1507110
      C(10)=112233.57                                                   A1507120
      C(11)=-75925.226                                                  A1507130
      C(12)=29059.715                                                   A1507140
      A=0.                                                              A1507150
      DO 401 J=1,12                                                     A1507160
401   A=A+C(J)                                                          A1507170
      C(13)=-A                                                          A1507180

      D(1)=0.023135193                                                  A1507190
      D(2)=-0.85604572                                                  A1507200
      D(3)=5.6068193                                                    A1507210
      D(4)=-69.462922                                                   A1507220
      D(5)=416.31118                                                    A1507230
      D(6)=-1254.6621                                                   A1507240
      D(7)=1238.783                                                     A1507250
      D(8)=3373.9172                                                    A1507260
      D(9)=-13041.151                                                   A1507270
      D(10)=19512.524                                                   A1507280
      DO 402 J=1,13                                                     A1507290
402   BM(J)=0.23162461+(J-1)                                            A1507300
      A=0.                                                              A1507310
      B=0.                                                              A1507320
      CC=0.                                                             A1507330
      DO 3 J=1,10                                                       A1507340
      A=A+D(J)/BM(J)**2                                                 A1507350
      B=B+D(J)                                                          A1507360
3     CC=CC+D(J)*BM(J)**2                                               A1507370
      D(11)=BM(11)**2/(BM(13)**2-BM(11)**2)/(BM(12)**2-BM(11)           A1507380
     ***2)*(-BM(12)**2*BM(13)**2*A+(BM(12)**2+BM(13)**2)*B-CC)          A1507390
      D(12)=BM(12)**2/(BM(11)**2-BM(12)**2)/(BM(13)**2-BM(12)           A1507400
     ***2)*(-BM(13)**2*BM(11)**2*A+(BM(13)**2+BM(11)**2)*B-CC)          A1507410
      D(13)=BM(13)**2/(BM(12)**2-BM(13)**2)/(BM(11)**2-BM(13)           A1507420
     ***2)*(-BM(11)**2*BM(12)**2*A+(BM(11)**2+BM(12)**2)*B-CC)          A1507430
      endif
      fd_Paris=(U(X/0.197328)**2+W(X/0.197328)**2)/0.197328**3          A1507440
      RETURN                                                            A1507450
      END                                                               A1507460

C ***** S PARTIAL WAVE ******                                           
      FUNCTION U(X)                                                     A1507480
      COMMON/PARIS/C(13),D(13),BM(13)                                   A1507490
      A=0.                                                              A1507500
      DO 1 J=1,13                                                       A1507510
1     A=C(J)/(X*X+BM(J)**2)+A                                           A1507520
      F=0.79788456                                                      A1507530
      U=A*F/SQRT(4.*3.14159265)                                         A1507540
      RETURN                                                            A1507550
      END                                                               A1507560
C  **** D PARTIAL WAVE *****                                            
      FUNCTION W(X)                                                     A1507580
      COMMON/PARIS/C(13),D(13),BM(13)                                   A1507590
      A=0.                                                              A1507600
      DO 1 J=1,13                                                       A1507610
1     A=D(J)/(X*X+BM(J)**2)+A                                           A1507620
      F=0.79788456                                                      A1507630
      W=A*F/SQRT(4.*3.14159265)                                         A1507640
      RETURN                                                            A1507650
      END                                                               A1507660




C ***** S' PARTIAL WAVE ******                                           
      FUNCTION U_A(X,XT)                                              
      COMMON/PARIS/C(13),D(13),BM(13)                                    
      U_A = 0.0
      A=0.                                                               
      DO 1 J=1,13                                                        
1     A=C(J)/(X*X+BM(J)**2)/sqrt(XT**2+BM(J)**2)+A                       
      F=0.79788456                                                       
      U_A=A*F/SQRT(4.*3.14159265)                                        
      RETURN                                                             
      END                                                                


C  **** D' 1  PARTIAL WAVE *****                                            
      FUNCTION W_A(X,XT)                                                     
      COMMON/PARIS/C(13),D(13),BM(13)                                   
      A=0.                                                              
      DO 1 J=1,13                                                       
1     A=D(J)/(X*X+BM(J)**2)/sqrt(XT**2+BM(J)**2)+ A                   
      F=0.79788456                                                      
      W_A=A*F/SQRT(4.*3.14159265)                                         
      RETURN                                                            
      END                                                               

C  **** D' 2  PARTIAL WAVE *****                                            
      FUNCTION W_AA(XT)                                                     
      COMMON/PARIS/C(13),D(13),BM(13)                                   
      A=0.                                                              
      DO 1 J=1,13                                                       
1     A=D(J)/(BM(J)**2)/sqrt(XT**2+BM(J)**2)+ A                   
      F=0.79788456                                                      
      W_AA = A*F/SQRT(4.*3.14159265)                                         
      RETURN                                                            
      END                                                               

C  **** D PARTIAL WAVE *****                                             
      FUNCTION W_Askzbi(X,XT)                                                 
      COMMON/PARIS/C(13),D(13),BM(13)                                    
      W_Askzbi = 0.0
      A=0.     
      IF(X.GT.0.0)THEN
      DO 1 J=1,13    
      TERM1 = (3.0*XT**2+2.0*BM(J)**2)/(X*X+BM(J)**2)
     &        /sqrt(XT**2+BM(J)**2)
      TERM2 = 0.0 !3.0*XT/X**2
      TERM  = TERM1   !- TERM2
      A   = D(J)/BM(J)**2 * TERM + A
1     CONTINUE
      ENDIF
      F=0.79788456                                                       
      W_Askzbi=A*F/SQRT(4.*3.14159265)                                        
      RETURN                                                             
      END                                                                

***************************************************************************************
* V18 wave function
***************************************************************************************
 
 
      FUNCTION fd_V18(X,i)                                            
C ************************************************                      
C *  DEUTRON WAVE FUNCTION WITH PARIS POTENTIAL  *                      
C ************************************************                      
      COMMON/V18/C(12),D(12),BM(12)                                   
      if(i.eq.1)then
      C(1)  =  0.706699E+00                                                
      C(2)  = -0.169743E+00                                             
      C(3)  =  0.112368E+01                                             
      C(4)  = -0.852995E+01                                             
      C(5)  =  0.195033E+02                                             
      C(6)  = -0.757831E+02                                             
      C(7)  =  0.283739E+03                                             
      C(8)  = -0.694734E+03                                             
      C(9)  =  0.885257E+03                                             
      C(10) = -0.720739E+03                                             
      C(11) =  0.412969E+03
      C(12) = -0.103336E+03                                             

      D(1)  =  0.176655E-01                                               
      D(2)  = -0.124551E+00                                            
      D(3)  = -0.108815E+01                                            
      D(4)  =  0.384848E+01                                            
      D(5)  = -0.852442E+01                                            
      D(6)  =  0.209435E+02                                            
      D(7)  = -0.490728E+02                                            
      D(8)  =  0.577382E+02                                            
      D(9)  = -0.127114E+01                                            
      D(10) = -0.628361E+02
      D(11) =  0.581016E+02
      D(12) = -0.177062E+02                                             

      BM(1)  = 0.2316 
      BM(2)  = 1.0
      BM(3)  = 1.5 
      BM(4)  = 2.0
      BM(5)  = 2.5
      BM(6)  = 3.5
      BM(7)  = 4.5
      BM(8)  = 5.5
      BM(9)  = 6.5
      BM(10) = 8.0
      BM(11) = 9.5
      BM(12) = 11.0
      fd_cebaf = 1.0
      fd_V18   = 0.0

      return
      else
      fd_V18=(U_V18(X/0.197328)**2+W_V18(X/0.197328)**2)/0.197328**3          
      endif
      RETURN                                                            
      END                                                               
C ***** S PARTIAL WAVE ******                                           
      FUNCTION U_V18(X)                                                     
      COMMON/V18/C(12),D(12),BM(12)                                   
      A=0.                                                              
      DO 1 J=1,12                                                       
1     A   = C(J)/(X*X+BM(J)**2) + A                                           
      F   = 1.0 !0.79788456          
      U_V18 = A*F/SQRT(4.*3.14159265)                                         
      RETURN                                                            
      END                                                               
C  **** D PARTIAL WAVE *****                                            
      FUNCTION W_V18(X)                                                     
      COMMON/V18/C(12),D(12),BM(12)                                   
      A = 0.0                                                              
      DO 1 J=1,12                                                       
1     A   = D(J)/(X*X+BM(J)**2) + A                                           
      F   = 1.0 !0.79788456 
      W_V18 = A*F/SQRT(4.*3.14159265)                                         
      RETURN                                                            
      END         


C ***** S' PARTIAL WAVE ******                                           
      FUNCTION U_A_V18(X,XT)                                              
      COMMON/V18/C(12),D(12),BM(12)                                    
      U_A_V18 = 0.0
      A=0.                                                               
      DO 1 J=1,12                                                        
1     A=C(J)/(X*X+BM(J)**2)/sqrt(XT**2+BM(J)**2)+A                       
      F=1.0 !0.79788456                                                       
      U_A_V18 = A*F/SQRT(4.*3.14159265)                                        
      RETURN                                                             
      END                                                                


C  **** D' 1  PARTIAL WAVE *****                                            
      FUNCTION W_A_V18(X,XT)                                                     
      COMMON/V18/C(12),D(12),BM(12)                                   
      A=0.                                                              
      DO 1 J=1,12                                                       
1     A=D(J)/(X*X+BM(J)**2)/sqrt(XT**2+BM(J)**2)+ A                   
      F=1.0 !0.79788456                                                      
      W_A_V18 = A*F/SQRT(4.*3.14159265)                                         
      RETURN                                                            
      END                                                               

C  **** D' 2  PARTIAL WAVE *****                                            
      FUNCTION W_AA_V18(XT)                                                     
      COMMON/V18/C(12),D(12),BM(12)                                   
      A=0.                                                              
      DO 1 J=1,12                                                       
1     A=D(J)/(BM(J)**2)/sqrt(XT**2+BM(J)**2)+ A                   
      F=1.0 !0.79788456                                                      
      W_AA_V18 = A*F/SQRT(4.*3.14159265)                                         
      RETURN                                                            
      END                                                               

                                                      
***************************************************************************************
* cd Bonn wave function
***************************************************************************************
 	
 	function fd_bonn(p,i)
	common/bonn/c(11), d(11) , om(11)
	if(i.eq.1)then
	c(1) =   0.88472985
	c(2) = - 0.26408759 
	c(3) = - 0.44114404e-01
	c(4) = - 0.14397512e+02
	c(5) =   0.85591256e+02
	c(6) = - 0.31876761e+03
	c(7) =   0.70336701e+03
	c(8) = - 0.90049586e+03
	c(9) =   0.66145441e+03
	c(10)= - 0.25958894e+03
	a = 0
	do k = 1,10
	a = a + c(k)
	enddo
	c(11)= - a

	d(1) =   0.22623762e-01
	d(2) = - 0.50471056e+00
	d(3) =   0.56278897e+00
	d(4) = - 0.16079764e+02
	d(5) =   0.11126803e+03
	d(6) = - 0.44667490e+03
	d(7) =   0.10985907e+04
	d(8) = - 0.16114995e+04

	do j = 1,11
	om(j) = 0.2315380 + float(j-1)*0.9
	enddo
	tm0 = 0.0
	tm1 = 0.0
	tm2 = 0.0
	do k = 1,8
	tm0 = tm0 + d(k)
	tm1 = tm1 + d(k) / om(k)**2
	tm2 = tm2 + d(k) * om(k)**2
	enddo
	d(9) = om(9)**2/(om(11)**2-om(9)**2)/(om(10)**2-om(9)**2) *
     &         ( -om(10)**2*om(11)**2*tm1 + (om(10)**2+om(11)**2)*tm0
     &           -tm2)
	d(10) = om(10)**2/(om(9)**2-om(10)**2)/(om(11)**2-om(10)**2) *
     &         ( -om(11)**2*om(9)**2*tm1 + (om(11)**2+om(9)**2)*tm0
     &           -tm2)
	d(11) = om(11)**2/(om(10)**2-om(11)**2)/(om(9)**2-om(11)**2) *
     &         ( -om(9)**2*om(10)**2*tm1 + (om(9)**2+om(10)**2)*tm0
     &           -tm2)
	else
	fd_bonn = 0.0 !uub(p)**2 + wwb(p)**2
	endif
	return
	end



********************************************************
*     Deuteron wave function with Bonn potential
*     uub - s wave
*     wwb - d - wave
*     p is in fm^-1
*     wf fm^{3/2}  
*
********************************************************

	function uub(q)
	common/bonn/c(11), d(11) , om(11)
	u = 0.0
	do k = 1,11
	u = u + c(k)/(q**2 + om(k)**2)
	enddo
	uub = 1.0/sqrt(2.0)/acos(-1.0) * u  
	return
	end	

	function wwb(q)
	common/bonn/c(11), d(11) , om(11)
	u = 0.0
	do k = 1,11
	u = u + d(k)/(q**2 + om(k)**2)
	enddo
	wwb = 1.0/sqrt(2.0)/acos(-1.0) * u   
	return
	end	





C ***** S' PARTIAL WAVE ******                                           
      FUNCTION uub_A(X,XT)                                              
      common/bonn/C(11),D(11),BM(11)                                    
      uub_A = 0.0
      A=0.                                                               
      DO 1 J=1,11                                                        
1     A=C(J)/(X*X+BM(J)**2)/sqrt(XT**2+BM(J)**2)+A                       
      F=1.0 !0.79788456                                                       
      uub_A = A*F/sqrt(2.0)/acos(-1.0)                                        
      RETURN                                                             
      END                                                                


C  **** D' 1  PARTIAL WAVE *****                                            
      FUNCTION wwb_A(X,XT)                                                     
      COMMON/bonn/C(11),D(11),BM(11)                                   
      A=0.                                                              
      DO 1 J=1,11                                                       
1     A=D(J)/(X*X+BM(J)**2)/sqrt(XT**2+BM(J)**2)+ A                   
      F=1.0 !0.79788456                                                      
      wwb_A = A*F/sqrt(2.0)/acos(-1.0)                                         
      RETURN                                                            
      END                                                               

C  **** D' 2  PARTIAL WAVE *****                                            
      FUNCTION wwb_AA(XT)                                                     
      COMMON/bonn/C(11),D(11),BM(11)                                   
      A=0.                                                              
      DO 1 J=1,11                                                       
1     A=D(J)/(BM(J)**2)/sqrt(XT**2+BM(J)**2)+ A                   
      F=1.0 !0.79788456                                                      
      wwb_AA = A*F/sqrt(2.0)/acos(-1.0)                                       
      RETURN                                                            
      END                                                               


***************************************************************************************
* V18sb wave function  - Sabina's version of V18 wave function
***************************************************************************************
      FUNCTION fd_V18sb(X,i)                                            
C ************************************************                      
C *  DEUTRON WAVE FUNCTION WITH PARIS POTENTIAL  *                      
C ************************************************                      
      COMMON/V18sb/C(12),D(12),AM(12)                                   

      if(i.eq.1)then
         nmax = 12
         pi = acos(-1.0)
         
      am(1) = 0.232500e+00
      am(2) = 0.500000e+00
      am(3) = 0.800000e+00
      am(4) = 0.120000e+01
      am(5) = 0.160000e+01
      am(6) = 0.200000e+01
      am(7) = 0.400000e+01
      am(8) = 0.600000e+01
      am(9) = 0.100000e+02
      am(10) = 0.140000e+02
      am(11) = 0.180000e+02
      am(12) = 0.220000e+02

      c(1) =  0.105252223e+02
      c(2) =  0.124352529e+02
      c(3) = -0.687541641e+02
      c(4) =  0.239111042e+03
      c(5) = -0.441014422e+03
      c(6) =  0.300140328e+03
      c(7) = -0.230639939e+03
      c(8) =  0.409671540e+03
      c(9) = -0.733453611e+03
      c(10)=  0.123506081e+04
      c(11)= -0.120520606e+04
 
      c(12) = 0.0
      mx1 = nmax-1
      do k=1,mx1
         c(12) = c(12)-c(k)
      enddo


      d(1) =  0.280995496e+00
      d(2) =  0.334117629e-01
      d(3) = -0.727192237e+00
      d(4) = -0.302809607e+01
      d(5) = -0.903824982e+01
      d(6) =  0.496045967e+01
      d(7) = -0.271985613e+02
      d(8) =  0.125334598e+03
      d(9) = -0.346742235e+03


      d(10) = 0.0
      d(11) = 0.0
      d(12) = 0.0

      n3 = nmax-3
      sp2 = 0.0
      sm = 0.0
      sm2 = 0.0
      do lp = 1,n3
         sp2 = sp2+d(lp)/am(lp)**2
         sm = sm+d(lp)
         sm2 = sm2+d(lp)*am(lp)**2
      enddo
 
      a = am(12)**2
      b = am(11)**2
      cc = am(10)**2
      d(10) = (cc/((a-cc)*(b-cc)))*(-b*a*sp2+(b+a)*sm-sm2)
 
      a = am(10)**2
      b = am(12)**2
      cc = am(11)**2
      d(11) = (cc/((a-cc)*(b-cc)))*(-b*a*sp2+(b+a)*sm-sm2)
 
      a = am(11)**2
      b = am(10)**2
      cc = am(12)**2
      d(12) = (cc/((a-cc)*(b-cc)))*(-b*a*sp2+(b+a)*sm-sm2)
 

      do nn = 1,12
         c(nn) = c(nn)/(4.0*pi) /sqrt(pi/2.0)
         d(nn) = d(nn)/(4.0*pi) /sqrt(pi/2.0)
      enddo

      fd_V18sb   = 0.0

      return
      else
      fg = 0.197328
      fd_V18sb=(U_V18sb(X/fg)**2+W_V18sb(X/fg)**2)/fg**3    
      endif
      RETURN                                                            
      END                                                               
C

C ***** S PARTIAL WAVE ******                                           
      FUNCTION U_V18sb(X)                                                     
      COMMON/V18sb/C(12),D(12),BM(12)                                   
      A=0.                                                              
      DO 1 J=1,12                                                       
1     A   = C(J)/(X*X+BM(J)**2) + A                                           
      F   = 1.0 !0.79788456          
      U_V18sb = A*F/SQRT(4.*3.14159265)                                         
      RETURN                                                            
      END                                                               
C  **** D PARTIAL WAVE *****                                            
      FUNCTION W_V18sb(X)                                                     
      COMMON/V18sb/C(12),D(12),BM(12)                                   
      A = 0.0                                                              
      DO 1 J=1,12                                                       
1     A   = D(J)/(X*X+BM(J)**2) + A                                           
      F   = 1.0 !0.79788456 
      W_V18sb = A*F/SQRT(4.*3.14159265)                                         
      RETURN                                                            
      END         


C ***** S' PARTIAL WAVE ******                                           
      FUNCTION U_A_V18sb(X,XT)                                              
      COMMON/V18sb/C(12),D(12),BM(12)                                    
      U_A_V18sb = 0.0
      A=0.                                                               
      DO 1 J=1,12                                                        
1     A=C(J)/(X*X+BM(J)**2)/sqrt(XT**2+BM(J)**2)+A                       
      F=1.0 !0.79788456                                                       
      U_A_V18sb = A*F/SQRT(4.*3.14159265)                                        
      RETURN                                                             
      END                                                                


C  **** D' 1  PARTIAL WAVE *****                                            
      FUNCTION W_A_V18sb(X,XT)                                                 
      COMMON/V18sb/C(12),D(12),BM(12)                                   
      A=0.                                                              
      DO 1 J=1,12                                                       
1     A=D(J)/(X*X+BM(J)**2)/sqrt(XT**2+BM(J)**2)+ A                   
      F=1.0 !0.79788456                                                      
      W_A_V18sb = A*F/SQRT(4.*3.14159265)                                        
      RETURN                                                            
      END                                                               

C  **** D' 2  PARTIAL WAVE *****                                            
      FUNCTION W_AA_V18sb(XT)                                                  
      COMMON/V18sb/C(12),D(12),BM(12)                                   
      A=0.                                                              
      DO 1 J=1,12                                                       
1     A=D(J)/(BM(J)**2)/sqrt(XT**2+BM(J)**2)+ A                   
      F=1.0 !0.79788456                                                      
      W_AA_V18sb = A*F/SQRT(4.*3.14159265)                                      
      RETURN                                                            
      END                                                               
