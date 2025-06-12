      subroutine bb_dg(mu,kjd,q2,s,q,ktf,p_f,ksf,sp_f,p_r,ksr,
     &           p_fpr,p_rpr,ivr,noff,moff,bb)

*****************************************************************
*       noff  - for off-shellness of electromagnetic current
*       moff  - for off-shellness of NN amplitude
*****************************************************************
      implicit none
      common/par/pi,pm,pmp,pmn,dm,eb
      real pi,pm,pmp,pmn,dm,eb
*      common/vphoton/q2,q0,qv
      real q2,q0,qv
      real, dimension(3) :: q,p_f,p_r,p_fpr,p_rpr,p_1
      real, dimension(3) :: sp_f,sp_1
      real p_fm,p_rm,p_fprm,p_rprm,p_1m 
      real e_f,E_r,e_fpr_on,e_fpr_off,E_rpr
      real s,t_on,t_off,delta_pf2
      real andam1,andam
      real f_re,f_im
      real e1_on,e1_off,til_m2
      real p1m,th1,phi1
      real N
      integer ktf,kjd,ksf,ksr,ks1
      integer ivr,noff,moff
      complex, dimension(2,2) :: jNpr
      complex, dimension(0:3,2,2):: jNprm,jN_on,jN_off
      complex fNN,bb
      complex wf
      integer sf,s1
      integer its,ini
      integer mu
      qv     = sqrt(q(1)**2 + q(2)**2 + q(3)**2)
      q0     = sqrt(qv**2 - q2)
      p_rm   = sqrt(  p_r(1)**2 +   p_r(2)**2 +   p_r(3)**2)
      p_rprm = sqrt(p_rpr(1)**2 + p_rpr(2)**2 + p_rpr(3)**2)
      E_r    = sqrt(pm**2 + p_rm**2) 
      E_rpr  = sqrt(pm**2 + p_rprm**2)
      andam1 = sqrt(E_rpr/E_r)*sqrt(s*(s-4.0*pm**2))
      andam  = andam1/(2.0*qv*E_rpr)/N(p_rprm)
*      andam  = 1.0
      p_1    = - p_rpr
      p_1m   = p_rprm
*      print *,p_1, p_rpr
      
*****************************************************************
*  Definition of rescattering amplitudes
*****************************************************************
	p_fm   = sqrt(p_f(1)**2   + p_f(2)**2   + p_f(3)**2)
	e_f    = sqrt(pm**2 + p_fm**2)


	p_fprm    = sqrt(p_fpr(1)**2 + p_fpr(2)**2 + p_fpr(3)**2) 
	e_fpr_on  =      sqrt(pm**2+p_fprm**2)
	e_fpr_off = dm - sqrt(pm**2+p_rprm**2) + q0

	delta_pf2 = (p_f(1)-p_fpr(1))**2 + (p_f(2)-p_fpr(2))**2
     &                                   + (p_f(3)-p_fpr(3))**2

        t_on  =  (e_f - e_fpr_on)**2   - delta_pf2
*        t_on = -delta_pf2
*        t_off =  (e_f - e_fpr_off)**2  - delta_pf2 
*        t_off =  -delta_pf2 - (p_r(3)-p_rpr(3))**2
        

        til_m2 = e_fpr_off**2-p_fprm**2
*        if(til_m.gt.pm**2.or.til_m.lt.0.0)return
*        if(til_m.gt.pm**2)return
        t_off = t_on -  abs((pm**2 - til_m2))

*        print *,"tonof",t_on,t_off,til_m
 

        f_re = 0.0
        f_im = 0.0
        its = 21
	ini = 1
        if(moff.eq.0)then
	if(ivr.eq.1)call f_NN_on(s,t_on,f_re,f_im,its)
	if(ivr.eq.2)call f_NN_sd(ksf,ksr,ksf,ksr,s,t_on,f_re,f_im,
     &              its,ini)
*	fNN_on = cmplx(f_re,f_im)
        else
 
        f_re = 0.0
        f_im = 0.0

*        til_m2 = e_fpr_off**2-p_fprm**2

*        if(til_m2.lt.pm**2.and.til_m2.gt.0.0)then

*        print *,"tonof",t_on,t_off,til_m
*        t_off = t_on - (pm**2-til_m2)
        if(ivr.eq.1)call f_NN_off(s,t_off,f_re,f_im,its)
        if(ivr.eq.2)call f_NN_off(s,t_off,f_re,f_im,its)
*        endif

*	fNN_off = cmplx(f_re,f_im)
        endif
        fNN = cmplx(f_re,f_im)

**********************************************************************


      	
***********************************************************************
* Calculating nucleon's electromagnetic current
***********************************************************************
	e1_on  = sqrt(pm**2 + p_1m**2)
	e1_off = dm - e1_on
    

	sp_1(1) = 0.0; sp_1(2) = 0.0; sp_1(3) = 1.0

*        call  j_N(ktf,q2,qph,p_fpr,p_1,jNpr,mu)  
*             call  j_Nm(ktf,q2,qph,p_fpr,p_1,jNprm)   
        call j_N_matrix_s_onof(ktf,q2,q,p_fpr,sp_f,p_1,sp_1,
     &      jN_on,jN_off,e1_on,e1_off,noff)

	if(noff.eq.0)jNprm = jN_on
	if(noff.gt.0)jNprm = jN_on + jN_off

	jNpr = jNprm(mu,1:2,1:2)
***********************************************
*  Spin of the knocked-out nucleon
***********************************************
	             sf  = 1
	if(ksf.eq.-1)sf  = 2

	bb  = (0.0,0.0)

*********************************************
* Summing by the spin of the initial nucleon
*********************************************
      do ks1 = -1,1,2      
*******************************
	              s1 = 1
	 if(ks1.eq.-1)s1 = 2

************************************************
*  deuteron wave function
************************************************
  	call polar(p_1(1),p_1(2),p_1(3),p1m,th1,phi1)

        call dfunction(kjd,ktf,ks1,ksr,p1m,th1,phi1,wf)
        
        bb = bb + andam * fNN * jNpr(sf,s1)* wf
        enddo
*        bb = (1.0,1.0)
        return
        end

      subroutine bb_ndg(mu,kjd,q2,s,q,ktf,p_f,ksf,sp_f,p_r,ksr,
     &           p_fpr,p_rpr,noff,moff,bb)

*****************************************************************
*       noff  - for off-shellness of electromagnetic current
*       moff  - for off-shellness of NN amplitude
*****************************************************************
      implicit none
      common/par/pi,pm,pmp,pmn,dm,eb
      real pi,pm,pmp,pmn,dm,eb
*      common/vphoton/q2,q0,qv
      real q2,q0,qv
      real, dimension(3) :: q,p_f,p_r,p_fpr,p_rpr,p_1
      real, dimension(3) :: sp_f,sp_1
      real p_fm,p_rm,p_fprm,p_rprm,p_1m 
      real e_f,E_r,e_fpr_on,e_fpr_off,E_rpr
      real s,t_on,t_off,delta_pf2,til_m2
      real andam1,andam
      real f_re,f_im
      real e1_on,e1_off
      real p1m,th1,phi1
      real N
      integer ktf,kjd,ksf,ksr,ks1
      integer ksfpr,ksrpr
      integer noff,moff
      complex, dimension(2,2) :: jNpr
      complex, dimension(0:3,2,2):: jNprm,jN_on,jN_off
      complex fNN,bb
      complex wf
      integer s1,sfpr
      integer its,ini
      integer mu
      qv     = sqrt(q(1)**2 + q(2)**2 + q(3)**2)
      q0     = sqrt(qv**2 - q2)
      p_rm   = sqrt(  p_r(1)**2 +   p_r(2)**2 +   p_r(3)**2)
      p_rprm = sqrt(p_rpr(1)**2 + p_rpr(2)**2 + p_rpr(3)**2)
      E_r    = sqrt(pm**2 + p_rm**2) 
      E_rpr  = sqrt(pm**2 + p_rprm**2)
      andam1 = sqrt(E_rpr/E_r)*sqrt(s*(s-4.0*pm**2))
      andam  = andam1/(2.0*qv*E_rpr)/N(p_rprm)
      p_1    = - p_rpr
      p_1m   = p_rprm
*      print *,p_1, p_rpr
      
*****************************************************************
*  Definition of rescattering amplitudes
*****************************************************************
	p_fm   = sqrt(p_f(1)**2   + p_f(2)**2   + p_f(3)**2)
	e_f    = sqrt(pm**2 + p_fm**2)

	p_fprm    = sqrt(p_fpr(1)**2 + p_fpr(2)**2 + p_fpr(3)**2) 
	e_fpr_on  =      sqrt(pm**2+p_fprm**2)
	e_fpr_off = dm - sqrt(pm**2+p_rprm**2) + q0

	delta_pf2 = (p_f(1)-p_fpr(1))**2 + (p_f(2)-p_fpr(2))**2
     &                                   + (p_f(3)-p_fpr(3))**2

        t_on  =  (e_f - e_fpr_on)**2   - delta_pf2
*        t_on  = -delta_pf2
        t_off =  (e_f - e_fpr_off)**2  - delta_pf2
        til_m2 = e_fpr_off**2-p_fprm**2
*        if(til_m.gt.pm**2.or.til_m.lt.0.0)return
*        if(til_m.gt.pm**2)return
        t_off = t_on -  abs((pm**2 - til_m2))


*        print *,"tonof",t_on,t_off  	
***********************************************************************
* Calculating nucleon's electromagnetic current
*********************************************************************************
	e1_on  = sqrt(pm**2 + p_1m**2)
	e1_off = dm - e1_on
    
	sp_1(1) = 0.0; sp_1(2) = 0.0; sp_1(3) = 1.0

        call j_N_matrix_s_onof(ktf,q2,q,p_fpr,sp_f,p_1,sp_1,
     &                        jN_on,jN_off,e1_on,e1_off,noff)

	if(noff.eq.0)jNprm = jN_on
	if(noff.gt.0)jNprm = jN_on + jN_off

	jNpr = jNprm(mu,1:2,1:2)

*********************************************************************************
	bb  = (0.0,0.0)

*********************************************
* Summing by the spin of the initial nucleon
*********************************************
      do ks1 = -1,1,2      
*******************************
	              s1 = 1
	 if(ks1.eq.-1)s1 = 2
*****************************************************************************
* Summing by the spin of the nucleon after  interaction with virtual photon (f\prime) 
* knocked-out nucleon in the intermediate state before FSI 
*****************************************************************************
      do ksfpr = -1,1,2      
*****************************************************************************
	                sfpr = 1
	 if(ksfpr.eq.-1)sfpr = 2

*****************************************************************************
* Summing by the spin of the recoil nucleon in the intermediate state before FSI 
*****************************************************************************
      do ksrpr = -1,1,2      
******************************************************************************
************************************************
*  deuteron wave function
************************************************
  	call polar(p_1(1),p_1(2),p_1(3),p1m,th1,phi1)
        call dfunction(kjd,ktf,ks1,ksrpr,p1m,th1,phi1,wf)
        f_re = 0.0
        f_im = 0.0
        its = 21
	ini = 1
        if(moff.eq.0)then
        
    	call f_NN_sd(ksfpr,ksrpr,ksf,ksr,s,t_on,f_re,f_im,its,ini)
        if(ksfpr.ne.ksf.or.ksrpr.ne.ksr)then
        f_im = f_im/4.0
        f_re = f_re/4.0
        endif
        else
        f_re = 0.0
        f_im = 0.0
*        if(ksfpr.eq.ksf.and.ksrpr.eq.ksr)then
*        print *,"aha"
*        call f_NN_off(s,t_off,f_re,f_im,its)
    	call f_NN_sd(ksfpr,ksrpr,ksf,ksr,s,t_off,f_re,f_im,its,ini)
*        endif
        endif
        fNN = cmplx(f_re,f_im)

**********************************************************************
        
        bb = bb + andam * fNN * jNpr(sfpr,s1)* wf
        enddo
        enddo
        enddo
        return
        end

      subroutine bb_chex(mu,kjd,q2,s,q,ktf,p_f,ksf,sp_f,p_r,ksr,
     &           p_fpr,p_rpr,ivr,noff,moff,bb)

*****************************************************************
*       ivr - chex version
*       noff  - for off-shellness of electromagnetic current
*       moff  - for off-shellness of NN amplitude
*****************************************************************
      implicit none
      common/par/pi,pm,pmp,pmn,dm,eb
      real pi,pm,pmp,pmn,dm,eb
*      common/vphoton/q2,q0,qv
      real q2,q0,qv
      real, dimension(3) :: q,p_f,p_r,p_fpr,p_rpr,p_1
      real, dimension(3) :: sp_f,sp_1
      real p_fm,p_rm,p_fprm,p_rprm,p_1m 
      real e_f,E_r,e_fpr_on,e_fpr_off,E_rpr
      real s,t_on,t_off,u_on,u_off,delta_pf2,delta_pfr2,til_m2
      real andam1,andam
      real f_re,f_im
      real e1_on,e1_off
      real p1m,th1,phi1
      real N
      integer ktf,ktr,kjd,ksf,ksr,ks1
      integer ksfpr,ksrpr
      integer ivr,noff,moff
      complex, dimension(2,2) :: jNpr
      complex, dimension(0:3,2,2):: jNprm,jN_on,jN_off
      complex fNN,bb
      complex wf
      integer s1,sfpr
      integer its,ini
      integer mu
      ktr = -ktf
      qv     = sqrt(q(1)**2 + q(2)**2 + q(3)**2)
      q0     = sqrt(qv**2 - q2)
      p_rm   = sqrt(  p_r(1)**2 +   p_r(2)**2 +   p_r(3)**2)
      p_rprm = sqrt(p_rpr(1)**2 + p_rpr(2)**2 + p_rpr(3)**2)
      E_r    = sqrt(pm**2 + p_rm**2) 
      E_rpr  = sqrt(pm**2 + p_rprm**2)
      andam1 = sqrt(E_rpr/E_r)*sqrt(s*(s-4.0*pm**2))
      andam  = andam1/(2.0*qv*E_rpr)/N(p_rprm)
      p_1    = - p_rpr
      p_1m   = p_rprm
*      print *,p_1, p_rpr
      
*****************************************************************
*  Definition of rescattering amplitudes
*****************************************************************
	p_fm   = sqrt(p_f(1)**2   + p_f(2)**2   + p_f(3)**2)
	e_f    = sqrt(pm**2 + p_fm**2)

	p_fprm    = sqrt(p_fpr(1)**2 + p_fpr(2)**2 + p_fpr(3)**2) 
	e_fpr_on  =      sqrt(pm**2+p_fprm**2)
	e_fpr_off = dm - sqrt(pm**2+p_rprm**2) + q0

	delta_pf2 = (p_f(1)-p_fpr(1))**2 + (p_f(2)-p_fpr(2))**2
     &                                   + (p_f(3)-p_fpr(3))**2

        t_on  =  (e_f - e_fpr_on)**2   - delta_pf2
        t_off =  (e_f - e_fpr_off)**2  - delta_pf2
*        print *,"tonof",t_on,t_off  	

	delta_pfr2 = (p_r(1)-p_fpr(1))**2 + (p_r(2)-p_fpr(2))**2
     &                                    + (p_r(3)-p_fpr(3))**2

	u_on  = (e_r - e_fpr_on )**2 - delta_pfr2
	u_off = (e_r - e_fpr_off)**2 - delta_pfr2

        til_m2 = e_fpr_off**2-p_fprm**2
*        if(til_m.gt.pm**2.or.til_m.lt.0.0)return
*        if(til_m.gt.pm**2)return
        u_off = u_on -  abs((pm**2 - til_m2))

***********************************************************************
* Calculating nucleon's electromagnetic current
*********************************************************************************
	e1_on  = sqrt(pm**2 + p_1m**2)
	e1_off = dm - e1_on
    
	sp_1(1) = 0.0; sp_1(2) = 0.0; sp_1(3) = 1.0

        call j_N_matrix_s_onof(ktr,q2,q,p_fpr,sp_f,p_1,sp_1,
     &                        jN_on,jN_off,e1_on,e1_off,noff)

	if(noff.eq.0)jNprm = jN_on
	if(noff.gt.0)jNprm = jN_on + jN_off

	jNpr = jNprm(mu,1:2,1:2)

*********************************************************************************
	bb  = (0.0,0.0)

*********************************************
* Summing by the spin of the initial nucleon
*********************************************
      do ks1 = -1,1,2      
*******************************
	              s1 = 1
	 if(ks1.eq.-1)s1 = 2
*****************************************************************************
* Summing by the spin of the nucleon after  interaction with virtual photon (f\prime) 
* knocked-out nucleon in the intermediate state before FSI 
*****************************************************************************
      do ksfpr = -1,1,2      
*****************************************************************************
	                sfpr = 1
	 if(ksfpr.eq.-1)sfpr = 2

*****************************************************************************
* Summing by the spin of the recoil nucleon in the intermediate state before FSI 
*****************************************************************************
      do ksrpr = -1,1,2      
******************************************************************************
************************************************
*  deuteron wave function
************************************************
  	call polar(p_1(1),p_1(2),p_1(3),p1m,th1,phi1)
        call dfunction(kjd,ktr,ks1,ksrpr,p1m,th1,phi1,wf)
        f_re = 0.0
        f_im = 0.0
        its = 21
	ini = 1
        if(moff.eq.0)then
    	call f_pn_chex(ksfpr,ksrpr,ksr,ksf,s,u_on,f_re,f_im,ivr)
        else
    	call f_pn_chex(ksfpr,ksrpr,ksf,ksr,s,u_off,f_re,f_im,ivr)
        endif
        fNN = cmplx(f_re,f_im)

**********************************************************************
        
        bb = bb + andam * fNN * jNpr(sfpr,s1)* wf
        enddo
        enddo
        enddo
        return
        end

