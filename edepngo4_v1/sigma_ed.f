*      program edepn_v2
*      common/par/pi,pm,pmp,pmn,dm,eb
*      icon = -1
*      call dmatrix(q0,qv,kj,ktm,ktr,pr,thr,phir,spec,icon)


**************************************
*       INPUT   PARAMETERS   
**************************************
*        ktr = -1
*        ktm =  1
*        ei   = 5.76          ! initial electron energy
*        pr   = 0.4         ! recoil nucleon momentum GeV/C
*        thr  = 71.*pi/180.   ! polar angle of the recoil nucleon
*        phir = pi            ! azimuthal angle of the recoil nucleon

*        q2 = 4.0
*        x  = 1.0
*        q0 = q2/2.0/pm/x
*        qv = sqrt(q2 + qv**2)

*        icon = 12

*        sum = 0.0
*        do kj = -1,1,1
*        call dmatrix(q0,qv,kj,ktm,ktr,pr,thr,phir,spec,icon)
*        sum = sum + spec
*        enddo
*        sum_av = sum/3.0
*        write(6,*)sum_av,fd(pr,0)
*        end

      subroutine sigma_ed(ei,q2,pr,thr,phir,ktf,sp_f,noff,npv,npro,icon,
     &                iv1,iv2,isd,ics,iw,ipp,crs_pol,crs,crs_ten,ifl,q0)
************************************************************************************************
*  Code calculates exclusive d(e,e'N_f)N_f reaction in unfactorized approximation
************************************************************************************************
*         icon = 0 initialization, 1 - PWIA, 2 - forward rescattering, 3 - charge exchange rescattering 
*                12 - PWIA+FSI, 13 - PWIA+CHEX, 123 - PWIA+FSI+CHEX
*         iv2  - version of parameterization of CHEX amplitude - 
*                10 - Gibbs-Loiseau (has only real part), 21 - SAID, 1021 - SAID+GL
*                10021 - SAID+GL (with Im part at p>p0 modelled using SAID's value
*                -21 SAID   - For details look at fpn_chex.f
*         iv1  - version of FSI; 1- diagonal FSI using fpn.f, 2- diagonal FSI using f_NN_sd.f 
*                which uses SAID parameterization for plab<3.68, 3 - nondiagonal approximation 
*                using f_NN_sd.f. At p>1.4 it uses diffractive parametrization like fpn.f
*         isd  - isd = 1 spectator mechanism, isd=2 direct mechanism isd =0 both   
*         ics  - 1 - SLAC, 2 Kelly, 3 - Bodek, BB Arrington for-factor paremeterizations
*         iw   - 1 - Paris 2 - V18 deuteron wave funcion, 3- cd Bonn, 4 - AV18sb
*         (ics and iw used only during the initialization, with icon=0)
*         ipp  - 1 - quasi "pp" wave functio - only S-parital wave, 0, both S and D waves
*         ei   - initial electron energy (GeV)
*         q2   - virtual photon's Q2
*         pr   - slow nucleon momentum
*         thr  - slow nucleons polar angle
*         phir - slow nucleons azimuthal angle
*         ktf  - type of the fast nucleon 1 - proton -1 - neutron
*         sp_f(3) - polarization axis of  fast nucleon, if not polarized (0,0,1)
*         noff - off shellness of elm current,  0 - no off shell, 1 - off-shell 
*         npv -  Principal Value part of rescattering 0 - not included 1 - included 
*         npro - 1 ds/dE'd\Omega d\Omega_f,  2 - ds/dE'd\Omega d\Omega_r
*         crs_pol(j)- cross section for given polarization, j=-1,0,1
*         crs  - cross section in nb/GeV/str^2
*         crs_ten - cross section from scattering of tensor polarized deuteron
*         ifl  - flag, if 1- it is knematically forbidden
*         q0   - virtual photon energy calculated from QE condition
*         z || q
*
* version-3
* 12/09/09 npv is added. 
* version -4
* 05/04/23 - tensor polarization is added
*****************************************************************
	implicit none
	real ei,pr,thr,phir
	integer ktf,noff,npro,icon,iv1,iv2,isd,iw,ipp,ifl
        real crs0,crs
        common/par/pi,pm,pmp,pmn,dm,eb
	real pi,pm,pmp,pmn,dm,eb
        common/parms/pie,pme,dme  ! this goes to the calculation of j_N
	real pie,pme,dme
	common/alpha/alpha
	real alpha
*	common/MASS/PMe,DMe,PIe 
        common/electron/eil,esc,ue
	real eil,esc,ue
        common/photon/q2e,q0e,qve
	real q2e,q0e,qve

	real,    dimension(3):: p_f,p_r
	real,    dimension(3)   :: k1,k2
	complex, dimension(2,2) :: je,je_0,je_1,je_2,je_3
	complex, dimension(0:3,-1:1,2,2) :: jd, jd_s, jd_d
	complex jejd
        complex, dimension(0:3,2,2):: jN, jN_off
        real,    dimension(3) :: qii,p2ii,p1ii,sp_f,sp_r,sp_1,sp_2        
        common/pp_or_pn/ipp0
	integer ipp0
	common/isd/isd0
	integer isd0
        common/which_wave/iw0
	integer iw0
	common/spfmag/spfmag
	real spfmag
	real q2,qv,q0
	real p_i_x,p_i_y,p_i_z,p_i
	real p_f_x,p_f_y,p_f_z,pf,th_f,phi_f,e_f
	real p_r_x,p_r_y,p_r_z,er,e_r
	real sf0,sf
	real arg,sn_th_eq,cs_th_eq
        real sn_th_escq, cs_th_escq 
	real th_eq,phi_eq,phi_escq
	real vl_a,vl_b,s_un,sum,fkin,fc,fc1,absM2
	real e_off,e_on
	real sgn
	integer k1s,k2s
	integer kjd,ksf,ksr,ktr
	integer nuc,mu,nu
        complex, dimension(0:3,0:3) :: wa_s,wa_pol_d_mn1,wa_pol_d_0
        complex, dimension(0:3,0:3) :: wa_sum,wa_pol_d_pl1,wa_ten
        complex, dimension(0:3,0:3) :: wa_sum_pol
	real v_l,v_t,v_tl,v_tt
	real w_l,w_t,w_tl,w_tt
	real, dimension(-1:1):: w_l_pol,w_t_pol,w_tl_pol,w_tt_pol
	real, dimension(-1:1):: a_pol,crs_sf_pol,crs_pol
! common block to makes responses accessible
        common /responses/ w_l, w_t, w_tl, w_tt
 

	real sfc,a_un,crs_sf,sumpolmn1,sumpol0,sumpolpl1,s_ten
*    real sumpolmn1 

	integer npv,npv0
	common/npv/npv0
        integer ics0,ics
        common/formfcase/ics0
*tensor related definitions     

*    complex, dimension(0:3,0:3) :: wa_pol_d_mn1
*    complex, dimension(0:3,0:3) :: wa_pol_d_0
*    complex, dimension(0:3,0:3) :: wa_pol_d_pl1
*    complex, dimension(0:3,0:3) :: wa_ten
*    real sumpol0,sumpolpl1,s_ten   
 
 	real w_l_ten,w_t_ten,w_tl_ten,w_tt_ten
 	real a_ten,crs_ten,crs_sf_ten
   
 
***********************************************
*    passing flags to common space
***********************************************
	ipp0 = ipp
	isd0 = isd
	iw0  = iw
	npv0 = npv
*	common complex im
*	im = (0.,1.)

        crs0 = 0.0
        crs  = 0.0

	if(icon.eq.0)then
        ics0 = ics
        call j_d(q0,qv,p_f,ktf,sp_f,p_r,ktr,jd,icon,iv1,iv2,noff)
*        call dmatrix(q0,qv,kjd,ktm,ktr,pr,thr,phir,spec,icon)

************************************************
*   passing variables to the common space for 
*   the electromagnetic currents
************************************************
        pie = pi
	alpha = 1.0/137
        pme = pm
        dme = dm

        nuc = 0 ! initialization of electroganetic current
*	print *, "aha j_Nm"
	q2 = 0.0
	qii = 0.0
	p2ii = 0.0
	sp_2 = 0.0
	p1ii = 0.0
	sp_1 = 0.0
	e_on = 0.0
	e_off = 0.0
*         call  j_Nm(nuc,q2,qii,p2ii,p1ii,jN)
        call j_N_matrix_s_onof(nuc,q2,qii,p2ii,sp_2,p1ii,sp_1,jN,jN_off,
     &  e_on,e_off,noff)

*	print *,"jN,pie,pme,dme,nuc"
	crs = 0.0
        return
        endif
*********************************************************************
*  calculation of the magnitude of the  polarization vector
**********************************************************************
	spfmag = sp_f(1)**2 + sp_f(2)**2 + sp_f(3)**2

*************************************************
* passing electron energy to the common space
*************************************************
        eil = ei
*************************************************
*         definitions for slow nucleon
*************************************************
        ktr    = -ktf
	    p_r_x  =  pr  * sin(thr)*cos(phir) 
        p_r_y  =  pr  * sin(thr)*sin(phir)
        p_r_z  =  pr  * cos(thr)
*        print *,"aha",p_r_x,p_r_y,p_r_z,phir
        e_r = sqrt(pm**2 + pr**2)
	p_r(1) = p_r_x
	p_r(2) = p_r_y
	p_r(3) = p_r_z

************************************************
*       definition q0 by the QE condition
************************************************ 
	ifl    = 0
 	call kine_q0(q2,pr,thr,q0,ifl)
	if(ifl.eq.1)return
 	esc  = ei - q0         !scattered electron energy
	arg = sqrt(q2/4.0/ei/esc)
	if(arg**2.gt.1.0)then
	ifl = 1
	return 
	endif

	ue  = 2.0*asin(arg)

 	qv  = sqrt(q2 + q0**2)

********************************************************
* passing virtual photon momentum to the common space
********************************************************
        q2e = q2
        q0e = q0
        qve = qv

************************************** 
*       initial  nucleon
*       with respect to vrt. photon
**************************************
*        ktm    =  ktf 
	p_i_x  = -p_r_x 	
	p_i_y  = -p_r_y 	
        p_i_z  = -p_r_z
	p_i    = sqrt(p_i_x**2 + p_i_y**2 + p_i_z**2)

************************************** 
*       fast nucleon
*       with respect to vrt. photon
**************************************

        p_f_x = p_i_x
        p_f_y = p_i_y
        p_f_z = p_i_z + qv
        pf = sqrt(p_f_x**2 + p_f_y**2 + p_f_z**2)
        e_f = sqrt(pm**2+pf**2)
  	call polar(p_f_x,p_f_y,p_f_z,pf,th_f,phi_f)
	p_f(1) = p_f_x
	p_f(2) = p_f_y
	p_f(3) = p_f_z
**************************************
*        print *,"thf",th_f,phi_f
***********************************************
*         Definition of Electron Current
*         in the reference frame in which q||z
***********************************************
*	qv = sqrt(ei**2 - 2.0*ei*esc*cos(ue)+esc**2)

	cs_th_eq = (ei - esc*cos(ue))/qv
	sn_th_eq =       esc*sin(ue) /qv
	phi_eq = 0.0
	k1(1) = 0.0; k1(2) = 0.0; k1(3) = ei

	k1(1) = ei*sn_th_eq*cos(phi_eq)
	k1(2) = ei*sn_th_eq*sin(phi_eq)
	k1(3) = ei*cs_th_eq

	th_eq = acos(cs_th_eq)

	cs_th_escq = (ei*cos(ue)-esc)/qv
	sn_th_escq =  ei*sin(ue)/qv
	phi_escq   = 0.0
	if((th_eq + ue).ge.pi)then
	phi_escq   = pi
	endif

	if(phi_escq.gt.0.0)then
	print *,"attention"
	stop
	endif


	k2(1) = esc*sin(ue); k2(2) = 0.0; k2(3) = esc*cos(ue) 

	k2(1) = esc*sn_th_escq*cos(phi_escq)
	k2(2) = esc*sn_th_escq*sin(phi_escq)
	k2(3) = esc*cs_th_escq
 
*	sc_k2k1_a = ei*esc*cos(ue)
*	sc_k2k1_b = k2(1)*k1(1) + k2(2)*k1(2) + k2(3)*k1(3)
*	write(25,*)sc_k2k1_a,sc_k2k1_b
************************************************
*         \mu components of electron current
************************************************
	do mu = 0,3
	call  j_e(k2,k1,je,mu)
	if(mu.eq.0)je_0 = je
	if(mu.eq.1)je_1 = je
	if(mu.eq.2)je_2 = je
	if(mu.eq.3)je_3 = je
	enddo 

*****************************
*        Checking
***************************
	vl_a = 4.0*(2.0*esc*ei-q2/2.0)
	vl_b = 0.0
	do k2s=1,2
	do k1s=1,2
	vl_b = vl_b + je_0(k2s,k1s)*conjg(je_0(k2s,k1s))
	enddo
	enddo
*	print *,vl_a,vl_b,vl_a/vl_b
******************************

***************************************************************************
* Calculating Electromagnetic Current of the Deuteron Electrodisintegration
****************************************************************************
******************************************************************************
*  spectator mechanism  - in this case fast nucleon is the knocked-out nucleon
******************************************************************************
	jd_s = 0.0
	if(isd.eq.0.or.isd.eq.1)then
	ktr = -ktf
	if(ipp.eq.1)ktr = 1  ! case of pp production immitation
        call j_d(q0,qv,p_f,ktf,sp_f,p_r,ktr,jd_s,icon,iv1,iv2,noff)
	endif

****************************************************************
* direct mechanism  - in this case fast nucleons is a spectator
****************************************************************
	 jd_d =0.0
	 if(isd.eq.0.or.isd.eq.2)then
	 ktr = -ktf
	 if(ipp.eq.1)ktr = 1  ! case of pp production immitation
         call j_d(q0,qv,p_r,ktr,sp_f,p_f,ktf,jd_d,icon,iv1,iv2,noff)
	 endif

	             sgn =  1.0
	 if(ipp.eq.1)sgn = -1.0

	 do mu = 0,3
	 do kjd = -1,1,1
	 do ksf = 1,2
	 do ksr = 1,2
	 jd(mu,kjd,ksf,ksr) = jd_s(mu,kjd,ksf,ksr) 
     &                 +  sgn*jd_d(mu,kjd,ksr,ksf)
	 enddo
	 enddo
	 enddo
	 enddo


************************************************
* calculating unpolorized cross section
************************************************
	 sum = 0.0

*********************************************************************

         do kjd = -1,1,1  ! spin of the deuteron
	 do ksf = 1,2     ! spin of the fast nucleon  1 is up 2 is down 
	 do ksr = 1,2     ! spin of the slow nucleon  1 is up 2 is down	
	 do k2s = 1,2     ! spin of the scattered electron
	 do k1s = 1,2     ! spin of the initial electron

***************************************************************************
      jejd    =  je_0(k2s,k1s)*jd(0,kjd,ksf,ksr) 
     &          -je_1(k2s,k1s)*jd(1,kjd,ksf,ksr) 
     &          -je_2(k2s,k1s)*jd(2,kjd,ksf,ksr) 
     &          -je_3(k2s,k1s)*jd(3,kjd,ksf,ksr)
	 
	 absM2 = abs(jejd)**2
****************************************************************************
	 sum   = sum + absM2
	 enddo 
	 enddo  
	 enddo 
	 enddo 
	 if(kjd.eq.-1)sumpolmn1 = sum
	 if(kjd.eq.0) sumpol0   = sum
	 if(kjd.eq.1) sumpolpl1 = sum	 
	 enddo 
        s_un = (1.0/2.0)*(1.0/3.0)*sum
        s_ten = (1.0/2.0)*(1.0/3.0)*(sumpolmn1+sumpolpl1 - 2.0*sumpol0)
*****************************************************************************
*       print *,sum



         wa_sum   = (0.0,0.0)

***********************************************************************
*  Calculating electromagnetic tensor of he3 nucleus  W^\mu^\nu
***********************************************************************
         do kjd = -1,1,1  ! spin of the deuteron
         wa_sum_pol   = (0.0,0.0)

	 do ksf = 1,2     ! spin of the fast nucleon  1 is up 2 is down 
	 do ksr = 1,2     ! spin of the slow nucleon  1 is up 2 is down	



         do mu = 0,3
         do nu = 0,3
         wa_s(mu,nu) = conjg(jd(mu,kjd,ksf,ksr))*jd(nu,kjd,ksf,ksr)
         enddo
         enddo

         wa_sum   = wa_sum + wa_s
         wa_sum_pol = wa_sum_pol+wa_s
	 enddo 
	 enddo
	 if(kjd.eq.-1)then
	 wa_pol_d_mn1 = wa_sum_pol
*	  w_l_pol(kjd)  = Real(wa_sum_pol(0,0))
*	  w_t_pol(kjd)  = Real(wa_sum_pol(2,2) + wa_sum_pol(1,1))  
*      w_tt_pol(kjd) = Real(wa_sum_pol(1,1) - wa_sum_pol(2,2))
*      w_tl_pol(kjd) = -2.0d+0*Real(wa_sum_pol(0,1))

	 elseif(kjd.eq.0)then
	  wa_pol_d_0   = wa_sum_pol
*      w_l_pol(kjd)  = Real(wa_sum_pol(0,0))
*	  w_t_pol(kjd)  = Real(wa_sum_pol(2,2) + wa_sum_pol(1,1))  
*      w_tt_pol(kjd) = Real(wa_sum_pol(1,1) - wa_sum_pol(2,2))
*      w_tl_pol(kjd) = -2.0d+0*Real(wa_sum_pol(0,1))
	 elseif(kjd.eq.1)then
	  wa_pol_d_pl1 = wa_sum_pol
*      w_l_pol(kjd)  = Real(wa_sum_pol(0,0))
*	  w_t_pol(kjd)  = Real(wa_sum_pol(2,2) + wa_sum_pol(1,1))  
*      w_tt_pol(kjd) = Real(wa_sum_pol(1,1) - wa_sum_pol(2,2))
*      w_tl_pol(kjd) = -2.0d+0*Real(wa_sum_pol(0,1))
	endif
      w_l_pol(kjd)  = Real(wa_sum_pol(0,0))
	  w_t_pol(kjd)  = Real(wa_sum_pol(2,2) + wa_sum_pol(1,1))  
      w_tt_pol(kjd) = Real(wa_sum_pol(1,1) - wa_sum_pol(2,2))
      w_tl_pol(kjd) = -2.0d+0*Real(wa_sum_pol(0,1))




	 enddo
       
	 wa_ten = wa_pol_d_mn1+wa_pol_d_pl1 - 2.0*wa_pol_d_0
	 	 
	 
*	 print *,wa_sum(0,0)

         v_l  = q2**2/qv**4
         v_t  = q2/2.0/qv**2 + tan(ue/2.0)**2
         v_tt = q2/2.0/qv**2
         v_tl = q2/qv**2*sqrt(q2/qv**2 + tan(ue/2.0)**2)

         w_l  =         Real(wa_sum(0,0))
         w_t  =         Real(wa_sum(2,2) + wa_sum(1,1))
         w_tt =         Real(wa_sum(1,1) - wa_sum(2,2))
         w_tl = -2.0d+0*Real(wa_sum(0,1))

         w_l_ten  =         Real(wa_ten(0,0))
         w_t_ten  =         Real(wa_ten(2,2) + wa_ten(1,1))
         w_tt_ten =         Real(wa_ten(1,1) - wa_ten(2,2))
         w_tl_ten = -2.0d+0*Real(wa_ten(0,1))

*	     w_l_pol(-1:1),w_t_pol(-1:1),w_tt_pol(-1:1),w_tl_pol(-1:1),a_pol
* 
	     
	      
	     
	     
	     

*         print *,"wtt",w_l,w_t,w_tt,w_tl

         sfc  = 4.0*esc*ei*cos(ue/2.0)**2 /3.0

        a_un =   sfc*
     &          (v_l*w_l + v_t*w_t + v_tl*w_tl + v_tt*w_tt)

        a_ten =   sfc*
     &    (v_l*w_l_ten + v_t*w_t_ten + v_tl*w_tl_ten + v_tt*w_tt_ten)

	   do kjd = -1,1,1
	   a_pol(kjd) =   3.0*sfc*
     &    (v_l*w_l_pol(kjd)  + v_t*w_t_pol(kjd) + 
     &    v_tl*w_tl_pol(kjd) + v_tt*w_tt_pol(kjd))
       enddo

*	   print *,"aha",a_un,(a_pol(-1)+a_pol(0)+a_pol(1))/3.0

*        a_un =   sfc*
*     &          (v_l*w_l + v_t*w_t ) ! for Kim's kinematics



************************************** 
*       Kinematic factor
**************************************
	if(npro.eq.1)then
 	fc1 = (pf - qv*cos(th_f))/e_r
 	fc = abs(pf/e_f + fc1)
	fkin = pf**2 / fc
	endif
	if(npro.eq.2)then
 	fc1 = (pr - qv*cos(thr))/e_f
 	fc = abs(pr/e_r + fc1)
	fkin = pr**2 / fc
	endif
	sf0 = alpha**2/q2**2 * esc/ei *1.0/(2.0*dm*e_f) 
*         Above derivation is from page 8 of folder
*	sf0 = alpha**2/q2**2 * esc/ei *1.0/(4.0*e_r*e_f) 
	sf  = sf0*0.389385*1000.*1000.
*	relf = dm/(2.*(dm - e_r))

**************************************
*       Cross section 
**************************************
	crs  =   sf*fkin*s_un
	crs_sf = sf*fkin*a_un
	crs = crs_sf

	crs_ten    = sf*fkin*s_ten
	crs_sf_ten = sf*fkin*a_ten
	crs_ten    =  crs_sf_ten
	
	
	do kjd = -1,1,1
	crs_sf_pol(kjd) = sf*fkin*a_pol(kjd)
	crs_pol(kjd) = crs_sf_pol(kjd)
	enddo
*	print *,crs_pol
*	print *,a_pol
	
	
*	print *,v_l,v_t,v_tl,v_tt
*	print *,w_l,w_t,w_tl,w_tt
        return
        end


       subroutine j_d(q0,qv,p_f,ktf,sp_f,p_r,ktr,jd,icon,iv1,iv2,noff)

**********************************************************************
* j_d - v.3.00
* Electromagnetic Current of the Deuteron 
*
* 
* mu  -0,1,2,3, component of the current 
* q0     - virtual photon energy
* qv     - virtual photon momentum
* p_f(3) - three vector of knocked-out nucleon 
* ktf    - isospin of the knocked-out nucleon 1-proton, -1 neutron
* ksf    - spin of the knocked-out nucleon 1-up, -1 down
* sp_f(3)- polarization of the knocked out nucleon
* p_r(3) - three vector of recoil nucleon
* ktr    - isospin of recoil nucleon 1-proton -1 neutron 
* ksr    - spin of recoil nucleon 1-up -1 down 
* kjd    - total spin projectsion of d2 target 1, 0, -1             
* jd(mu,kjd,kf,ks)  - the value of the deuteron electromagnetic current components mu=0,1,2,3- complex
*          kjd = -1,-0,1, kf = 1,2 ks = 1,2  1- up, 2 - down
* icon   - initialization (0),  PWIA (1), FSI (2), CHEX (3),  PWIA+FSI (12), PWIA+CHEX(13), PWIA+FSI+CHEX (123)
* iv1    - version of diagonal FSI see explanation on the top of the code
* iv2    - version of charge-exchange parameterization, see explanation on the top of the code
*
*
* The azimuthal angles are defined in the reference frame where 
* q || z and xz is the (eq), e is the initial electron
*
*jd is 4d matrix
* 
* Misak Sargsian
* FIU, July, 2005
* FIU, July, 2009
*
********************************************************************** 
	implicit none
	complex, dimension(0:3,-1:1,2,2) :: jd
        common/par/pi,pm,pmp,pmn,dm,eb
	real pi,pm,pmp,pmn,dm,eb
        real, dimension(3) :: p_f,p_r,p_m,q
	real prm,thr,phir,pmv,thm,phim
	real e_on,e_off
*        common/sfun/ wa_re,  wa_im
        common/Bystr/ibys !new flag,when ibys = 1 reads said ampl.from data file
        integer ibys
	complex wf0,wf1,wf 
        complex, dimension(2,2) :: jN
        complex, dimension(0:3,2,2):: jNm,jN_on,jN_off
 	complex, dimension(-1:1) :: A0pol_d       
	real, dimension(3) :: sp_f,sp_r,sp_i
	real, dimension(10) :: BB
*	complex, dimension(-1:1) :: A0pol_d
	real q2,qv,q0
	complex spec,A0,A1_dg,A1_chx
        common/which_wave/iw
	integer mu,ktf,ktr,kjd,ksf,sf,ksr,sr,ksm,s1,iw
	real N
	integer icon,iv1,iv2
	real xxx,fd_paris,fd_v18,fd_bonn,fd_v18sb
	integer noff
        ibys = 1

        jd  = (0.0,0.0)
        
 
 	if(icon.eq.0)then
****************************************************************
*       initialization 
****************************************************************
        pi  = acos(-1.0)       
        pm  = 0.938279         
        dm  = 1.875628
        pmp = pm
        pmn = 0.939565  
	eb  = 0.00221
        if(iw.eq.1)xxx = fd_paris(0.,1)  ! initialize the Paris WF
        if(iw.eq.2)xxx = fd_v18(0.,1)    ! initialize the V18 WF
        if(iw.eq.3)xxx = fd_bonn(0.,1)    ! initialize the V18 WF
        if(iw.eq.4)xxx = fd_v18sb(0.,1)    ! initialize the V18 WF
        call Bystricky_dat(0.01,0.,-1,BB) ! initialization
*	print *,"aha_aha"
        return
        endif
****************************************************************
*       end initialization 
****************************************************************

****************************************************************
*       INPUT      
****************************************************************
	q2   = qv**2 - q0**2
	q(1) = 0.0
	q(2) = 0.0
	q(3) = qv
************************************** 
*       spectator/recoil nucleon
*       angles are  with respect to q
**************************************
  	call polar(p_r(1),p_r(2),p_r(3),prm,thr,phir)
	p_m(1) = - p_r(1)
	p_m(2) = - p_r(2)
	p_m(3) = - p_r(3)
	
	
        pmv  = prm
        thm  = pi - thr
        phim = pi + phir
	e_on  = sqrt(pm**2 + pmv**2)
	e_off = dm - e_on
    
	sp_i(1) = 0.0; sp_i(2) = 0.0; sp_i(3) = 1.0

*       call  j_Nm(ktf,q2,q,p_f,p_m,jNm)
      call j_N_matrix_s_onof(ktf,q2,q,p_f,sp_f,p_m,sp_i,jN_on,jN_off,
     & e_on,e_off,noff)
*	write(24,*)"ktf=",ktf,"q2=",q2,"q=",q,"p_f=",p_f,
*     &    "sp_f=",sp_f,"p_m=",p_m,"sp_i=",sp_i,"j_on=",jN_on

	if(noff.eq.0)jNm = jN_on
	if(noff.gt.0)jNm = jN_on + jN_off

	A0pol_d = (0.0,0.0)
	do mu=0,3      ! running \mu components 
        do kjd = -1,1,1  ! spin of the deuteron
	do ksf = -1,1,2  ! spin of the knocked out nucleon  1 is up 2 is down 
	do ksr = -1,1,2  ! spin of the recoil  nucleon  1 is up 2 is down

***********************************************
*  Spin of the knocked-out nucleon
***********************************************
	             sf  = 1
	if(ksf.eq.-1)sf  = 2
***********************************************
*  Spin of the recoil nucleon
***********************************************
	             sr  = 1
	if(ksr.eq.-1)sr  = 2



       A0      = (0.0,0.0)
       A1_dg   = (0.0,0.0)
       A1_chx  = (0.0,0.0)         
    

      if(icon.eq.1.or.icon.ge.12)then    ! A0  - contribution
 
***********************************************
*  Spin of the initial nucleon
***********************************************
*******************************
      do ksm = -1,1,2      
*******************************
	              s1 = 1
	 if(ksm.eq.-1)s1 = 2
      wf0 = (0.0,0.0)
      call dfunction(kjd,ktf,ksm,ksr,pmv,thm,phim,wf0)	
      A0 = A0 +   jNm(mu,sf,s1)*wf0/N(pmv)
*************************************************************
* calculating PWIA amplitude for given deuteron polarization
*      if(kjd.eq.-1)A0pmn1 =  A0pmn1 + jNm(mu,sf,s1)*wf0/N(pmv)
*       if(kjd.eq.0)A0p0   =  A0p0   + jNm(mu,sf,s1)*wf0/N(pmv)
*       if(kjd.eq.1)A0ppl1 =  A0ppl1  + jNm(mu,sf,s1)*wf0/N(pmv)
*       A0pol_d(kjd) = A0pol_d(kjd) + jNm(mu,sf,s1)*wf0/N(pmv)
      enddo  !ksm
      endif



*	print *, A0
      if(icon.eq.2.or.icon.eq.12.or.icon.eq.123)then    ! A1_dg - contribution
      call A1(mu,q0,qv,ktf,ksf,sp_f,ktr,ksr,p_r,kjd,iv1,noff,A1_dg)
      endif
*	print *,A1_dg,q0,qv,p_r,iv1


      if(icon.eq.3.or.icon.eq.13.or.icon.eq.123)then    ! A1_chex - contribution
*      call      A1(mu,q0,qv,ktf,ksf,ktr,ksr,p_r,kjd,A1_dg)
      call A1_chex(mu,q0,qv,ktf,ksf,sp_f,ktr,ksr,p_r,kjd,iv2,noff,
     &             A1_chx)
      endif

*	write(18,*)A1_dg
      jd(mu,kjd,sf,sr) = A0 + A1_dg + A1_chx
************************************************
*      Applying current conservation
************************************************
*      jd(3) = jd(0)*q0/qv


*      jd(mu) = A0 
      enddo  ! do ksr
      enddo  ! do ksf
      enddo	 ! do kjd
      enddo  ! do mu
      return
      end

	function N(pr)
	implicit none
	common/par/pi,pm,pmp,pmn,dm,eb
	real pi,pm,pmp,pmn,dm,eb
	real er,pr
	real N
	N = 1.0
	er = sqrt(pm**2 + pr**2)
	if(er.ge.dm)return
	N = sqrt(2.0*(dm - er)/dm)
	N = 1.0
	return
	end

      subroutine A1(mu,q0,qv,ktf,ksf,spf,ktr,ksr,p_r,kjd,ivr,nf,a1_fsi)
******************************************************************
*  The subroutine calculates the diagonal FSI
* inputs 
*  mu    - \pu 0,1,2,3
* q0     - virtual photon energy
* qv     - virtual photon momentum
* ktf    - isospin projection of knocked-out nucleon 1 - proton, -1 neutron
* ksf    - spin projection of knocked-out nucleon 1 -> (1/2)  -1 ->(-1/2)
* ktr    - isospin  projection of recoil nucleon 1 - proton, -1 neutron
* ksr    - spin projection of recoil nucleon   1 -> (1/2)  -1 ->(-1/2)
* p_r(3) - x,y,z components of recoil nucleon momentum
* kjd    - spin projection of d2 -1,0,1 
* ivr    - version
*
* wf1_re
* wf1_im
*******************************************************************
	implicit none
	integer mu,ktf,ksf,ktr,ksr,kjd,ivr,nf
	real q0,qv
        external phl, phu, un_1_dg, un_1_dg_ndg
	real     phl, phu, un_1_dg, un_1_dg_ndg
        real, dimension(3) :: p_r,spf
        common/par/pi,pm,pmp,pmn,dm,eb
	real       pi,pm,pmp,pmn,dm,eb
	common/vphoton/Q2,q0s,qvs
	real           Q2,q0s,qvs
	common/knocked/  p_f_x,  p_f_y,  p_f_z
	real             p_f_x,  p_f_y,  p_f_z

        common/missing/  p_m_x,  p_m_y,  p_m_z
	real             p_m_x,  p_m_y,  p_m_z
        common/recoil_r/ p_r_x,  p_r_y,  p_r_z
        real             p_r_x,  p_r_y,  p_r_z

        common/tmomentums/    q(3),p_rr(3),p_f(3),sp_f(3)
	real                  q,   p_rr,   p_f,   sp_f
        common/conf/itf,isf,itm,ism,itr,isr,ij
	integer     itf,isf,itm,ism,itr,isr,ij
        common/enginv/s
	real          s
        common/delta0/delta_0
	real          delta_0
	common/ireim/ireim
        integer      ireim
	common/comp/nu
	integer     nu
	common/fsi_version/kvr
	integer            kvr
	common/noff/noff
	integer     noff
	complex a1_fsi
	real qa,qb,sum_re,sum_im,eps
	real T_r,E_r,p_rm,p_fm,E_f

	itf    = ktf
	isf    = ksf
        itm    = ktf
        ism    = ksf
        itr    = ktr
        isr    = ksr
        ij     = kjd
	nu     = mu
	kvr    = ivr

	noff   = nf
	a1_fsi    = (0.0,0.0)
*

	q2  = qv**2 - q0**2
	q0s = q0
	qvs = qv

        p_r_x =  p_r(1)
        p_r_y =  p_r(2)
        p_r_z =  p_r(3)

        p_m_x =  -p_r_x
        p_m_y =  -p_r_y
        p_m_z =  -p_r_z

        p_f_x = p_m_x
        p_f_y = p_m_y
        p_f_z = p_m_z + qv

*********************************************************
*        passing three momentums intro the common space
*********************************************************
	p_rr = p_r
	q(1) = 0.0; q(2) = 0.0;	q(3) = qv
	p_f  = - p_r + q
	sp_f = spf

********************************************************************
*         calculating deltra factor
********************************************************************	
        p_fm = sqrt(p_f_x**2 + p_f_y**2 + p_f_z**2)
        p_rm = sqrt(p_r_x**2 + p_r_y**2 + p_r_z**2)

        E_f = sqrt(pm**2+p_fm**2)
        E_r = sqrt(pm**2+p_rm**2)
	
        T_r = E_r - pm
*        delta_0 =  T_r*(dm+q0)/qv !( *(q0/qv*(T_r + eb)
       delta_0 =  q0/qv*(T_r + eb) + dm/qv * T_r

**************************************************
* This block calculates the s of the final NN system
**************************************************
        s = pmp**2+pmn**2+2.0*(E_f*E_r - 
     &                p_f_x*p_r_x -p_f_y*p_r_y -p_f_z*p_r_z)
***************************************************************

	ireim = 1
        eps = 0.0001
        qa  = 0.0
        qb  = 1.0
        sum_re = 0.0
	if(ivr.le.2)call gadap2(qa,qb,phl,phu,un_1_dg,eps,sum_re) 
	if(ivr.eq.3)call gadap2(qa,qb,phl,phu,un_1_dg_ndg,eps,sum_re) 

	ireim = -1
        eps = 0.0001
        qa  = 0.0
        qb  = 1.0
        sum_im = 0.0
  	if(ivr.le.2)call gadap2(qa,qb,phl,phu,un_1_dg,eps,sum_im)        
	if(ivr.eq.3)call gadap2(qa,qb,phl,phu,un_1_dg_ndg,eps,sum_im)
        a1_fsi = cmplx(sum_re,sum_im)

*        print *, sum_re,sum_im,ivr
        return
        end

        function phl(q)
	phl = 0.0
	return
	end
	function phu(q)
	phu = 2.0*acos(-1.0)
	return
	end

	function  un_1_dg(qq,phiq)
	implicit none
	real un_1_dg,qq,phiq
	external under_prz,under_prz2
	real     under_prz,under_prz2
 
        common/par/pi,pm,pmp,pmn,dm,eb
        real       pi,pm,pmp,pmn,dm,eb

	common/vphoton/Q2,q0,qv
	real           Q2,q0,qv

	common/knocked/  p_f_x,  p_f_y,  p_f_z
	real             p_f_x,  p_f_y,  p_f_z

        common/missing/  p_m_x,  p_m_y,  p_m_z
        real             p_m_x,  p_m_y,  p_m_z

        common/recoil_r/ p_r_x,  p_r_y,  p_r_z
        real             p_r_x,  p_r_y,  p_r_z

        common/tmomentums/ qph(3),p_r(3),p_f(3),sp_f(3)
        real               qph   ,p_r   ,p_f   ,sp_f

	common/tilda_primes/ p_rpr_til(3), p_fpr_til(3)
	real                 p_rpr_til   , p_fpr_til

        common/conf/ktf,ksf,ktm,ksm,ktr,ksr,kjd
	integer     ktf,ksf,ktm,ksm,ktr,ksr,kjd
        common/enginv/s
	real          s
        common/delta0/delta_0
	real          delta_0
	common/ireim/ireim
	integer      ireim
	common/comp/mu
	integer     mu
	common/fsi_version/ivr
	integer            ivr
	common/noff/noff
	integer     noff
	common/iver_pv/iver_pv
	integer        iver_pv
	complex im,bb_on
	complex un_1_pole
        real    un_1_pv
        real, dimension(3):: p_fpr,p_rpr,p_1,sp_1
	integer sf,s1
	real q_x,q_y
	real p_1_x,p_1_y,p_1_z,p1t,th1t,phi1t
        real epsu,prz_min,prz_max
	real virt,sum
	integer moff
	integer npv
	common/npv/npv

	im      = (0.0,1.0)

        un_1_dg = 0.0

	q_x      = qq*cos(phiq)
	q_y      = qq*sin(phiq)

	p_1_y   = p_m_y + q_y
	p_1_x   = p_m_x + q_x

  	call polar(p_1_x,p_1_y,p_m_z,p1t,th1t,phi1t)

	virt = ((dm - sqrt(pm**2+p1t**2))**2-p1t**2 - pm**2)/2.0/qv

 	p_1_z   = p_m_z + delta_0 !+ virt

*******************************************************************************
*         defining prime momenta
*******************************************************************************
        p_fpr(1) =  p_1_x; p_fpr(2) =  p_1_y; p_fpr(3) =  p_1_z + qv
	p_rpr(1) = -p_1_x; p_rpr(2) = -p_1_y; p_rpr(3) = -p_1_z

	p_fpr_til = p_fpr
	p_rpr_til = p_rpr
**************************************************
* calling B function for the on-shell amplitude
**************************************************
	un_1_pole = (0.0,0.0)
	moff = 0
	call bb_dg(mu,kjd,q2,s,qph,ktf,p_f,ksf,sp_f,p_r,ksr,
     &           p_fpr,p_rpr,ivr,noff,moff,bb_on)


	un_1_pole = im/4.0*bb_on


************************************************************************
* principal value calculation
* if ireim=1 it calculates the real part
* if ireim=-1 it calculates the imaginary part
************************************************************************	
	un_1_pv = 0.0
	sum = 0.0
	iver_pv = 1
	epsu = 0.01
	prz_min = -1.0
	prz_max =  1.0
*o	call gadapu(prz_min,prz_max,under_prz,epsu,sum)


	sum = 0.0
	if(npv.eq.1)then
	epsu = 0.01
	prz_min = -1.0
	prz_max = 1.0
	call gadapu(prz_min,prz_max,under_prz2,epsu,sum)
	endif

	un_1_pv   = -1.0/2.0 * sum 

*	print *, "un_1_pv",un_1_pv,ireim,un_1_pole
	if(ireim.eq.1)then
        un_1_dg = qq*(real(un_1_pole) + un_1_pv)/(2.0*pi)**2 
	else
	un_1_dg = qq*(imag(un_1_pole) + un_1_pv)/(2.0*pi)**2 
	endif
        return
	end



	function under_prz2(p_rpr_z)
	implicit none
	real under_prz2, p_rpr_z

        common/par/pi,pm,pmp,pmn,dm,eb
        real       pi,pm,pmp,pmn,dm,eb

	common/vphoton/Q2,q0,qv
	real           Q2,q0,qv

        common/tmomentums/ qph(3),p_r(3),p_f(3),sp_f(3)
        real               qph   ,p_r   ,p_f   ,sp_f

	common/tilda_primes/ p_rpr_til(3), p_fpr_til(3)
        real                 p_rpr_til   , p_fpr_til

        common/conf/ktf,ksf,ktm,ksm,ktr,ksr,kjd
        integer     ktf,ksf,ktm,ksm,ktr,ksr,kjd


        common/enginv/s
        real          s
	common/ireim/ireim
	integer      ireim

	common/comp/mu
	integer     mu

	common/fsi_version/ivr
	integer            ivr
	common/chex_version/ivr_chx
	integer             ivr_chx
	common/noff/noff
	integer     noff
	common/iver_pv/iver_pv
	integer        iver_pv
	real, dimension(3) :: p_rpr ,p_fpr
	
	complex bb,hama
	real sing,pokr
	real hayt
	integer moff

	under_prz2 = 0.0
	pokr = 0.0001
*	sing = p_rpr_z**2 - p_rpr_til(3)**2
	hayt = p_rpr_z - p_rpr_til(3)
	moff = 1
	if(abs(hayt).le.pokr)return
	p_rpr(1)=p_rpr_til(1); p_rpr(2)=p_rpr_til(2); p_rpr(3)=p_rpr_z
	p_fpr   = -p_rpr + qph

	if(iver_pv.eq.1.or.iver_pv.eq.5)then     ! diagonal FSI
	 call bb_dg(mu,kjd,q2,s,qph,ktf,p_f,ksf,sp_f,p_r,ksr,
     &           p_fpr,p_rpr,ivr,noff,moff,bb)
	elseif(iver_pv.eq.2)then   ! nondiagonal FSI
	 call bb_ndg(mu,kjd,q2,s,qph,ktf,p_f,ksf,sp_f,p_r,ksr,
     &           p_fpr,p_rpr,noff,moff,bb)
*	 call  bb_dg(mu,kjd,q2,s,qph,ktf,p_f,ksf,sp_f,p_r,ksr,
*     &           p_fpr,p_rpr,2,noff,moff,bb)
	elseif(iver_pv.eq.3)then   ! chex
	 call bb_chex(mu,kjd,q2,s,qph,ktf,p_f,ksf,sp_f,p_r,ksr,
     &           p_fpr,p_rpr,ivr_chx,noff,moff,bb)
	endif
*	hama = (p_rpr_z + p_rpr_til(3))*bb_p + 
*     &        (-p_rpr_z + p_rpr_til(3))*bb_m
	hama = bb
	if(ireim.eq.1)then
	under_prz2 = real(hama)/hayt/(2.0*pi)
	else
	under_prz2 = imag(hama)/hayt/(2.0*pi)
	endif
*	print *,"under_prz",under_prz
	return
	end



	function under_prz(p_rpr_z)
        common/par/pi,pm,pmp,pmn,dm,eb
	common/vphoton/Q2,q0,qv
        common/tmomentums/ qph(3),p_r(3),p_f(3),sp_f(3)
	common/tilda_primes/ p_rpr_til(3), p_fpr_til(3)
        common/conf/ktf,ksf,ktm,ksm,ktr,ksr,kjd
        common/enginv/s
	common/ireim/ireim
	common/comp/mu
	common/fsi_version/ivr
	common/noff/noff
	common/iver_pv/iver_pv
	real, dimension(3) :: p_rpr ,p_fpr
	complex bb,bb_til,hama
	real p_rpr_z
	real sing,pokr
	under_prz = 0.0
	moff = 1
	pokr = 0.0001
	sing = p_rpr_z**2 - p_rpr_til(3)**2
	if(sing.le.pokr)return
	p_rpr(1)=p_rpr_til(1); p_rpr(2)=p_rpr_til(2); p_rpr(3)=p_rpr_z
	p_fpr   = -p_rpr + qph

	if(iver_pv.eq.1)then     ! diagonal FSI
	 call bb_dg(mu,kjd,q2,s,qph,ktf,p_f,ksf,sp_f,p_r,ksr,
     &           p_fpr,p_rpr,ivr,noff,moff,bb)

	 call bb_dg(mu,kjd,q2,s,qph,ktf,p_f,ksf,sp_f,p_r,ksr,
     &           p_fpr_til,p_rpr_til,ivr,noff,moff,bb_til)
*	 print *,"aha"
	elseif(iver_pv.eq.2)then   ! nondiagonal FSI
	 call bb_ndg(mu,kjd,q2,s,qph,ktf,p_f,ksf,sp_f,p_r,ksr,
     &           p_fpr,p_rpr,noff,moff,bb)

	 call bb_ndg(mu,kjd,q2,s,qph,ktf,p_f,ksf,sp_f,p_r,ksr,
     &           p_fpr_til,p_rpr_til,noff,moff,bb_til)
*	 print *,"aha2"

	elseif(iver_pv.eq.3)then   ! chex
	 call bb_chex(mu,kjd,q2,s,qph,ktf,p_f,ksf,sp_f,p_r,ksr,
     &           p_fpr,p_rpr,ivr,noff,moff,bb)

	 call bb_chex(mu,kjd,q2,s,qph,ktf,p_f,ksf,sp_f,p_r,ksr,
     &           p_fpr_til,p_rpr_til,ivr,noff,moff,bb_til)
	endif
	hama = (p_rpr_z + p_rpr_til(3))*bb -2.0*p_rpr_til(3)*bb_til
	if(ireim.eq.1)then
	under_prz = real(hama)/sing/(2.0*pi)
	else
	under_prz = imag(hama)/sing/(2.0*pi)
	endif
*	print *,"under_prz",under_prz
	return
	end


	function  un_1_dg_ndg(qq,phiq)
	implicit none
	real      un_1_dg_ndg,qq,phiq
	external under_prz2
	real     under_prz2
        common/par/pi,pm,pmp,pmn,dm,eb
        real       pi,pm,pmp,pmn,dm,eb
	common/vphoton/Q2,q0,qv
	real           Q2,q0,qv
	common/knocked/  p_f_x,  p_f_y,  p_f_z
	real             p_f_x,  p_f_y,  p_f_z
        common/missing/  p_m_x,  p_m_y,  p_m_z
        real              p_m_x,  p_m_y,  p_m_z
        common/recoil_r/ p_r_x,  p_r_y,  p_r_z
        real             p_r_x,  p_r_y,  p_r_z
        common/tmomentums/ qph(3),p_r(3),p_f(3),sp_f(3)
        real               qph   ,p_r   ,p_f   ,sp_f
	common/tilda_primes/ p_rpr_til(3), p_fpr_til(3)
	real                 p_rpr_til   , p_fpr_til
        common/conf/ktf,ksf,ktm,ksm,ktr,ksr,kjd
        integer     ktf,ksf,ktm,ksm,ktr,ksr,kjd
        common/enginv/s
	real          s
        common/delta0/delta_0
	real          delta_0
	common/ireim/ireim
        integer      ireim
	common/comp/mu
	integer     mu
	common/noff/noff
	integer     noff
	common/iver_pv/iver_pv
	integer        iver_pv
	complex im,bb_on
	complex un_1_pole
        real    un_1_pv
        real, dimension(3):: p_fpr,p_rpr,p_1,sp_1
	integer sf,s1,sfpr
	real q_x,q_y
	real p_1_x,p_1_y,p_1_z,p1t,th1t,phi1t
	real virt,sum,prz_min,prz_max,epsu
	integer moff
	integer npv
	common/npv/npv

	im      = (0.0,1.0)
        un_1_dg_ndg = 0.0
	q_x      = qq*cos(phiq)
	q_y      = qq*sin(phiq)

	p_1_y   = p_m_y + q_y
	p_1_x   = p_m_x + q_x

  	call polar(p_1_x,p_1_y,p_m_z,p1t,th1t,phi1t)

	virt = ((dm - sqrt(pm**2+p1t**2))**2-p1t**2 - pm**2)/2.0/qv

 	p_1_z   = p_m_z + delta_0 !+ virt

*******************************************************************************
*         defining prime momenta
*******************************************************************************
        p_fpr(1) =  p_1_x; p_fpr(2) =  p_1_y; p_fpr(3) =  p_1_z + qv
	p_rpr(1) = -p_1_x; p_rpr(2) = -p_1_y; p_rpr(3) = -p_1_z

	p_fpr_til = p_fpr
	p_rpr_til = p_rpr

	
**************************************************
* calling B function for the on-shell amplitude
**************************************************
	un_1_pole = (0.0,0.0)
	moff = 0
	call bb_ndg(mu,kjd,q2,s,qph,ktf,p_f,ksf,sp_f,p_r,ksr,
     &           p_fpr,p_rpr,noff,moff,bb_on)


	un_1_pole = im/4.0*bb_on

************************************************************************
* principal value calculation
* if ireim=1 it calculates the real part
* if ireim=-1 it calculates the imaginary part
************************************************************************	
	sum = 0.0
	if(npv.eq.1)then
	un_1_pv = 0.0
	iver_pv = 2
	epsu = 0.01
	prz_min = -1.0
	prz_max =  1.0
	call gadapu(prz_min,prz_max,under_prz2,epsu,sum)
	endif

	un_1_pv   = -1.0/2.0 * sum 
*	print *,sum

	if(ireim.eq.1)then
        un_1_dg_ndg = qq*(real(un_1_pole) + un_1_pv)/(2.0*pi)**2 
	else
	un_1_dg_ndg = qq*(imag(un_1_pole) + un_1_pv)/(2.0*pi)**2 
	endif
        return
	end



      subroutine A1_chex(mu,q0,qv,ktf,ksf,spf,ktr,ksr,p_r,kjd,ivr,nf,
     &                   a1chex)
******************************************************************
*  The subroutine calculates the charge-exchange FSI
* inputs 
*  mu    - \pu 0,1,2,3
* q0     - virtual photon energy
* qv     - virtual photon momentum
* ktf    - isospin projection of knocked-out nucleon 1 - proton, -1 neutron
* ksf    - spin projection of knocked-out nucleon 1 -> (1/2)  -1 ->(-1/2)
* ktr    - isospin  projection of recoil nucleon 1 - proton, -1 neutron
* ksr    - spin projection of recoil nucleon   1 -> (1/2)  -1 ->(-1/2)
* p_r(3) - x,y,z components of recoil nucleon momentum
* kjd    - spin projection of d2 -1,0,1 
* ivr    - CHEX version 21, 10, 1021, -21, 10021
* wf1_re
* wf1_im
*******************************************************************
	implicit none
	integer mu,ktf,ksf,ktr,ksr,kjd,ivr,nf
	real q0,qv
        external phl, phu, un_1_chex
	real     phl, phu, un_1_chex
        real, dimension(3):: p_r,spf
        common/par/pi,pm,pmp,pmn,dm,eb
        real       pi,pm,pmp,pmn,dm,eb

	common/vphoton/Q2,q0s,qvs
	real           Q2,q0s,qvs

	common/knocked/  p_f_x,  p_f_y,  p_f_z
	real             p_f_x,  p_f_y,  p_f_z

        common/missing/  p_m_x,  p_m_y,  p_m_z
        real             p_m_x,  p_m_y,  p_m_z



        common/recoil_r/ p_r_x,  p_r_y,  p_r_z
        real             p_r_x,  p_r_y,  p_r_z



        common/tmomentums/ q(3),p_rr(3),p_f(3),sp_f(3)
	real               q   ,p_rr   ,p_f   ,sp_f

        common/conf/itf,isf,itm,ism,itr,isr,ij
	integer     itf,isf,itm,ism,itr,isr,ij


        common/enginv/s
	real          s

        common/delta0/delta_0
	real          delta_0

	common/ireim/ireim
	integer      ireim
	common/comp/nu
	integer     nu
	common/chex_version/kvr
	integer             kvr
	common/noff/noff
	integer     noff
	complex a1chex
	real p_fm,E_f
	real p_rm,E_r,T_r

	real qa,qb,eps,sum_re,sum_im


	itf    = ktf
	isf    = ksf
        itm    = ktr
        ism    = ksr
        itr    = ktr
        isr    = ksr
        ij     = kjd
	nu     = mu
	kvr    = ivr

	noff   = nf
	a1chex    = (0.0,0.0)
*

	q2  = qv**2 - q0**2
	q0s = q0
	qvs = qv

        p_r_x =  p_r(1)
        p_r_y =  p_r(2)
        p_r_z =  p_r(3)

        p_m_x =  -p_r_x
        p_m_y =  -p_r_y
        p_m_z =  -p_r_z

        p_f_x = p_m_x
        p_f_y = p_m_y
        p_f_z = p_m_z + qv

*********************************************************
*        passing three momentums intro the common space
*********************************************************
	p_rr = p_r
	q(1) = 0.0; q(2) = 0.0;	q(3) = qv
	p_f  = - p_r + q
	sp_f = spf


********************************************************************
*         calculating delta factor
********************************************************************	

        p_fm = sqrt(p_f_x**2 + p_f_y**2 + p_f_z**2)
        p_rm = sqrt(p_r_x**2 + p_r_y**2 + p_r_z**2)
        E_f = sqrt(pm**2+p_fm**2)
        E_r = sqrt(pm**2+p_rm**2)
	

        T_r = E_r - pm
*        delta_0 =  T_r*(dm+q0)/qv !( *(q0/qv*(T_r + eb)
       delta_0 =  q0/qv*(T_r + eb) + dm/qv * T_r

**************************************************
* This block calculates the s of the final NN system
**************************************************
        s = pmp**2+pmn**2+2.0*(E_f*E_r - 
     &                p_f_x*p_r_x -p_f_y*p_r_y -p_f_z*p_r_z)
***************************************************************

	ireim = 1
        eps = 0.0001
        qa  = 0.0
        qb  = 1.0
        sum_re = 0.0
        call gadap2(qa,qb,phl,phu,un_1_chex,eps,sum_re) 
	ireim = -1
        eps = 0.0001
        qa  = 0.0
        qb  = 1.0
        sum_im = 0.0
  	call gadap2(qa,qb,phl,phu,un_1_chex,eps,sum_im)                

        a1chex  = cmplx(sum_re,sum_im)

*        write(12,*)sum_re,sum_im
        return
        end

	function  un_1_chex(qq,phiq)
	implicit none
	real      un_1_chex,qq,phiq
	external under_prz2
	real     under_prz2
        common/par/pi,pm,pmp,pmn,dm,eb
        real       pi,pm,pmp,pmn,dm,eb
	common/vphoton/Q2,q0,qv
	real           Q2,q0,qv
	common/knocked/  p_f_x,  p_f_y,  p_f_z
	real             p_f_x,  p_f_y,  p_f_z
        common/missing/  p_m_x,  p_m_y,  p_m_z
        real             p_m_x,  p_m_y,  p_m_z
        common/recoil_r/ p_r_x,  p_r_y,  p_r_z
        real             p_r_x,  p_r_y,  p_r_z
        common/tmomentums/ qph(3),p_r(3),p_f(3),sp_f(3)
	real               qph   ,p_r   ,p_f   ,sp_f
	common/tilda_primes/ p_rpr_til(3), p_fpr_til(3)
	real                 p_rpr_til   , p_fpr_til
        common/conf/ktf,ksf,ktm,ksm,ktr,ksr,kjd
	integer     ktf,ksf,ktm,ksm,ktr,ksr,kjd
        common/enginv/s
	real          s
        common/delta0/delta_0
	real          delta_0
	common/ireim/ireim
	integer      ireim
	common/comp/mu
	integer     mu
	common/chex_version/ivr
	integer             ivr
	common/noff/noff
	integer     noff
	common/iver_pv/iver_pv
	integer        iver_pv
	complex im,bb_on
	complex un_1_pole
	real    un_1_pv
        real, dimension(3):: p_fpr,p_rpr,p_1,sp_1
	integer sf,s1,sfpr
	real q_x,q_y
	real p_1_x,p_1_y,p_1_z,p1t,th1t,phi1t
	real virt
	real prz_min,prz_max,epsu,sum
	integer moff
	integer npv
	common/npv/npv

	im      = (0.0,1.0)

        un_1_chex = 0.0

	q_x      = qq*cos(phiq)
	q_y      = qq*sin(phiq)

	p_1_y   = p_m_y + q_y
	p_1_x   = p_m_x + q_x

  	call polar(p_1_x,p_1_y,p_m_z,p1t,th1t,phi1t)

	virt = ((dm - sqrt(pm**2+p1t**2))**2-p1t**2 - pm**2)/2.0/qv

 	p_1_z   = p_m_z + delta_0 !+ virt

*******************************************************************************
*         defining prime momenta
*******************************************************************************
        p_fpr(1) =  p_1_x; p_fpr(2) =  p_1_y; p_fpr(3) =  p_1_z + qv
	p_rpr(1) = -p_1_x; p_rpr(2) = -p_1_y; p_rpr(3) = -p_1_z

	p_fpr_til = p_fpr
	p_rpr_til = p_rpr

**************************************************
* calling B function for the on-shell amplitude
**************************************************
	un_1_pole = (0.0,0.0)
	moff = 0
	call bb_chex(mu,kjd,q2,s,qph,ktf,p_f,ksf,sp_f,p_r,ksr,
     &           p_fpr,p_rpr,ivr,noff,moff,bb_on)


	un_1_pole = im/4.0*bb_on
	
************************************************************************
* principal value calculation
* if ireim=1 it calculates the real part
* if ireim=-1 it calculates the imaginary part
************************************************************************	
	sum = 0.0
	if(npv.eq.1)then
	un_1_pv = 0.0
	iver_pv = 3
	epsu = 0.01
	prz_min = -1.0
	prz_max =  1.0
	call gadapu(prz_min,prz_max,under_prz2,epsu,sum)
	endif

	un_1_pv   = -1.0/2.0 * sum 


	if(ireim.eq.1)then
        un_1_chex =  qq*(real(un_1_pole) + un_1_pv)/(2.0*pi)**2 
	else
	un_1_chex =  qq*(imag(un_1_pole) + un_1_pv)/(2.0*pi)**2 
	endif
        return
	end



	subroutine kine_get_thrq(q2,q0,qv,pr,thrq,ifl)
	implicit none
	real q2,q0,qv,pr,thrq
	integer ifl
        common/par/pi,pm,pmp,pmn,dm,eb
        real       pi,pm,pmp,pmn,dm,eb
	real e_r,d,cs_thrq

	ifl = 1
	thrq = 0.0
	e_r = sqrt(pm**2 + pr**2)
	d = -q2 + 2.0*q0*dm  + dm**2
	cs_thrq = (2.0*e_r*(q0+dm)-d)/(2.0*pr*qv)
*	print *, "*",cs_thrq,pm,dm,q0
	if(cs_thrq**2.gt.1.0)return
	ifl  = 0
	thrq = acos(cs_thrq)
	return
	end



	subroutine kine_q0(q2,pr,thr,q0,ifl)
	implicit none     
	real               q2,pr,thr,q0
	integer ifl
        common/par/pi,pm,pmp,pmn,dm,eb
        real       pi,pm,pmp,pmn,dm,eb
        common/electron/eil,erl,ue
        real            eil,erl,ue
	integer i_test
	real er,arm_a,arm_b,arm_g
	real A,B,C,x1,x2
	real test

 	ifl   = 1
	er    = sqrt(pm**2 + pr**2)
   	arm_a = -q2 + dm**2 - 2.0*er*dm
	arm_b = 2.0*(dm-er)
	arm_g = -2.0*pr*cos(thr)

*	if(arm_g.eq.0.0)then
	if(arm_g**2.le.0.00001)then
	q0    = - arm_a/arm_b
	if(q0.lt.0.0.or.q0.gt.eil)return
	ifl   = 0.0
	return
	endif
	A     = arm_b**2 - arm_g**2
	B     = 2.0*arm_a*arm_b
	C     = arm_a**2 - q2*arm_g**2
	if(A.eq.0.0)then
	q0 = -C/B
	if(q0.lt.0.or.q0.gt.eil)return
 	else
   	call  squar_eq(a,b,c,x1,x2,i_test)
	if(i_test.eq.1)return
*	write(6,*)'xx',thr*180.0/pi,a,b,c,x1,x2
    	q0 = x2
	test  = (arm_a + q0*arm_b)/arm_g  
	if(test.lt.0.0.or.q0.lt.0.0.or.q0.gt.eil)then
     	q0 = x1
	test  = (arm_a + q0*arm_b)/arm_g  
	if(test.lt.0.0.or.q0.lt.0.0.or.q0.gt.eil)return
	endif
	endif
	ifl   = 0	
 	return
 	end



        SUBROUTINE SQUAR_EQ(A,B,C,X1,X2,I_TEST)
	implicit none
	real a,b,c,x1,x2,dskm
	integer i_test
        I_TEST = 0
        DSKM  = B**2 - 4.0*A*C
        IF(DSKM.LT.0.0)THEN
                IF(DSKM.GT.-1.0E-4)THEN
                DSKM = 0.0
                GOTO 111
                ELSE
                I_TEST = 1
                ENDIF
        RETURN
        ENDIF
111     X1 = ( -B + SQRT(DSKM) ) / 2.0 / A
        X2 = ( -B - SQRT(DSKM) ) / 2.0 / A
        RETURN
        END



 	subroutine polar(p_x,p_y,p_z,p,th,phi)
	implicit none
	real p_x,p_y,p_z,p,th,phi,arg,sn_th
        p = sqrt(p_x**2 + p_y**2 + p_z**2)
        th = 0.0
        phi = 0.0
        if(p.eq.0.0)return
	arg = p_z/p
	if(arg.gt. 1.0)arg =  1.0
	if(arg.lt.-1.0)arg = -1.0
	th  = acos(arg)
	sn_th = sin(th)
	phi = 0.0
        if(p_y.eq.0.0.and.p_x.lt.0.0)then
        phi = acos(-1.0)
        return
        elseif(p_y.eq.0.0.and.p_x.gt.0.0)then
        phi = 0.0
        return
        endif
	if(sn_th.ne.0.0)then
	arg = p_x/p/sn_th
	if(arg.gt. 1.0)arg =  1.0
	if(arg.lt.-1.0)arg = -1.0
	phi = acos(arg)
	if((p_y/sn_th).lt.0.0)phi = 2.0*acos(-1.0) - phi
	endif
	return
	end


