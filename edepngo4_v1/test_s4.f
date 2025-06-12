        program test
        dimension sp_f(3)
        real, dimension(-1:1):: crs_pol0,crs_pol12

        pi = acos(-1.0)       
        pm = 0.938279         
        dm = 1.875628 

*     Initialization

	icon = 0
        iw   = 1  ! the case for the deuteron wave function 
                  !1-Paris, 2-V18, 3- CD Bonn
        ics  = 1  ! case for the elastic form factor parameterization
                  ! 1 - SLAC, 2-Kelly, 3-Bodek,Arrington
        call sigma_ed(ei,q2,pr,thr,phir,ktf,sp_f,noff,npv,npro,icon,
     &            iv1,iv2,isd,ics,iw,ipp,crs_pol0,crs,cts_ten,ifl,q0)

* sigma_ed(ei,q2,pr,thr,phir,ktf,sp_f,noff,npv,npro,icon,iv1,iv2,isd,ics,iw,ipp,crs,crs_ten,ifl,q0))
**************************************
*       INPUT   PARAMETERS   
**************************************
        ei   = 11.00         ! initial electron energy
        ktf  = 1             ! the knock-out nucleon is proton(1), neutron(-1)
        noff = 1             ! Off-Shell effects included (1) not included (0)
        npv  = 0             ! only pole term in FSI (0), pole+PV (1) 
        npro = 1             ! it calculates d\sigma/dE'_e d\Omega_e d\Omega_pf
        icon = 0             ! calculate (initialization parameter)
        isd  = 1             ! direct mechanism
        ipp  = 0             ! S and D component in deuteron

* for details on parameters see subroutine: 

        q2 =   4.00
                
 
*********************************************************************
*         fixing knock-out nucleons polarization axis along z
*********************************************************************
         sp_f(1) = 0.0; sp_f(2) = 0.0; sp_f(3) = 1.0


         pr   = 0.1           ! recoil nucleon momentum GeV/C
         thr0 = 30.0  ! recoil nucleon polar angle
         thr  = thr0*pi/180.0

 

  
         icon=1                ! PWIA
        call sigma_ed(ei,q2,pr,thr,phir,ktf,sp_f,noff,npv,npro,icon,
     &         iv1,iv2,isd,ics,iw,ipp,crs_pol0,crs0,crs_ten0,ifl,q0)

 
 	     crs12 = 1.0
 	     crs_ten12 = 1.0
 
         iv1 = 1            ! FSI with diffractive parameterization for f_nn
         iv2 = 10           ! Charge Exchange Gibbs_Luiseu parameterisation
         icon=12!!2!3           ! PWI+FSI !+CHEX
        call sigma_ed(ei,q2,pr,thr,phir,ktf,sp_f,noff,npv,npro,icon,
     &        iv1,iv2,isd,ics,iw,ipp,crs_pol12,crs12,crs_ten12,ifl,q0)

       if(ifl.eq.1)goto 1      ! kinematically forbidden
        x = q2/(2.0*pm*q0)
        ratio = 0
        if(crs0.ne.0.0)ratio=crs12/crs0
 
 	    asym0 =  crs_ten0 /crs0    ! tensor asymmetry for PWIA
	    asym12 = crs_ten12/crs12   ! tensor asymmetry with FSI
	    
*******************************************************************************
* crs_pol0(j)  - cross section for deuteron polarization j= -1,0,1, in PWIA
* crs_pol12(j) - cross section for deuteron polarization j= -1,0,1, in PWIA+FSI
********************************************************************************
* making unpolarized cross section from polarized components in PWIA
*****************************************************************************
	    crs0make = (crs_pol0(-1) + crs_pol0(0) + crs_pol0(1))/3.0 
*****************************************************************************
* making tensor unpolarized cross section from polarized components in PWIA
*****************************************************************************
        crsten0make = 
     &  (crs_pol0(-1) + crs_pol0(1) - 2.0*crs_pol0(0))/3.0

        write(6,*)q2,pr,thr0,x,crs0,crs0make,crs12,ratio,
     &  crs_ten0,crsten0make,
     &  crs_ten12,
     &  asym0,asym12
         
        write(12,*)thr0,pr,ei,q2,x,plab,crs0,crs12,ratio

1       continue
  
         
  	end






