program run_sigma_ed4
  implicit None
  real, dimension(3) :: sp_f
  real, dimension(-1:1):: crs_pol, crs_pol0, crs_pol12


  real:: ei, q2, pr, thr, phir, phir0

  integer :: ktf, noff, npv, npro, icon, iv1, iv2, isd, ics, iw, ipp, ifl
  real :: crs, crs_ten, crs_ten0 ,q0  
  real :: asym0, asym12
  real :: crs0, crs12, crs_ten12, eil, elab, esc, plab
  real :: ratio, s, theta_e, phr0, thr0, tiv, x
  real :: w_l, w_t, w_lt, w_tt
  real :: w_l_pwia, w_t_pwia, w_lt_pwia, w_tt_pwia
  
  character *80 ID, comment
  
  ! electron kinematics 
  common/electron/eil,esc,theta_e
  
  ! response function
  common /responses/ w_l, w_t, w_lt, w_tt
  
  real,parameter ::  pi = acos(-1.0)       
  real,parameter ::  pm = 0.938279         
  real,parameter ::  dm = 1.875628 
  real,parameter ::  dtr = pi/180.

  integer :: icon_pwia

  integer :: nkin, nk
  
  ! open input file      
  open(15, file='sigma_ed.kin', status = 'old')
  open(12, file='sigma_ed4.data')
  open(26, file='sigma_ed4.log')

  ! write header into output file
  write(12,*) "#! thr0[f,0]/ phr0[f,1]/ pr[f,2]/ ei[f,3]/ th_e[f,4]/ q2[f,5]/ x[f,6]/ plab[f,7]/ "// &
       "crs0[f,8]/ crs12[f,9]/  ratio[f,10]/ crs_ten0[f,11]/ crs_ten12[f,12]/  "//&
       "crs_pol0_-1[f,13]/ crs_pol0_0[f,14]/ crs_pol0_1[f,15]/ "// &
       "crs_pol12_-1[f,16]/ crs_pol12_0[f,17]/ crs_pol12_1[f,18]/ "// &
       "asym0[f,19]/ asym12[f,20]/ ID[s,21]/ "// &
       "w_l[f,22]/ w_t[f,23]/ w_lt[f,24]/ w_tt[f,25]/ "// &
       "w_l_pwia[f,26]/ w_t_pwia[f,27]/ w_lt_pwia[f,28]/ w_tt_pwia[f,29]/ "      
  
  ! number of kinematics, wave function and formfactor parameterisation
  ! read comments and initial kin. data line
  read(15,*) comment
  read(15,*) nkin, iw, ics
  
  !     Initialization
  icon = 0
  write (6,*) "init--> using wave function type : ", iw
  !      call sigma_ed(ei,q2,pr,thr,phir,ktf,sp_f,noff,npv,npro,icon,
  !     &     iv1,iv2,isd,ics,iw,ipp,crs,cts_ten,ifl,q0)

  call sigma_ed(ei,q2,pr,thr,phir,ktf,sp_f,noff,npv,npro,icon,iv1,iv2,isd,ics,iw,ipp,crs_pol,crs,crs_ten,ifl,q0)

  write (6,*) "init-->done"
  
  
  !     
  !*************************************
  !     INPUT   PARAMETERS   
  !*************************************
  !      ei   = 11.00              ! initial electron energy
  !      pr   = 0.1                ! recoil nucleon momentum GeV/C
  !      thr  = 71.*pi/180.        ! polar angle of the recoil nucleon
  !      phir = pi                 ! azimuthal angle of the recoil nucleon
  !      ktf  = 1                  ! the knock-out nucleon is proton(1), neutron(-1)
  !      noff = 1                  ! Off-Shell effects included (1) not included (0)
  !      npv  = 0                  ! only pole term in FSI (0), pole+PV (1) 
  !      npro = 1                  ! it calculates d\sigma/dE'_e d\Omega_e d\Omega_pf
  !      icon = 0                  ! 0=calculate (initialization parameter)
  !      isd  = 1                  ! spectator mechanism
  !      ipp  = 0                  ! S and D component in deuteron
  !      iv1 = 1                   ! FSI with diffractive parameterization for f_nn
  !      iv2 = 10                  ! Charge Exchange Gibbs_Luiseu parameterisation
  !      icon=12                   !3           ! PWI+FSI !+CHEX
  
  
  ! loop over kinematics
  read(15,*) comment
  
  do nk = 1, nkin
     read(15,*) ei, pr, thr0, phir0, ktf, noff, npv, npro, icon, isd, ipp, iv1, iv2, q2 , ID
     ! convert to radians         
     thr = thr0*dtr
     phir = phir0*dtr
     
     !*********************************************************************
     !     fixing knock-out nucleons polarization axis along z
     !********************************************************************
     sp_f(1) = 0.0; sp_f(2) = 0.0; sp_f(3) = 1.0
     ! first PWIA calc
     icon_pwia = 1                 ! PWIA
     print *, 'start PWIA calc'
     call sigma_ed(ei,q2,pr,thr,phir,ktf,sp_f,noff,npv,npro,icon_pwia,iv1,iv2,isd,ics,iw,ipp,crs_pol0, crs0, crs_ten0, ifl,q0) 
     ! save pwia response functions
     w_l_pwia = w_l
     w_t_pwia = w_t
     w_lt_pwia = w_lt
     w_tt_pwia = w_tt   
     if (icon .ne. icon_pwia) then  ! do pwia calculation if desired
         print *, 'start FSI calc'
         ! now calcualtions with selected FSI
         call sigma_ed(ei,q2,pr,thr,phir,ktf,sp_f,noff,npv,npro,icon, iv1,iv2,isd,ics,iw,ipp,crs_pol12,crs12,crs_ten12,ifl,q0)
     endif
     
     write (26,*) "Kinematics : ", nk, " theta, phi = ", thr0, phir0, crs_pol0, crs_pol12
     call flush
     
     if(ifl .ne. 1) then
        x = q2/(2.0*pm*q0)
        ratio = 0
        if(crs0.ne.0.0)ratio=crs12/crs0
        
        tiv = 0.389385*1000.*1000.
        s = dm**2 + 2.0*dm*q0 - q2
        elab = (s-2.0*pm**2)/(2.0*pm)
        plab = sqrt(elab**2-pm**2)
        asym0 = crs_ten0 /crs0      ! asymmetry for PWIA
        asym12 = crs_ten12/crs12   ! asymmetry with FSI
        
        write(6,*) q2, pr, thr0, x, plab, crs0, crs12,ratio
        
        write(12,*) thr0, phir0, pr, ei, theta_e/dtr, q2, x, plab, crs0, crs12, ratio, & 
             crs_ten0, crs_ten12, & 
             crs_pol0(-1), crs_pol0(0), crs_pol0(1), &
             crs_pol12(-1), crs_pol12(0), crs_pol12(1), &
             asym0, asym12,  ID, &
             w_l, w_t, w_lt, w_tt, &
             w_l_pwia, w_t_pwia, w_lt_pwia, w_tt_pwia
     else
        write(26, *) "Kinematics : ", nk, " theta, phi = ", thr0, phir0, " is bad !"
     endif
  enddo
  close(12)
  close(15)
  close(26)
end program run_sigma_ed4




