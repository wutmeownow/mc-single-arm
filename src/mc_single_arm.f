	program mc_single_arm

C+______________________________________________________________________________
!
! Monte-Carlo of SHMS spectrometer using uniform illumination.
!   This version uses TRANSPORT right-handed coordinate system.
C-______________________________________________________________________________

	implicit none

	include 'hms/struct_hms.inc'
	include 'shms/struct_shms.inc'
	include 'spectrometers.inc'
	include 'constants.inc'
	include 'hbook.inc'

c Vector (real*4) for hut ntuples - needs to match dimension of variables
	real*8		shms_hut(37)
	real*8          shms_spec(59)

	real*8          hms_hut(37)
c
	real*8 xs_num,ys_num,xc_sieve,yc_sieve
	real*8 xsfr_num,ysfr_num,xc_frsieve,yc_frsieve
        logical use_front_sieve /.false./
        logical use_sieve /.false./           
c
        common /sieve_info/  xs_num,ys_num,xc_sieve,yc_sieve
     > ,xsfr_num,ysfr_num,xc_frsieve,yc_frsieve,use_sieve, use_front_sieve


C Local declarations.
	integer*4	i,
     >			chanin	/1/,
     >			chanout	/2/,
     >			n_trials,trial,
     >			tmp_int

	integer*4 Itrial                        ! TH - add this for gfortran: forces integer type cast
	logical*4	iss

	real*8 th_nsig_max                      ! TH - add this for gfortran
	parameter(th_nsig_max=3.0d0)            !max #/sigma for gaussian ran #s

C Event limits, topdrawer limits, physics quantities
	real*8 gen_lim(8)			!M.C. phase space limits.
	real*8 gen_lim_up(3)
	real*8 gen_lim_down(3)
        real*8 gen_mom                         !local variable for track momentum in elastic event

	real*8 cut_dpp,cut_dth,cut_dph,cut_z	!cuts on reconstructed quantities
	real*8 xoff,yoff,zoff                   !Beam offsets
        real*8 spec_xoff,spec_yoff,spec_zoff    !Spectrometer offsets
	real*8 spec_xpoff, spec_ypoff           !Spectrometer angle offsets
	real*8 th_ev,cos_ev,sin_ev		!cos and sin of event angle

	real*8 cos_ts,sin_ts			!cos and sin of spectrometer angle
	real*8 x,y,z,dpp,dxdz,dydz,t1,t2,t3,t4	!temporaries

	real*8 x_a,y_a,z_a,dydz_a,dif_a,dydz_aa,dif_aa   ! TH - for target aperture check
	real*8 musc_targ_len			!target length for multiple scattering
        real*8 foil_nm, foil_tk                 !multifoil target
	real*8 foil_zcent
        parameter (foil_tk=0.02)
	real*8 m2				!particle mass squared.
	real*8 rad_len_cm			!conversion r.l. to cm for target
	real*8 pathlen				!path length through spectrometer.
	logical*4 ok_spec			!indicates whether event makes it in MC
	integer*4 hit_calo                      !flag for hitting the calorimeter
	integer*4 armSTOP_successes,armSTOP_trials
        real *8 beam_energy, el_energy, theta_sc!elastic calibration
        real *8 tar_mass, tar_atom_num          !elastic calibration
        real *8 ebeam_model            !beam energy (GeV) for F1F2IN21
        real *8 Z_tar                  !target charge for F1F2IN21
		integer sf_model_flag           !0=F1F2IN21(default), 1=3He table fit
		integer sf_fit_imod             !fit variant for GETSF_F1F2fit
        real *8 Eprime, Q2_model, W2_model, nu_model
        real *8 W_model, xbj_model, theta_model
        real *8 F1_model, F2_model, FL_model
		logical SF_STAT
        real *8 Mp_GeV

C --- Added for 3He = 2p + n construction ---
	real*8 F1_p, F2_p
	real*8 F1_n, F2_n
	real*8 Z_p, A_p
	real*8 Z_n, A_n	

C --- Added: cross section (optional) absolute rate, analogous to old script ---
	real*8 p_accept, th_accept, ph_accept
	real*8 mott_nb, tan2, w1_model, w2_inel_model
	real*8 sigma_f1f2, sigma_weight, rate_f1f2
	real*8 th2, p_spec_GeV, targ_len_m, jac_ep
	real*8 n_areal_m2, lumi_per_C
	real*8 beam_current_uA, target_dens_m3
	real*8 alpha_em, gev2_to_nb, nb_to_m2
	real*8 echarge, echarge_uC
	
C=======================================================================
C New variable declarations for inclusive-model phase-space weighting
C=======================================================================
	double precision p0_GeV	! Central spectrometer momentum [GeV/c]
	double precision pprime_GeV ! Event scattered momentum p' [GeV/c]
	double precision jac_E_delta ! dE'/d(delta), delta fractional
	double precision jac_xpy ! dOmega/(dxptar dyptar)
	double precision ddelta_width ! Full generated width in delta
	double precision dxptar_width ! Full generated width in xptar [rad]
	double precision dyptar_width ! Full generated width in yptar [rad]
	
	parameter (alpha_em   = 7.2973525693d-3)      ! fine-structure constant
	parameter (gev2_to_nb = 3.89379338d5) ! 1 GeV^-2 = 3.89379338e5 nb
	parameter (nb_to_m2 = 1.0d-37) ! 1 GeV^-2 = 3.89379338e5 nb
	parameter (echarge    = 1.602176634d-19) ! Coulomb
	
C Initial and reconstructed track quantities.
	real*8 dpp_init,dth_init,dph_init,xtar_init,ytar_init,ztar_init
	real*8 dpp_recon,dth_recon,dph_recon,ztar_recon,ytar_recon
	real*8 x_fp,y_fp,dx_fp,dy_fp		!at focal plane
	real*8 fry,fr1,fr2
	real*8 p_spec,th_spec			!spectrometer setting
	real*8 resmult

C Circular raster
	double precision urad, uphi, rr, phi
	real*8 ax_raster, ay_raster	

C Control flags (from input file)
	integer*4 ispec
	integer*4 p_flag			!particle identification
	logical*4 ms_flag
	logical*4 wcs_flag
	logical*4 store_all

c	common /hutflag/ cer_flag,vac_flag
C Hardwired control flags.
	logical*4 hut_ntuple	/.true./
        logical*4 spec_ntuple   /.false./
	logical*4 decay_flag	/.false./

	real*8	dpp_var(2),dth_var(2),dph_var(2),ztg_var(2)
	integer*8	stime,etime

	character*132	str_line
C Local  spectrometer varibales
	real*8 x_s,y_s,z_s
	real*8 dxdz_s,dydz_s,dpp_s

C Function definitions.

	integer*4	last_char
	logical*4	rd_int,rd_real
	real*8          grnd,gauss1
	INTEGER irnd
	REAL rnd(99)
        integer      itime,ij
        character	timestring*30

        character*80 rawname, filename, hbook_filename
	real*4  secnds,zero

	parameter(zero=0.0)

	integer ivar

	save		!Remember it all!

C ================================ Executable Code =============================

C Initialize
C using SIMC unstructured version
C
C SHMS
	shmsSTOP_trials	= 0
        shmsSTOP_targ_hor	= 0
        shmsSTOP_targ_vert	= 0
        shmsSTOP_targ_oct	= 0
        shmsSTOP_FRONTSLIT_hor	= 0
        shmsSTOP_FRONTSLIT_vert	= 0
	shmsSTOP_HB_in	= 0
        shmsSTOP_HB_men = 0
	shmsSTOP_HB_mex	= 0
	shmsSTOP_HB_out	= 0	
	shmsSTOP_DOWNSLIT	= 0
	shmsSTOP_COLL_hor	= 0
	shmsSTOP_COLL_vert	= 0
	shmsSTOP_COLL_oct	= 0
	shmsSTOP_Q1_in	= 0
	shmsSTOP_Q1_men	= 0
	shmsSTOP_Q1_mid	= 0
	shmsSTOP_Q1_mex	= 0
	shmsSTOP_Q1_out	= 0
	shmsSTOP_Q2_in	= 0
	shmsSTOP_Q2_men	= 0
	shmsSTOP_Q2_mid	= 0
	shmsSTOP_Q2_mex	= 0
	shmsSTOP_Q2_out	= 0
	shmsSTOP_Q3_in	= 0
	shmsSTOP_Q3_men	= 0
	shmsSTOP_Q3_mid	= 0
	shmsSTOP_Q3_mex	= 0
	shmsSTOP_Q3_out	= 0
c	shmsSTOP_Q3_out1	= 0
c	shmsSTOP_Q3_out2	= 0
c	shmsSTOP_Q3_out3	= 0
c	shmsSTOP_Q3_out4	= 0
c	shmsSTOP_Q3_out5	= 0
c	shmsSTOP_Q3_out6	= 0
	shmsSTOP_D1_in	= 0
        shmsSTOP_D1_flr = 0
	shmsSTOP_D1_men = 0
	shmsSTOP_D1_mid1 = 0
	shmsSTOP_D1_mid2 = 0
	shmsSTOP_D1_mid3 = 0
	shmsSTOP_D1_mid4 = 0
	shmsSTOP_D1_mid5 = 0
	shmsSTOP_D1_mid6 = 0
	shmsSTOP_D1_mid7 = 0
	shmsSTOP_D1_mex = 0
	shmsSTOP_D1_out	= 0
	shmsSTOP_BP_in  = 0
	shmsSTOP_BP_out = 0
	shmsSTOP_hut	= 0
	shmsSTOP_dc1	= 0
	shmsSTOP_dc2	= 0
	shmsSTOP_s1	= 0
	shmsSTOP_s2	= 0
	shmsSTOP_s3	= 0
	shmsSTOP_cal	= 0
	shmsSTOP_successes	= 0
	shmsSTOP_id = 0
C HMS
	hSTOP_trials	= 0
	hSTOP_fAper_hor	= 0
	hSTOP_fAper_vert= 0
	hSTOP_fAper_oct	= 0
	hSTOP_bAper_hor	= 0
	hSTOP_bAper_vert= 0
	hSTOP_bAper_oct	= 0
	hSTOP_slit	= 0
	hSTOP_Q1_in	= 0
	hSTOP_Q1_mid	= 0
	hSTOP_Q1_out	= 0
	hSTOP_Q2_in	= 0
	hSTOP_Q2_mid	= 0
	hSTOP_Q2_out	= 0
	hSTOP_Q3_in	= 0
	hSTOP_Q3_mid	= 0
	hSTOP_Q3_out	= 0
	hSTOP_D1_in	= 0
	hSTOP_D1_out	= 0
	hSTOP_VP1	= 0
	hSTOP_VP2	= 0
	hSTOP_VP3	= 0
	hSTOP_VP4	= 0
	hSTOP_hut	= 0
	hSTOP_dc1	= 0
	hSTOP_dc2	= 0
	hSTOP_scin	= 0
	hSTOP_cal	= 0
	hSTOP_successes	= 0
	hSTOP_id        = 0

C Open setup file.

	write(*,*)'Enter input filename (assumed to be in infiles dir)'
	read(*,1968) rawname
 1968	format(a)
	filename = '../infiles/'//rawname(1:last_char(rawname))//'.inp'
	print *,filename,'opened'
	open(unit=chanin,status='old',file=filename)

C Define HBOOK/NTUPLE filename if used.
	if (hut_ntuple) then
	  hbook_filename = '../worksim/'//rawname(1:last_char(rawname))//'.bin'
	endif
C Open Output file.
	filename = '../outfiles/'//rawname(1:last_char(rawname))//'.out'
	open (unit=chanout,status='unknown',file=filename)

C Read in real*8's from setup file

	str_line = '!'

C Strip off header

	do while (str_line(1:1).eq.'!')
	  read (chanin,1001) str_line
	enddo

! Read data lines.

	write(*,*),str_line(1:last_char(str_line))
	iss = rd_int(str_line,n_trials)
	if (.not.iss) stop 'ERROR (ntrials) in setup!'

! Spectrometer flag:
	read (chanin,1001) str_line
	write(*,*),str_line(1:last_char(str_line))
	iss = rd_int(str_line,ispec)
	if (.not.iss) stop 'ERROR (Spectrometer selection) in setup!'
! Open HBOOK/NTUPLE file here
!	if(hut_ntuple) then
!	   if(ispec.eq.2) then
!	      call shms_hbook_init(hbook_filename,spec_ntuple)
!	   elseif(ispec.eq.1) then
!	      call hms_hbook_init(hbook_filename,spec_ntuple)
!	   else
!	      write(6,*) 'Uknown spectrometer, stopping.'
!	      stop
!	   endif
!	endif

! Spectrometer momentum:
	read (chanin,1001) str_line
	write(*,*),str_line(1:last_char(str_line))
	iss = rd_real(str_line,p_spec)
	if (.not.iss) stop 'ERROR (Spec momentum) in setup!'

! Spectrometer angle:
	read (chanin,1001) str_line
	write(*,*),str_line(1:last_char(str_line))
	iss = rd_real(str_line,th_spec)
	if (.not.iss) stop 'ERROR (Spec theta) in setup!'
	th_spec = abs(th_spec) / degrad
	cos_ts = cos(th_spec)
	sin_ts = sin(th_spec)

! M.C. limits (half width's for dp,th,ph, full width's for x,y,z)
	do i=1,3
	  read (chanin,1001) str_line
	  write(*,*),str_line(1:last_char(str_line))
	  iss = rd_real(str_line,gen_lim(i))
	  if (.not.iss) stop 'ERROR (M.C. limits) in setup!'
	  gen_lim_down(i) = gen_lim(i)
	  read (chanin,1001) str_line
	  write(*,*),str_line(1:last_char(str_line))
	  iss = rd_real(str_line,gen_lim(i))
	  if (.not.iss) stop 'ERROR (M.C. limits) in setup!'
	  gen_lim_up(i) = gen_lim(i)
	enddo

C Acceptance (constant for the run), analogous to old script
	p_accept  = (gen_lim_up(1) - gen_lim_down(1))/100.d0     ! dp/p (fraction)
	th_accept = (gen_lim_up(2) - gen_lim_down(2))/1000.d0    ! rad (from mrad)
	ph_accept = (gen_lim_up(3) - gen_lim_down(3))/1000.d0    ! rad (from mrad)	

	do i = 4,6
	  read (chanin,1001) str_line
	  write(*,*),str_line(1:last_char(str_line))
	  iss = rd_real(str_line,gen_lim(i))
	  if (.not.iss) stop 'ERROR (M.C. limits) in setup!'
	enddo

! Raster size
	do i=7,8
	   read (chanin,1001) str_line
	   write(*,*),str_line(1:last_char(str_line))
	   iss = rd_real(str_line,gen_lim(i))
	   if (.not.iss) stop 'ERROR (Fast Raster) in setup'
	enddo

! Cuts on reconstructed quantities
	read (chanin,1001) str_line
	write(*,*),str_line(1:last_char(str_line))
	if (.not.rd_real(str_line,cut_dpp)) 
     > stop 'ERROR (CUT_DPP) in setup!'

	read (chanin,1001) str_line
	write(*,*),str_line(1:last_char(str_line))
	if (.not.rd_real(str_line,cut_dth)) 
     > stop 'ERROR (CUT_DTH) in setup!'

	read (chanin,1001) str_line
	write(*,*),str_line(1:last_char(str_line))
	if (.not.rd_real(str_line,cut_dph)) 
     > stop 'ERROR (CUT_DPH) in setup!'

	read (chanin,1001) str_line
	write(*,*),str_line(1:last_char(str_line))
	if (.not.rd_real(str_line,cut_z)) 
     > stop 'ERROR (CUT_Z) in setup!'

! Read in radiation length of target material in cm
	read (chanin,1001) str_line
	write(*,*),str_line(1:last_char(str_line))
	if (.not.rd_real(str_line,rad_len_cm)) 
     > stop 'ERROR (RAD_LEN_CM) in setup!'

! Beam and target offsets
	read (chanin, 1001) str_line
	write(*,*),str_line(1:last_char(str_line))
	iss = rd_real(str_line,xoff)
	if(.not.iss) stop 'ERROR (xoff) in setup!'

	read (chanin, 1001) str_line
	write(*,*),str_line(1:last_char(str_line))
	iss = rd_real(str_line,yoff)
	if(.not.iss) stop 'ERROR (yoff) in setup!'

	read (chanin, 1001) str_line
	write(*,*),str_line(1:last_char(str_line))
	iss = rd_real(str_line,zoff)
	if(.not.iss) stop 'ERROR (zoff) in setup!'

! Spectrometer offsets
	read (chanin, 1001) str_line
	write(*,*),str_line(1:last_char(str_line))
	iss = rd_real(str_line,spec_xoff)
	if(.not.iss) stop 'ERROR (spect. xoff) in setup!'

	read (chanin, 1001) str_line
	write(*,*),str_line(1:last_char(str_line))
	iss = rd_real(str_line,spec_yoff)
	if(.not.iss) stop 'ERROR (spect. yoff) in setup!'

	read (chanin, 1001) str_line
	write(*,*),str_line(1:last_char(str_line))
	iss = rd_real(str_line,spec_zoff)
	if(.not.iss) stop 'ERROR (spect. zoff) in setup!'

	read (chanin, 1001) str_line
	write(*,*),str_line(1:last_char(str_line))
	iss = rd_real(str_line,spec_xpoff)
	if(.not.iss) stop 'ERROR (spect. xpoff) in setup!'

	read (chanin, 1001) str_line
	write(*,*),str_line(1:last_char(str_line))
	iss = rd_real(str_line,spec_ypoff)
	if(.not.iss) stop 'ERROR (spect. ypoff) in setup!'

! read in flag for particle type.
	read (chanin,1001) str_line
	write(*,*),str_line(1:last_char(str_line))
	if (.not.rd_int(str_line,p_flag)) 
     > stop 'ERROR: p_flag in setup file!'


! Read in flag for multiple scattering.
	read (chanin,1001) str_line
	write(*,*),str_line(1:last_char(str_line))
	if (.not.rd_int(str_line,tmp_int)) 
     > stop 'ERROR: ms_flag in setup file!'
	if (tmp_int.eq.1) ms_flag = .true.

! Read in flag for wire chamber smearing.
	read (chanin,1001) str_line
	write(*,*),str_line(1:last_char(str_line))
	if (.not.rd_int(str_line,tmp_int)) 
     > stop 'ERROR: wcs_flag in setup file!'
	if (tmp_int.eq.1) wcs_flag = .true.

! Read in flag to keep all events - success or not
	read (chanin,1001) str_line
	write(*,*),str_line(1:last_char(str_line))
	if (.not.rd_int(str_line,tmp_int)) 
     > stop 'ERROR: store_all in setup file!'
	if (tmp_int.eq.1) store_all = .true.

!     Read in flag for 'beam energy(MeV)' to trigger on elastic event if present
      beam_energy=-0.1  !by default do not use elastic event generator
      tar_atom_num=12.  !by default it is carbon
      ebeam_model=-1.0d0      !default: disable F1F2IN21 unless set
      Z_tar = 1.0d0           !default: proton
	  sf_model_flag = 0        !default model: F1F2IN21
      beam_current_uA = 0.d0  !optional (uA); if <=0, absolute rate disabled
      target_dens_m3  = 0.d0  !optional (#/m^3); if <=0, absolute rate disabled	
      n_areal_m2 = 0.d0
      lumi_per_C = 0.d0
      echarge_uC = echarge * 1.0d6 ! uC

C=======================================================================
C Initialization
C=======================================================================
      p0_GeV       = 0.d0
      pprime_GeV   = 0.d0
      jac_E_delta  = 0.d0
      jac_xpy      = 0.d0
      ddelta_width = 0.d0
      dxptar_width = 0.d0
      dyptar_width = 0.d0
	
      read (chanin,1001,end=1000,err=1000) str_line
      write(*,*),str_line(1:last_char(str_line))
      iss = rd_real(str_line,beam_energy)

	
! Read in flag to use sieve
        read (chanin,1001,end=1000,err=1000) str_line
        write(*,*),str_line(1:last_char(str_line))
        if (.not.rd_int(str_line,tmp_int)) 
     > stop 'ERROR: use_sieve in setup file!'
        if (tmp_int.eq.1) then
          if (ispec.eq.1) use_sieve=.true.
          if (ispec.eq.2) use_sieve=.true.
          if (ispec.eq.2) use_front_sieve=.false.
        endif

!     Read in flag for 'target atomic number (Z+N)' for elastic event if present
      read (chanin,1001,end=1000,err=1000) str_line
      write(*,*),str_line(1:last_char(str_line))
      iss = rd_real(str_line,tar_atom_num)

!     Read in beam energy (GeV) for F1F2IN21 model (optional, <=0 disables call)
      read (chanin,1001,end=1000,err=1000) str_line
      write(*,*),str_line(1:last_char(str_line))
      iss = rd_real(str_line,ebeam_model)

!     Read in target charge Z for F1F2IN21 model (optional)
      read (chanin,1001,end=1000,err=1000) str_line
      write(*,*),str_line(1:last_char(str_line))
      iss = rd_real(str_line,Z_tar)

!     Read in SF model selector (optional): 0=F1F2IN21, 1=3He fit table
      read (chanin,1001,end=1000,err=1000) str_line
      write(*,*) str_line(1:last_char(str_line))
      if (rd_int(str_line,sf_model_flag)) then
         if (sf_model_flag.ne.0 .and. sf_model_flag.ne.1) sf_model_flag=0
      else
         sf_model_flag=0
      endif


!     Optional: beam current (uA) and target number density (#/m^3)
!     Add TWO extra lines at end of your setup file if you want absolute rates.
!
!     Backward compatibility note:
!     some setup files still carry two legacy lines after Z_tar
!     ("Multiple scattering type" and "Enegy Loss"). If present,
!     skip them before reading luminosity inputs.
      read (chanin,1001,end=1000,err=1000) str_line
      write(*,*),str_line(1:last_char(str_line))
      if (index(str_line,'Multiple scattering type').gt.0 .or.
     >    index(str_line,'Enegy Loss').gt.0) then
        read (chanin,1001,end=1000,err=1000) str_line
        write(*,*),str_line(1:last_char(str_line))

        read (chanin,1001,end=1000,err=1000) str_line
        write(*,*),str_line(1:last_char(str_line))
      endif
      iss = rd_real(str_line,beam_current_uA)

      if (iss) then
        read (chanin,1001,end=1000,err=1000) str_line
        write(*,*),str_line(1:last_char(str_line))
        iss = rd_real(str_line,target_dens_m3)
      endif

 1000  continue
      Mp_GeV = 0.93827208d0


	print *, 'ebeam_model=', ebeam_model
	print *, 'beam_energy=', beam_energy
	print *, 'Z_tar=', Z_tar
	print *, 'tar_atom_num=', tar_atom_num
	print *, 'sf_model_flag=', sf_model_flag
	print *, 'beam_current_uA=', beam_current_uA
	print *, 'target_dens_m3=', target_dens_m3

!	Open HBOOK/NTUPLE file here
	if(hut_ntuple) then
		if (sf_model_flag.eq.1) then
		      hbook_filename = hbook_filename(1:last_char(hbook_filename))//
     >                     '3HeFit.bin'
		else
			hbook_filename = hbook_filename(1:last_char(hbook_filename))//'.bin'
		endif
		if(ispec.eq.2) then
			call shms_hbook_init(hbook_filename,spec_ntuple)
		elseif(ispec.eq.1) then
			call hms_hbook_init(hbook_filename,spec_ntuple)
		else
			write(6,*) 'Uknown spectrometer, stopping.'
			stop
		endif
	endif
C	pause

	

C Set particle masses.
	m2 = me2			!default to electron
	if(p_flag.eq.0) then
	  m2 = me2
	else if(p_flag.eq.1) then
	  m2 = mp2
	else if(p_flag.eq.2) then
	  m2 = md2
	else if(p_flag.eq.3) then
	  m2 = mpi2
	else if(p_flag.eq.4) then
	  m2 = mk2
	endif

C------------------------------------------------------------------------------C
C                           Top of Monte-Carlo loop                            C
C------------------------------------------------------------------------------C

	stime = secnds(zero)

! TH - use "Itrial" instead of "trial" for gfortran. Somehow the stringlib.f
! function does not type cast string to integer otherwise.
          itime=time8()
   	  call ctime(itime,timestring)
	  write(6,*) 'Using random seed based on clock time'
          write(6,*) 'Starting random number seed: ',itime
C DJG - If you want to use default (fixed) seed, comment out the line below
          call sgrnd(itime)

	do Itrial = 1,n_trials
	   if(ispec.eq.1) then
	      armSTOP_successes=hSTOP_successes
	   elseif(ispec.eq.2) then
	      armSTOP_successes=shmsSTOP_successes
	   endif
	  if(mod(Itrial,200000).eq.0) write(*,*)'event #: ',
     >Itrial,'       successes: ',armSTOP_successes


	  irnd=Itrial

C Pick starting point within target. Z is picked uniformly, X and Y are
C chosen as truncated gaussians, using numbers picked above.
C Units are cm.

! TH - use a double precision for random number generation here.
	  x = gauss1(th_nsig_max) * gen_lim(4) / 6.0	!beam width
	  y = gauss1(th_nsig_max) * gen_lim(5) / 6.0	!beam height

          if(gen_lim(6).gt.0) then                      
	     z = (grnd() - 0.5) * gen_lim(6)		!along target

          elseif(gen_lim(6).eq.-3) then                 !optics1: three foils
             foil_nm=3*grnd()-1.5                       !20um foils;  z=0, +/- 10cm
             foil_nm=anint(foil_nm)                     != -1, 0, 1
	     foil_zcent = foil_nm * 10
	     z = (grnd() - 0.5) * foil_tk+ foil_nm * 10

          elseif(gen_lim(6).eq.-2) then                 !optics2: two foils
             foil_nm=grnd()                             !20um foils; z= +/- 5cm
             foil_nm=anint(foil_nm)                     != 0, 1
	     foil_zcent = foil_nm * 5
	     z = (grnd() - 0.5) * foil_tk - 5+ foil_nm * 10
	  elseif(gen_lim(6).eq.-5) then
            foil_nm=5*grnd()-2.5                       ! pol target optics
             foil_nm=anint(foil_nm)                     !=
	     if (foil_nm .eq. -2) foil_zcent = 20.
	     if (foil_nm .eq. -1) foil_zcent = 13.34
	     if (foil_nm .eq. 0)  foil_zcent = 0.0
	     if (foil_nm .eq. 1) foil_zcent = -20.
	     if (foil_nm .eq. 2) foil_zcent = -30.
             z= (grnd() - 0.5) * foil_tk + foil_zcent

          endif
C DJG Assume flat raster
c	  fr1 = (grnd() - 0.5) * gen_lim(7)   !raster x
c	  fr2 = (grnd() - 0.5) * gen_lim(8)   !raster y

c	  fry = -fr2  !+y = up, but fry needs to be positive when pointing down

c	  x = x + fr1
c	  y = y + fr2

c	  x = x + xoff
c	  y = y + yoff
c	  z = z + zoff

c=======================================================================
c Generate raster offset as a UNIFORM FILLED DISK / ELLIPSE
c
c Old code:
c   fr1 = (grnd()-0.5d0) * gen_lim(7)
c   fr2 = (grnd()-0.5d0) * gen_lim(8)
c
c That produces a UNIFORM SQUARE / RECTANGLE in (x,y), not a circular
c raster. The replacement below generates points uniformly in area inside
c a unit disk, then scales them to the requested raster full widths.
c
c Input convention is preserved:
c   gen_lim(7) = raster full width in x  (cm)
c   gen_lim(8) = raster full width in y  (cm)
c
c For equal widths:
c   radius = gen_lim(7)/2 = gen_lim(8)/2
c   --> uniform circular raster
c
c For unequal widths:
c   semi-axis x = gen_lim(7)/2
c   semi-axis y = gen_lim(8)/2
c   --> uniform elliptical raster
c
c IMPORTANT:
c   The factor sqrt(rnd) is required for UNIFORM AREA density.
c   Using r = rnd would overpopulate the center.
c=======================================================================

c----- Semi-axes of the raster footprint (cm)
	  ax_raster = 0.5d0 * gen_lim(7)
	  ay_raster = 0.5d0 * gen_lim(8)

c----- Draw two independent random numbers in [0,1)
	  urad = grnd()
	  uphi = grnd()

c----- Convert to polar coordinates for uniform area in a disk:
c      rr  in [0,1] with p(rr) = 2*rr  --> rr = sqrt(U)
c      phi in [0,2*pi)
	  rr  = dsqrt(urad)
	  phi = twopi * uphi

c----- Scale unit disk point to requested raster size
	  fr1 = ax_raster * rr * dcos(phi)
	  fr2 = ay_raster * rr * dsin(phi)

c----- Historical sign convention used elsewhere in the code
	  fry = -fr2

c----- Apply raster offsets to beam position at the target
	  x = x + fr1
	  y = y + fr2
	  
C Pick scattering angles and DPP from independent, uniform distributions.
C dxdz and dydz in HMS TRANSPORT coordinates.

	  dpp  = grnd()*(gen_lim_up(1)-gen_lim_down(1))
     &             + gen_lim_down(1)
	  dydz = grnd()*(gen_lim_up(2)-gen_lim_down(2))
     &          /1000.   + gen_lim_down(2)/1000.
	  dxdz = grnd()*(gen_lim_up(3)-gen_lim_down(3))
     &          /1000.   + gen_lim_down(3)/1000.

C Calculate for the elastic energy calibration using the beam energy.
	  if(beam_energy.gt.0) then
	     if(ispec.eq.2) then ! SHMS
		theta_sc = acos((cos_ts-dydz*sin_ts)/sqrt(1. + dxdz**2. + dydz**2.))
	     elseif(ispec.eq.1) then ! HMS
		theta_sc = acos((cos_ts+dydz*sin_ts)/sqrt(1. + dxdz**2. + dydz**2.))
	     else
		write(6,*) 'Elastic scattering not set up for your spectrometer' 
		STOP
	     endif
	     tar_mass = tar_atom_num*931.5 !carbon
	     el_energy = tar_mass*beam_energy/(tar_mass+2.*beam_energy*(sin(theta_sc/2.))**2)
	     dpp = (el_energy-p_spec)/p_spec*100.
	  endif


	  if(ispec.eq.2) then ! SHMS
C Transform from target to SHMS (TRANSPORT) coordinates.
C Version for a spectrometer on the left-hand side: (i.e. SHMS)
	     x_s    = -y
	     y_s    = x * cos_ts - z * sin_ts
	     z_s    = z * cos_ts + x * sin_ts
	  elseif(ispec.eq.1) then ! HMS
C Below assumes that HMS is on the right-hand side of the beam
C line (looking downstream).
	     x_s    = -y
	     y_s    = x * cos_ts + z * sin_ts
	     z_s    = z * cos_ts - x * sin_ts
	  else
	     write(6,*) 'unknown spectrometer: stopping'
	     stop
	  endif

C DJG Apply spectrometer offsets
C DJG If the spectrometer if too low (positive x offset) a particle
C DJG at "x=0" will appear in the spectrometer to be a little high
C DJG so we should subtract the offset

	  x_s = x_s - spec_xoff
	  y_s = y_s - spec_yoff
	  z_s = z_s - spec_zoff

	  dpp_s  = dpp
	  dxdz_s = dxdz
	  dydz_s = dydz

C DJG Apply spectrometer angle offsets
	  dxdz_s = dxdz_s - spec_xpoff/1000.0
	  dydz_s = dydz_s - spec_ypoff/1000.0

C Drift back to zs = 0, the plane through the target center
	  x_s = x_s - z_s * dxdz_s
	  y_s = y_s - z_s * dydz_s
	  z_s = 0.0

C Save init values for later.
	  xtar_init = x_s
	  ytar_init = y_s
	  ztar_init = z
	  dpp_init = dpp
	  dth_init = dydz_s*1000.		!mr
	  dph_init = dxdz_s*1000.		!mr

C Calculate multiple scattering length of target
	  if(ispec.eq.1) then ! spectrometer on right
	     cos_ev = (cos_ts+dydz_s*sin_ts)/sqrt(1+dydz_s**2+dxdz_s**2)
	  elseif(ispec.eq.2) then ! spectrometer on left
	     cos_ev = (cos_ts-dydz_s*sin_ts)/sqrt(1+dydz_s**2+dxdz_s**2)
	  endif
	  th_ev = acos(cos_ev)
	  sin_ev = sin(th_ev)

C Initialize F1F2IN21-derived quantities for this event
		Eprime       = 0.d0
		Q2_model     = -1.d0
		W2_model     = -1.d0
		nu_model     = -1.d0
		W_model      = -1.d0
		xbj_model    = -1.d0
		theta_model  = -1.d0
		F1_model     = 0.d0
		F2_model     = 0.d0
		FL_model     = 0.d0
		SF_STAT      = .false.
		mott_nb      = 0.d0
		sigma_f1f2   = 0.d0
		sigma_weight = 0.d0
		rate_f1f2    = 0.d0
	  
C Inclusive structure-function model (F1F2IN21) for acceptance weighting
      if (ebeam_model.gt.0.d0) then
         Eprime = p_spec*(1.d0 + dpp_init/100.d0)/1000.d0
	 if(ispec.eq.1) then	! spectrometer on right
	    theta_model = acos((cos_ts+(dth_init/1000.d0)*sin_ts)
     >	               /sqrt(1+(dth_init/1000.d0)**2
     >                 +(dph_init/1000.d0)**2))*degrad	    
	 elseif(ispec.eq.2) then ! spectrometer on left
	    theta_model = acos((cos_ts-(dth_init/1000.d0)*sin_ts)
     >	               /sqrt(1+(dth_init/1000.d0)**2
     >                 +(dph_init/1000.d0)**2))*degrad	    
	 endif	 
         Q2_model = 4.d0*ebeam_model*Eprime
     >	          *(sin((theta_model/degrad)/2.d0)**2)
         nu_model = ebeam_model - Eprime
         W2_model = Mp_GeV*Mp_GeV + 2.d0*Mp_GeV*nu_model - Q2_model
	 
	 if (W2_model.gt.0.d0) then
            W_model = sqrt(W2_model)
         else
            W_model = -1.d0
         endif
         if (nu_model.gt.0.d0) then
            xbj_model = Q2_model/(2.d0*Mp_GeV*nu_model)
         else
            xbj_model = -1.d0
         endif
	else
	   Eprime = 0.d0
	   Q2_model = -1.d0
	   W2_model = -1.d0
	   nu_model = -1.d0
	   W_model = -1.d0
	   xbj_model = -1.d0
	   theta_model = -1.d0
	   F1_model = 0.d0
	   F2_model = 0.d0
	endif
	 
	if (Q2_model.gt.0.d0 .and. W2_model.gt.0.d0 .and.
     >    theta_model.gt.0.d0) then

C        The SF model flag only selects F1_model/F2_model here.
C        Both branches use the common cross-section and weight code below.	 
C        Optional 3He structure-function fit table
         if (sf_model_flag.eq.1 .and. tar_atom_num.eq.3.d0
     >       .and. Z_tar.eq.2.d0) then
            call GETSF_F1F2fit(4,sf_fit_imod,xbj_model,Q2_model,
     >                         F1_model,F2_model,FL_model,SF_STAT)
            if (Itrial.eq.1 .and. SF_STAT) then
               write(6,*) 'SF model=1 uses GETSF_F1F2fit from'
               write(6,*) 'interp/sf_tables/Table_3He_F1F2_SF*.csv'
               write(6,*) 'isf=4 (x,z integrated), imod=', sf_fit_imod
            endif
            if (.not.SF_STAT) then
               write(6,*) 'ERROR: GETSF_F1F2fit failed with sf_model_flag=1'
               write(6,*) 'Check src/interp/sf_tables/Table_3He_F1F2_SF*.csv'
               stop 'SF table lookup failed'
            endif

C       Legacy/default 3He behavior: 2p + n from nucleon F1F2IN21
        else if (tar_atom_num.eq.3.d0 .and. Z_tar.eq.2.d0) then

C           Define proton/neutron A,Z using the same type as model inputs
        	Z_p = 1.d0
            A_p = 1.d0
            Z_n = 0.d0
            A_n = 1.d0

C           Free proton
            call F1F2IN21(Z_p,A_p,Q2_model,W2_model,F1_p,F2_p)

C           Free neutron
            call F1F2IN21(Z_n,A_n,Q2_model,W2_model,F1_n,F2_n)

C           3He = 2p + n
            F1_model = 2.d0*F1_p + F1_n
            F2_model = 2.d0*F2_p + F2_n

        else

C           Default behavior for all other targets
            call F1F2IN21(Z_tar,tar_atom_num,Q2_model,W2_model,
     >                    F1_model,F2_model)

        endif

      else

        F1_p     = 0.d0
        F2_p     = 0.d0
        F1_n     = 0.d0
        F2_n     = 0.d0
        F1_model = 0.d0
        F2_model = 0.d0

      endif
	 

C --- Cross section and (optional) absolute rate (analogous to old block) ---
       mott_nb      = 0.d0
       sigma_f1f2   = 0.d0
       sigma_weight = 0.d0
       rate_f1f2    = 0.d0
	
       if (ebeam_model.gt.0.d0 .and. Q2_model.gt.0.d0 .and.
     >     nu_model.gt.0.d0 .and. theta_model.gt.0.d0) then

          if (xbj_model.gt.0.99d0) then
             sigma_f1f2   = 0.d0
             sigma_weight = 0.d0
             rate_f1f2    = 0.d0
          else
             th2  = (theta_model/degrad)/2.d0
             tan2 = tan(th2)**2

C Mott in nb/sr (GeV^-2 converted to nb via gev2_to_nb)
             mott_nb = ((alpha_em*cos(th2)/
     >           (2.d0*ebeam_model*sin(th2)*sin(th2)))**2)*gev2_to_nb
	             w1_model      = F1_model/Mp_GeV
	             w2_inel_model = F2_model/nu_model
	             sigma_f1f2    = mott_nb*(w2_inel_model
     >                             + 2.d0*w1_model*tan2)

C=======================================================================
C Weight for generation in (delta, xptar, yptar)
C
C Assumptions:
C   sigma_f1f2 = d^2sigma / (dOmega dE')   [nb / sr / GeV]
C   dpp_init   = fractional delta = (p-p0)/p0
C   xptar_init = generated xptar [rad]
C   yptar_init = generated yptar [rad]
C   Eprime     already defined correctly earlier
C=======================================================================
		     p0_GeV     = p_spec/1000.d0
		     pprime_GeV = p0_GeV*(1.d0 + dpp_init/100.d0)

C dE'/d(delta) with delta fractional
		     jac_E_delta = (pprime_GeV/Eprime) * p0_GeV

C Exact solid-angle Jacobian for target slopes
		     jac_xpy = 1.d0 /
     >                         ((1.d0 
     >                         + (dth_init/1000.d0)*(dth_init/1000.d0)
     >                         + (dph_init/1000.d0)*(dph_init/1000.d0)
     >                         )**1.5d0)

C       Use the FULL generated widths here.
C If the generator throws from [-lim,+lim], then full width = 2*lim.
		     ddelta_width = p_accept
		     dxptar_width = th_accept
		     dyptar_width = ph_accept

C Per-trial cross section weight in nb
		     sigma_weight = sigma_f1f2 *
     >                              jac_E_delta *
     >                              jac_xpy *
     >                              ddelta_width *
     >                              dxptar_width * dyptar_width /
     >                              n_trials
C       Optional absolute rate (Hz): requires target_dens_m3 and beam_current_uA
             if (beam_current_uA.gt.0.d0 .and. target_dens_m3.gt.0.d0
     >           .and. gen_lim(6).ne.0.d0) then
C       --- after sigma_weight is computed (still in nb * phase-space) ---
		targ_len_m = abs(gen_lim(6))/100.d0 ! cm -> m
		n_areal_m2 = target_dens_m3 * targ_len_m ! (#/m^3)*m = #/m^2
		lumi_per_C   = n_areal_m2 / echarge_uC ! 1/(m^2*uC)

C       convert nb -> m^2 and normalize by charge contribution
C		sigma_weight = sigma_weight * nb_to_m2 * lumi_per_C
		sigma_weight = sigma_weight * nb_to_m2 ! do luminosity scaling in python
		rate_f1f2 = sigma_weight * target_dens_m3 * nb_to_m2 * lumi_per_C
     >                     * (beam_current_uA*1.d-6)
     >                     * targ_len_m / (echarge)
             endif
          endif
       endif

       if ((Itrial.le.1).or.(mod(Itrial,50000).le.3)) then
          write(*,'("trial #",i8," xsec(nb)=",G14.5," L/C=",G14.5,
     >      " weight=",G14.5," rate(Hz)=",G14.5," at x,Q2=",2F8.4)')
     >      Itrial, sigma_f1f2, lumi_per_C, sigma_weight, rate_f1f2,
     >      xbj_model, Q2_model
       endif
	
C Case 1 : extended cryo target:
C Choices: 
C 1. cryocylinder: Basic cylinder(2.65 inches diameter --> 3.37 cm radius) w/flat exit window (5 mil Al)
C 2. cryotarg2017: Cylinder (1.32 inches radisu)  with curved exit window (same radius) 5 mil sides/exit
C 3. Tuna can: shaped like a tuna can - 4 cm diameter (usually)  - 5 mil window. 
	  if (gen_lim(6).gt.3.) then ! anything longer than 3 cm assumed to be cryotarget
c	    call cryotuna(z,th_ev,rad_len_cm,gen_lim(6),musc_targ_len)
c       call cryocylinder(z,th_ev,rad_len_cm,gen_lim(6),musc_targ_len)
             call he3targ2019(z,th_ev,rad_len_cm, 40.0, musc_targ_len)  !specify 40 cm length all the time	     
c	     call cryotarg2017(z,th_ev,rad_len_cm,gen_lim(6),musc_targ_len)	    
C Simple solid target
	 else
	    if (gen_lim(6).gt.0) then
	    musc_targ_len = abs(gen_lim(6)/2. - z)/rad_len_cm/cos_ev
	    else ! using multifoil target
	    musc_targ_len = abs(foil_tk/2. - (z-foil_zcent))/rad_len_cm/cos_ev	       
	    endif
	 endif

C Scattering before magnets:  Approximate all scattering as occuring AT TARGET.
C SHMS
C  20 mil Al scattering chamber window (X0=8.89cm)
C  57.27 cm air (X0=30420cm)
C spectrometer entrance window
C  10 mil Al s (X0=8.89cm)

C  HMS
C  20 mil Al scattering chamber window (X0=8.89cm)
C  24.61 cm air (X0=30420cm)
C spectrometer entrance window
C  15 mil Kevlar (X0=74.6 cm)
C   5 mil Mylar (X0=28.7 cm)

	  if(ispec.eq.2) then
	     musc_targ_len = musc_targ_len + .020*2.54/8.89 +
     >          57.27/30420. +  .010*2.54/8.89
	  elseif(ispec.eq.1) then
	     musc_targ_len = musc_targ_len + .020*2.54/8.89 +
     >          24.61/30420. +  .015*2.54/74.6 + .005*2.54/28.7
	  endif

c
	  if (ms_flag ) call musc(m2,p_spec*(1.+dpp_s/100.),
     > musc_targ_len,dydz_s,dxdz_s)

!-----------------------------------------------------------------------------
! TH - START TARGET APERTURE TESTS
! ----------------------------------------------------------------------------
! This is for SHMS only
! Restore xs to values at pivot. 
!	   xs = x_transp
!	   ys = y_transp
	  x_a = 0
	  y_a = 2.99 !cm
	  z_a = 57.2 !cm

	  dydz_a = (y_a-ytar_init)/(z_a-ztar_init)
	  dydz_aa = atan(dydz_a)

! Check target aperture, at about 0.572 meter
! theta_a = lower limit of aperture window
! theta_s = scattering angle (=spectrometer angle + position)
! The difference between the scattering and the limiting angle of the
! window for a given central spectrometer angle.
	  dif_a = (th_spec*1000+dth_init-dydz_aa*1000)  ! mrad

! ----------------------------------------------------------------------------
	  if(ispec.eq.2) then

	     call mc_shms(p_spec, th_spec, dpp_s, x_s, y_s, z_s, 
     >          dxdz_s, dydz_s,
     >		x_fp, dx_fp, y_fp, dy_fp, m2, shms_spec,
     >		ms_flag, wcs_flag, decay_flag, resmult, xtar_init, ok_spec, 
     >          pathlen, 5)

	     if (spec_ntuple) then
		shms_spec(59) = shmsSTOP_id
		shms_spec(58)= ytar_init
		shms_spec(53)= dpp_init
		shms_spec(54)= dph_init ! dx/dz (mr)
		shms_spec(55)= dth_init ! dy/dz (mr)
c            if (ok_spec) spec(58) =1.
		do ivar=1,SpecNtupleSize
		   write(SpecNtupleIO) shms_spec(ivar)
		enddo
	     endif
	  elseif(ispec.eq.1) then
	     call mc_hms(p_spec, th_spec, dpp_s, x_s, y_s, z_s, 
     >          dxdz_s, dydz_s,
     >          x_fp, dx_fp, y_fp, dy_fp, m2,
     >          ms_flag, wcs_flag, decay_flag, resmult, xtar_init, ok_spec, 
     >          pathlen)
	  else
	     write(6,*) 'Unknown spectrometer! Stopping..'
	     stop
	  endif

	  if (ok_spec) then !Success, increment arrays
	    dpp_recon = dpp_s
            dth_recon = dydz_s*1000.			!mr
	    dph_recon = dxdz_s*1000. !mr
C       RLT: To keep RHCS consistency with HCANA, SHMS ztar must pick up minus sign	    
	    if(ispec.eq.2) then
	       ztar_recon = - y_s / sin_ts ! Flipped +ztar->-ztar to keep RHCS
	    elseif(ispec.eq.1) then
	       ztar_recon = + y_s / sin_ts ! Keep +ztar to keep RHCS (y_s->-y_s for HMS)
	    else
	       write(6,*) 'Unknown spectrometer! Stopping..'
	       stop
	    endif
            ytar_recon = y_s

C Compute sums for calculating reconstruction variances.
	    dpp_var(1) = dpp_var(1) + (dpp_recon - dpp_init)
	    dth_var(1) = dth_var(1) + (dth_recon - dth_init)
	    dph_var(1) = dph_var(1) + (dph_recon - dph_init)
	    ztg_var(1) = ztg_var(1) + (ztar_recon - ztar_init)

	    dpp_var(2) = dpp_var(2) + (dpp_recon - dpp_init)**2
	    dth_var(2) = dth_var(2) + (dth_recon - dth_init)**2
	    dph_var(2) = dph_var(2) + (dph_recon - dph_init)**2
	    ztg_var(2) = ztg_var(2) + (ztar_recon - ztar_init)**2
	 endif			!Incremented the arrays


C Output NTUPLE entry.
C This is ugly, but want the option to have different outputs
C for spectrometer ntuples
	 if(ispec.eq.2) then
	    if (store_all.OR.(hut_ntuple.AND.ok_spec)) then
	       shms_hut(1) = x_fp
	       shms_hut(2) = y_fp
	       shms_hut(3) = dx_fp
	       shms_hut(4) = dy_fp
	       shms_hut(5) = ztar_init
	       shms_hut(6) = ytar_init
	       shms_hut(7) = dpp_init
	       shms_hut(8) = dth_init/1000.
	       shms_hut(9) = dph_init/1000.
	       shms_hut(10) = ztar_recon
	       shms_hut(11) = ytar_recon
	       shms_hut(12)= dpp_recon
	       shms_hut(13)= dth_recon/1000.
	       shms_hut(14)= dph_recon/1000.
	       shms_hut(15)= xtar_init
	       shms_hut(16)= fry
	       if (use_front_sieve) then
		  shms_hut(17)= xsfr_num
		  shms_hut(18)= ysfr_num
		  shms_hut(19)= xc_frsieve
		  shms_hut(20)= yc_frsieve
	       endif
	       if (use_sieve) then
		  shms_hut(17)= xs_num
		  shms_hut(18)= ys_num
		  shms_hut(19)= xc_sieve
		  shms_hut(20)= yc_sieve
	       endif
	       shms_hut(21)= shmsSTOP_id
	       shms_hut(22)= x
	       shms_hut(23)= y
	       shms_hut(24)= ebeam_model
	       shms_hut(25)= Eprime
	       shms_hut(26)= Q2_model
	       shms_hut(27)= W2_model
	       shms_hut(28)= W_model
	       shms_hut(29)= nu_model
	       shms_hut(30)= xbj_model
	       shms_hut(31)= theta_model
	       shms_hut(32)= F1_model
	       shms_hut(33)= F2_model
 	       shms_hut(34)= mott_nb
 	       shms_hut(35)= sigma_f1f2
 	       shms_hut(36)= sigma_weight
 	       shms_hut(37)= rate_f1f2
 	       do ivar=1,37
		  write(NtupleIO) shms_hut(ivar)
	       enddo
	    endif
	 endif

	 if(ispec.eq.1) then
	    if (store_all.OR.(hut_ntuple.AND.ok_spec)) then
	       hms_hut(1) = x_fp
	       hms_hut(2) = y_fp
	       hms_hut(3) = dx_fp
	       hms_hut(4) = dy_fp
	       hms_hut(5) = xtar_init
	       hms_hut(6) = ytar_init
	       hms_hut(7) = dph_init/1000.
	       hms_hut(8) = dth_init/1000.
	       hms_hut(9) = ztar_init
	       hms_hut(10)= dpp_init
	       hms_hut(11)= ytar_recon
	       hms_hut(12)= dph_recon/1000.
	       hms_hut(13)= dth_recon/1000.
	       hms_hut(14)= ztar_recon
	       hms_hut(15)= dpp_recon
	       hms_hut(16)= fry
	       if (use_sieve) then
		  hms_hut(17)= xs_num
		  hms_hut(18)= ys_num
		  hms_hut(19)= xc_sieve
		  hms_hut(20)= yc_sieve
	       endif
               hms_hut(21)=hSTOP_id
	       hms_hut(22)= x
	       hms_hut(23)= y
	       hms_hut(24)= ebeam_model
	       hms_hut(25)= Eprime
	       hms_hut(26)= Q2_model
	       hms_hut(27)= W2_model
	       hms_hut(28)= W_model
	       hms_hut(29)= nu_model
	       hms_hut(30)= xbj_model
	       hms_hut(31)= theta_model
	       hms_hut(32)= F1_model
	       hms_hut(33)= F2_model
 	       hms_hut(34)= mott_nb
 	       hms_hut(35)= sigma_f1f2
 	       hms_hut(36)= sigma_weight
 	       hms_hut(37)= rate_f1f2
 	       do ivar=1,37
		  write(NtupleIO) hms_hut(ivar)
	       enddo
	    endif
	 endif

C We are done with this event, whether GOOD or BAD.
C Loop for remainder of trials.

500	  continue

	enddo				!End of M.C. loop

C------------------------------------------------------------------------------C
C                           End of Monte-Carlo loop                            C
C------------------------------------------------------------------------------C

C Close NTUPLE file.

	close(NtupleIO)
	if (spec_ntuple) close(SPecNtupleIO)


	write (chanout,1002)
	write (chanout,1003) p_spec,th_spec*degrad
        write (chanout,1004) (gen_lim(i),i=1,6)

	write (chanout,1005) n_trials

	if(ispec.eq.1) then
	   armSTOP_successes=hSTOP_successes
	   armSTOP_trials=hSTOP_trials
	elseif(ispec.eq.2) then
	   armSTOP_successes=shmsSTOP_successes
	   armSTOP_trials=shmsSTOP_trials
	endif

C Indicate where particles are lost in spectrometer.
	if(ispec.eq.2) then
	   write (chanout,1015)
     >        shmsSTOP_targ_hor,shmsSTOP_targ_vert,shmsSTOP_targ_oct,
     >        shmsSTOP_FRONTSLIT_hor,shmsSTOP_FRONTSLIT_vert,
     >        shmsSTOP_HB_in,shmsSTOP_HB_men,shmsSTOP_HB_mex,
     >        shmsSTOP_HB_out,shmsSTOP_DOWNSLIT,
     >        shmsSTOP_COLL_hor,shmsSTOP_COLL_vert,shmsSTOP_COLL_oct,
     >        shmsSTOP_Q1_in,shmsSTOP_Q1_men,
     >        shmsSTOP_Q1_mid,shmsSTOP_Q1_mex,shmsSTOP_Q1_out,
     >        shmsSTOP_Q2_in,shmsSTOP_Q2_men,shmsSTOP_Q2_mid,
     >        shmsSTOP_Q2_mex,shmsSTOP_Q2_out,
     >        shmsSTOP_Q3_in,shmsSTOP_Q3_men,shmsSTOP_Q3_mid,
     >        shmsSTOP_Q3_mex,shmsSTOP_Q3_out,
     >        shmsSTOP_D1_in,shmsSTOP_D1_flr,shmsSTOP_D1_men,
     >        shmsSTOP_D1_mid1,shmsSTOP_D1_mid2,shmsSTOP_D1_mid3,
     >        shmsSTOP_D1_mid4,shmsSTOP_D1_mid5,shmsSTOP_D1_mid6,
     >        shmsSTOP_D1_mid7,shmsSTOP_D1_mex,shmsSTOP_D1_out,
     >        shmsSTOP_BP_in, shmsSTOP_BP_out

	   write (chanout,1006)
     >	   shmsSTOP_trials,shmsSTOP_hut,shmsSTOP_dc1,shmsSTOP_dc2,
     >     shmsSTOP_s1,shmsSTOP_s2,shmsSTOP_s3,shmsSTOP_cal,
     >     shmsSTOP_successes,shmsSTOP_successes

	elseif(ispec.eq.1) then
	   write (chanout,1016)
     >        hSTOP_fAper_hor,hSTOP_fAper_vert,hSTOP_fAper_oct,
     >        hSTOP_bAper_hor,hSTOP_bAper_vert,hSTOP_bAper_oct,
     >        hSTOP_slit,
     >        hSTOP_Q1_in,hSTOP_Q1_mid,hSTOP_Q1_out,
     >        hSTOP_Q2_in,hSTOP_Q2_mid,hSTOP_Q2_out,
     >        hSTOP_Q3_in,hSTOP_Q3_mid,hSTOP_Q3_out,
     >        hSTOP_D1_in,hSTOP_D1_out,
     >        hSTOP_VP1,hSTOP_VP2,hSTOP_VP3,hSTOP_VP4

	   write (chanout,1007)
     >	   hSTOP_trials,hSTOP_hut,hSTOP_dc1,hSTOP_dc2,hSTOP_scin,hSTOP_cal,
     >     hSTOP_successes,hSTOP_successes
	endif


C Compute reconstruction resolutions.

	if (armSTOP_successes.eq.0) armSTOP_successes=1
	t1 = sqrt(max(0.,dpp_var(2)/armSTOP_successes 
     > - (dpp_var(1)/armSTOP_successes)**2))
	t2 = sqrt(max(0.,dth_var(2)/armSTOP_successes 
     > - (dth_var(1)/armSTOP_successes)**2))
	t3 = sqrt(max(0.,dph_var(2)/armSTOP_successes 
     > - (dph_var(1)/armSTOP_successes)**2))
	t4 = sqrt(max(0.,ztg_var(2)/armSTOP_successes 
     > - (ztg_var(1)/armSTOP_successes)**2))

	write (chanout,1011) dpp_var(1)/armSTOP_successes,t1,
     > dth_var(1)/armSTOP_successes,
     >		t2,dph_var(1)/armSTOP_successes,t3,
     > ztg_var(1)/armSTOP_successes,t4

	write(6,*) armSTOP_trials,' Trials',armSTOP_successes
     > ,' Successes'
	write (6,1011) dpp_var(1)/armSTOP_successes,t1,
     > dth_var(1)/armSTOP_successes,
     >		t2,dph_var(1)/armSTOP_successes,t3,
     > ztg_var(1)/armSTOP_successes,t4

C ALL done!

	stop ' '

C =============================== Format Statements ============================

1001	format(a)
1002	format('!',/,'! Uniform illumination Monte-Carlo results')
1003	format('!',/'! Spectrometer setting:',/,'!',/,
     >g11.5,' =  P  spect (MeV)',/,
     >g11.5,' =  TH spect (deg)')

1004	format('!',/'! Monte-Carlo limits:',/,'!',/,
     >  g11.5,'= GEN_LIM(1) - DP/P   (half width,% )',/,
     >  g11.5,'= GEN_LIM(2) - Theta  (half width,mr)',/,
     >  g11.5,'= GEN_LIM(3) - Phi    (half width,mr)',/,
     >  g11.5,'= GEN_LIM(4) - HORIZ (full width of 3 sigma cutoff,cm)',/,
     >  g11.5,'= GEN_LIM(5) - VERT  (full width of 3 sigma cutoff,cm)',/,
     >  g11.5,'= GEN_LIM(6) - Z      (Full width,cm)')

!inp     >	,/,
!inp     >	g18.8,' =  Hor. 1/2 gap size (cm)',/,
!inp     >	g18.8,' =  Vert. 1/2 gap size (cm)')

1005	format('!',/,'! Summary:',/,'!',/,
!     >	i,' Monte-Carlo trials:')
     >  i11,' Monte-Carlo trials:')

1006	format(i11,' Initial Trials',/
     >  i11,' Trials made it to the hut',/
     >  i11,' Trial cut in dc1',/
     >  i11,' Trial cut in dc2',/
     >  i11,' Trial cut in s1',/
     >  i11,' Trial cut in s2',/
     >  i11,' Trial cut in s3',/
     >  i11,' Trial cut in cal',/
     >  i11,' Trials made it thru the detectors and were reconstructed',/
     >  i11,' Trials passed all cuts and were histogrammed.',/
     >  )

 1007	format(i11,' Initial Trials',/
     >  i11,' Trials made it to the hut',/
     >  i11,' Trial cut in dc1',/
     >  i11,' Trial cut in dc2',/
     >  i11,' Trial cut in scin',/
     >  i11,' Trial cut in cal',/
     >  i11,' Trials made it thru the detectors and were reconstructed',/
     >  i11,' Trials passed all cuts and were histogrammed.',/
     >  )


!1008	format(8i)
!1009	format(1x,i4,g,i)
!1010	format(a,i)
1011	format(
     >  'DPP ave error, resolution = ',2g18.8,' %',/,
     >  'DTH ave error, resolution = ',2g18.8,' mr',/,
     >  'DPH ave error, resolution = ',2g18.8,' mr',/,
     >  'ZTG ave error, resolution = ',2g18.8,' cm')

1012	format(1x,16i4)

1015	format(/,
     >     i11,' stopped in the TARG APERT HOR',/
     >     i11,' stopped in the TARG APERT VERT',/
     >     i11,' stopped in the TARG APERT OCTAGON',/
     >     i11,' stopped in the FIXED FRONT SLIT HOR (id=-1)',/
     >     i11,' stopped in the FIXED FRONT SLIT VERT (id=-1)',/
     >     i11,' stopped in HB ENTRANCE (id=1)',/
     >     i11,' stopped in HB MAG ENTRANCE (id=2)',/
     >     i11,' stopped in HB MAG EXIT (id=3)',/
     >     i11,' stopped in HB EXIT (id=4)',/
     >     i11,' stopped in the DOWN SIEVE SLIT (id=99)',/
     >     i11,' stopped in the COLLIMATOR HOR (id=5)',/
     >     i11,' stopped in the COLLIMATOR VERT (id=5)',/
     >     i11,' stopped in the COLLIMATOR OCTAGON (id=5)',/
     >     i11,' stopped in Q1 ENTRANCE (id=6)',/
     >     i11,' stopped in Q1 MAG ENTRANCE (id=7)',/
     >     i11,' stopped in Q1 MIDPLANE (id=8)',/
     >     i11,' stopped in Q1 MAG EXIT (id=9)',/
     >     i11,' stopped in Q1 EXIT (id=10)',/
     >     i11,' stopped in Q2 ENTRANCE (id=11)',/
     >     i11,' stopped in Q2 MAG ENTRANCE (id=12)',/
     >     i11,' stopped in Q2 MIDPLANE (id=13)',/
     >     i11,' stopped in Q2 MAG EXIT (id=14)',/
     >     i11,' stopped in Q2 EXIT (id=15)',/
     >     i11,' stopped in Q3 ENTRANCE (id=16)',/
     >     i11,' stopped in Q3 MAG ENTRANCE (id=17)',/
     >     i11,' stopped in Q3 MIDPLANE (id=18)',/
     >     i11,' stopped in Q3 MAG EXIT (id=19)',/
     >     i11,' stopped in Q3 EXIT (id=20)',/
     >     i11,' stopped in D1 ENTRANCE (id=21)',/
     >     i11,' stopped in D1 FLARE (id=22)',/
     >     i11,' stopped in D1 MAG ENTRANCE (id=23)',/
     >     i11,' stopped in D1 MID-1 (id=24)',/
     >     i11,' stopped in D1 MID-2 (id=25)',/
     >     i11,' stopped in D1 MID-3 (id=26)',/
     >     i11,' stopped in D1 MID-4 (id=27)',/
     >     i11,' stopped in D1 MID-5 (id=28)',/
     >     i11,' stopped in D1 MID-6 (id=29)',/
     >     i11,' stopped in D1 MID-7 (id=30)',/
     >     i11,' stopped in D1 MAG EXIT (id=31)',/
     >     i11,' stopped in D1 EXIT (id=32)',/
     >     i11,' stopped in BP ENTRANCE',/
     >     i11,' stopped in BP EXIT',/
     >     )

 1016 format(/,
     >     i11,' stopped in the Front-end Aperture HOR (id=5)',/
     >     i11,' stopped in the Front-end Aperture VERT (id=5)',/
     >     i11,' stopped in the Front-end Aperture OCTAGON (id=5)',/
     >     i11,' stopped in the Back-end Aperture HOR (id=6)',/
     >     i11,' stopped in the Back-end Aperture VERT (id=6)',/
     >     i11,' stopped in the Back-end Aperture OCTAGON (id=6)',/
     >     i11,' stopped in Sieve Slit (id=99)',/
     >     i11,' stopped in Q1 MAG ENTRANCE (id=7)',/
     >     i11,' stopped in Q1 MIDPLANE (id=8)',/
     >     i11,' stopped in Q1 MAG EXIT (id=9)',/
     >     i11,' stopped in Q2 MAG ENTRANCE (id=12)',/
     >     i11,' stopped in Q2 MIDPLANE (id=13)',/
     >     i11,' stopped in Q2 MAG EXIT (id=14)',/
     >     i11,' stopped in Q3 MAG ENTRANCE (id=17)',/
     >     i11,' stopped in Q3 MIDPLANE (id=18)',/
     >     i11,' stopped in Q3 MAG EXIT (id=19)',/
     >     i11,' stopped in D1 ENTRANCE (id=23)',/
     >     i11,' stopped in D1 EXIT (id=31)',/
     >     i11,' stopped in Vacuum Pipe Plane-1 (id=32)',/
     >     i11,' stopped in Vacuum Pipe Plane-2 (id=33)',/
     >     i11,' stopped in Vacuum Pipe Plane-3 (id=34)',/
     >     i11,' stopped in Vacuum Pipe Plane-4 (id=35)',/
     >     )

1100	format('!',79('-'),/,'! ',a,/,'!')
1200	format(/,'! ',a,' Coefficients',/,/,
     >  (5(g18.8,','))
     >  )
1300	format(/,'! ',a,' Coefficient uncertainties',/,/,
     >  (5(g18.8,','))
     >  )

	end
