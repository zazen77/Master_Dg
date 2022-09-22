       Program bssa14_gm_loop_tmrv

! Program to compute ground motions from
! NGA-West2 provisional GMPE

! Uses control file:
!! Control file for program bssa14_gm_loop_tmrv.for
!!Revision of program involving a change in the parameter file on this date:
!   04/20/15
!!Header to add to output file (no "!" at beginning!)
! [blank]
!!file containing regression coefficients:
!C:\nga_w2\paper4eqspectra_bssa14\BSSA14_Coefficients_071314_Revisedf4_071514.csv
!!header lines to skip:
! 2
!!name of output file:
!  bssa14_gm_loop_tmrv.out
!! mech, region, z1, irelation:   
!! mech (mech= 0, 1, 2, 3 for unspecified, ss,ns,rs), 
!! region(0=unspecified,1=CA+Taiwan+New Zealand;2=China+Turkey;3=Italy,Japan), z1, irelation
!!   Sediment depth parameters:  The equations use the difference between the depth, in km,
!!     to a shear-wave velocity of 1 km/s and the depth for the specified Vs30.   This
!!     differential depth (delta_z1 in Boore et al, 2014) is not specified directly here,
!!     but is computed in the program from z1 and a specified relation between 
!!     depth and Vs30 (mu_z1 in equation 10 in Boore et al, 2014).
!     Two possibilities are included for this relation: 1=California, 2= Japan
!     As shown by the code fragment below, delta_z1 = 0 if z1 < 0 or if z1>= 0 but 
!     the mu1_vs30 relation (irelation in the code) is neither 1 or 2):
!        if (z1 < 0.0) then
!          delta_z1 = 0.0
!        else
!          if (irelation == 1 .or. irelation == 2) then
!            delta_z1 = z1 - mu1_vs30(v30, irelation)/1000.0
!          else
!            delta_z1 = 0.0
!          end if
!        end if
!  1       0	-1.0    1
!!Order of loops over T, M, R, and Vs (the first letter in the 4-character
!! string is the outer loop, etc):
!! The following are allowed as of 04/04/11: VTMR, VTRM, VMRT, and TMRV
! VTMR
!!log-spaced periods (0) or individual periods (1)
! 1
!!If log-spaced, compute pgv and pga also? (Y/N):
!! Note: need an entry (not used) even if individual periods
!  Y
!!if log-spaced periods: nper, per_start, per_stop:
!!if individual periods: nper, list of nper periods:
!! 200 0.01 100.0
! 5 0.01, 0.1, 0.2, 1.0, 2.0
!!linearly spaced magnitudes (0) or individual magnitudes (1)?
! 0
!!if individual magnitudes: nmag, list of nmag magnitudes:
!!if linearly spaced magnitudes: m_start, delta_m, m_stop
!! 3 4.0 6.0 7.0
! 4.0 4 8.0
!!log-spaced distances (0) or individual distances (1)
! 0
!!if log-spaced distances: nr, r_start, r_stop:
!!if individual distances: nr, list of nr distances:
! 20 1.0 400.0
!! 3 10.0 40.0 80.0
!!log-spaced V30 (0) or individual V30 (1)
! 0
!!if log-spaced V30: nv30, v30_start, v30_stop:
!!if individual V30: nv30, list of nv30 distances:
! 20 180.0 1100.0
!! 2 255.0 760.0
    
! Dates: 04/20/15 - written by D.M. Boore, based on bssa14_gm_tmr and tmrs_loop_rv_drvr
!        01/25/16 - Check for desired period being greater than the maximum period of the
!                   coefficient file and reset to the maximum period.
!        10/03/17 - More precision for period output, and add frequency to the output.
!        01/31/18 - Replaced "write(*, '(a\)')"  with "write(*, '(a)')" (previous version did not
!                   compile on a Mac, according to Larry Baker).


 
! change fm_gmpe, fd_gmpe to fe_gmpe, fp_gmpe, etc.


      implicit none

      integer, PARAMETER :: 
     :                      nper_max = 120,
     :                      ncoeff_stage1 = 3,
     :                      ncoeff_stage2 = 7,
     :                      nregions = 3
       
      real per_desired    
      

      character :: loop_order*10
      integer   :: nc_loop_order

      integer :: iperflag, nper, iper_offset, iper
      real    :: per_start, per_stop, dlogPer, per(1000)
      logical :: log_spaced_periods
      
      integer :: imagflag, nmag, imag
      real    :: amag_in(500), m_start, dm, m_stop
      
      integer :: idistflag, nr, ir
      real    :: r_in(1000), r_start, dlogr, r_stop
      logical :: log_spaced_distances
      
      integer :: iv30flag, nv30, iv30 
      real    :: v30_start, v30_stop, dlogv30, v30_in(1000), v30
      logical :: log_spaced_v30
      
      character buf*300, buf4*4, f_ctl*300, f_file*300, 
     :          f_coeff*300,
     :          f_out*300,
     :          header4output*80,  
     :          format_out*80, y_type*10
     
      
      character coeff_head(10)*3
      
      character cmnts2skip(200)*79
      
      character date_ctl_correct*8, date_ctl_in*30

      logical file_exist, compute_pgv_pga 
      
      real per_stage1(120), per_stage2(120)
      
      real per_gmpe(120), e_gmpe(0:6,120), amh_gmpe(120), c_gmpe(3,120), 
     :     amref_gmpe(120), rref_gmpe(120), h_gmpe(120)
      
      real delta_c3(0:nregions,120),
     :     clin(120), vclin(120), vref(120),
     :     f1(120), f3(120), f4(120), f5(120),
     :     f6(120), f7(120) 
     
      real r1(120), r2(120), delta_phiR(120), 
     :     delta_phiV(120), v1(120), v2(120), 
     :     phi1(120), phi2(120),
     :     tau1(120), tau2(120)
      
      real e_4func(0:6)
      
      real m, logy_pred, rjb, z1
      real lny, y 
      
      integer mech, iregion, number_of_headers2skip
 
      integer :: nc_date_ctl_in, 
     :  nc_date_ctl_correct, 
     :  nc_h4out, 
     :  nc_buf, 	   
     :  nc_f_coeff,
     :  i, 
     :  j, 
     :  k,
     :  i_read_status, 
     :  indx_permax, 
     :  index_pgv, 
     :  nc_f_out, 
     :  nu_out, 
     :  nc_buf4, 
     :  indx_per, 
     :  index_pga,
     :  iflag_per,
     :  indx1,  
     :  indx2, 
     :  nu_debug, nc_f_ctl, nu_ctl, nc_cmnts2skip
     
      integer :: 
     :  nu_coeff, nper_gmpe,
     :  irelation

      real r, r_pga4nl
      
      real :: sigt, expsigt,
     :  fs, fsb,
     :  fp, fpb,
     :  fe,
     :  amp_total, 
     :  amp_nl, 
     :  amp_lin, 
     :  fp_pga4nl, fpb_pga4nl,
     :  fe_pga4nl,
     :  gspread_pga4nl,
     :  pga4nl, 
     :  expsiglny,
     :  per1, 
     :  per2,
     :  per_max, 
     :  sigt_dummy, 
     :  slope_logy, 
     :  slope_phi, 
     :  slope_tau, 
     :  slope_sigma, 
     :  phi, tau, sigma,
     :  yg, 
     :  ycgs,
     :  pi, twopi,
     :  gspread,
     :  weight

      real :: 
     :  r_t1,
     :  y_t1,
     :  fsb_t1, 
     :  fs_t1, 
     :  fpb_t1,
     :  fp_t1,
     :  fe_t1,
     :  amp_total_t1, 
     :  amp_nl_t1, 
     :  amp_lin_t1, 
     :  phi_t1, tau_t1, sigma_t1,
     :  gspread_t1
     
      real :: 
     :  r_t2,
     :  y_t2,
     :  fsb_t2, 
     :  fs_t2, 
     :  fpb_t2,
     :  fp_t2,
     :  fe_t2,
     :  amp_total_t2, 
     :  amp_nl_t2, 
     :  amp_lin_t2, 
     :  phi_t2, tau_t2, sigma_t2,
     :  gspread_t2
     
      real ::
     :  slope_fe, 
     :  slope_fpb,
     :  slope_fp,
     :  slope_fs,
     :  slope_fsb,
     :  slope_logamp_lin,
     :  slope_logamp_nl,
     :  slope_logamp_total,
     :  slope_logr,
     :  slope_gspread

      pi = 4.0*atan(1.0)
      twopi = 2.0 * pi
  
!DEBUG
      call get_lun(nu_debug)
      open(unit=nu_debug,file='debug.out',status='unknown')
!DEBUG

      c_gmpe = 0.0
       
      e_gmpe = 0.0
      
      delta_c3 = 0.0
       
      file_exist = .false.
      do while (.not. file_exist)
        f_ctl = ' '
        write(*, '(a)') 
     :  ' Enter name of control file '//
     :               '("Enter" = bssa14_gm_loop_tmrv.ctl): '
        read(*, '(a)') f_ctl
        if (f_ctl == ' ') f_ctl = 'bssa14_gm_loop_tmrv.ctl'
        call trim_c(f_ctl,nc_f_ctl)
        inquire(file=f_ctl(1:nc_f_ctl), exist=file_exist)
        if (.not. file_exist) then
          write(*,'(a)') ' ******* FILE DOES NOT EXIST ******* '
        end if
      end do
      
      call get_lun(nu_ctl)
      open(unit=nu_ctl,file=f_ctl(1:nc_f_ctl),status='unknown')

!Check date of control file:
      call skipcmnt(nu_ctl,cmnts2skip, nc_cmnts2skip)
      date_ctl_in = ' '
      read(nu_ctl,'(a)') date_ctl_in
      call trim_c(date_ctl_in,nc_date_ctl_in)
      
      date_ctl_correct = ' '
      date_ctl_correct = '04/20/15'
      call trim_c(date_ctl_correct,nc_date_ctl_correct)


      if (date_ctl_correct(1:nc_date_ctl_correct) /= 
     :    date_ctl_in(1:nc_date_ctl_in)) then
        write(*,'(a)') 
     :     ' The control file has the wrong date; update your '//
     :       'control file and rerun the program!'
        write(*,'(a)') 
     :     ' The date for the control file you used is '
     :       //date_ctl_in(1:nc_date_ctl_in)
        write(*,'(a)') 
     :     ' The correct date is '
     :       //date_ctl_correct(1:nc_date_ctl_correct)
        close(nu_ctl)
        stop
      end if
!Check date of control file

      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)  
      header4output = ' '
      read(nu_ctl,'(a)') header4output
      call trim_c(header4output,nc_h4out)
      
      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)      

      buf = ' '
      read(nu_ctl,'(a)') buf
      call trim_c(buf,nc_buf)
      f_coeff  = ' '
      f_coeff  = buf(1:nc_buf)
      call trim_c(f_coeff, nc_f_coeff)
 
      file_exist = .false.
      inquire(file=f_coeff(1:nc_f_coeff), exist=file_exist)
      if (.not. file_exist) then
        write(*,'(a)') ' ******* Coefficient file '//
     :         f_coeff(1:nc_f_coeff)//
     :        ' does not exist: QUITTING ******* '
        close(nu_ctl)
        stop
      end if

      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)      
      read(nu_ctl,*) number_of_headers2skip
      
! Open coefficient files and
! read coefficients for regression coefficients:

      print *,' f_coeff  = '//
     :          f_coeff(1:nc_f_coeff)
      
      call get_lun(nu_coeff)
      open(unit=nu_coeff,file=f_coeff(1:nc_f_coeff),
     :     status='unknown')

      call skip(nu_coeff,number_of_headers2skip)  ! The column header row is NOT preceded with '!'
                             ! because might want to open the coefficient file
                             ! in Excel.

 
      i = 0
      
      read_coefficients: DO

!Coefficient file:
! 1	         2	3	4	5	6	7	8	9	10	11	12	13	14	          15	        16	         17	                 18	         19	     20	      21	22	23	24	25	26	        27	   28	29	30	31	32	33	34	35	36	37
!Period(sec)	e0	e1	e2	e3	e4	e5	e6	Mh	c1	c2	c3	Mref	Rref(km)	h(km)	Dc3(globalCATWNZ)	Dc3(ChinaTurkey)	Dc3(ItalyJapan)	  c	Vc(m/s)	Vref(m/s)	f1	f3	f4	f5	f6(1/km)	f7	R1(km)	R2(km)	DfR	DfV	V1(m/s)	V2(m/s)	f1	f2	t1	t2
!period	        e0	e1	e2	e3	e4	e5	e6	Mh	c1	c2	c3	Mref	Rref	        h	Dc3.globalCATWNZ	Dc3.ChinaTurkey	        Dc3.ItalyJapan	clin	Vc	Vref	        f1	f3	f4	f5	f6	        f7	R1	R2	DfR	DfV	V1	V2	f1	f2	t1	t2

!     :                      ncoeff_stage1 = 3,
!     :                      ncoeff_stage2 = 7,
!     :                      nregions = 3

        read(nu_coeff,*,iostat=i_read_status) 
     :      per_gmpe(i+1), 
     :      (e_gmpe(j,i+1), j = 0, ncoeff_stage2-1),
     :      amh_gmpe(i+1),
     :      (c_gmpe(j,i+1), j = 1, ncoeff_stage1), 
     :      amref_gmpe(i+1), rref_gmpe(i+1), h_gmpe(i+1), 
     :      (delta_c3(j,i+1), j = 1, nregions),     
     :      clin(i+1), vclin(i+1), vref(i+1),
     :      f1(i+1), f3(i+1), f4(i+1), f5(i+1),f6(i+1), f7(i+1),
     :      r1(i+1), r2(i+1), delta_phiR(i+1), delta_phiV(i+1),
     :      v1(i+1), v2(i+1),
     :      phi1(i+1), phi2(i+1), tau1(i+1), tau2(i+1)
     
     
        if (i_read_status /= 0) exit
        
        i = i + 1
      
      END DO read_coefficients
      
      close(nu_coeff)
      
      nper_gmpe = i
      
      
 
!DEBUG
      write(nu_debug,*)
      DO i = 1, nper_gmpe
        WRITE(nu_debug,*) 
     :   per_gmpe(i), 
     :   (e_gmpe(j,i), j = 0, ncoeff_stage2-1),
     :   amh_gmpe(i), (c_gmpe(j,i), j = 1, ncoeff_stage1), 
     :   amref_gmpe(i), rref_gmpe(i), h_gmpe(i)
      END DO
      write(nu_debug,*)
!DEBUG
        
! Find the index of the maximum period
      per_max = -1.0
      DO i = 1, nper_gmpe
        IF (per_gmpe(i) > per_max) THEN
          per_max = per_gmpe(i)
          indx_permax = i
        END IF
      END DO
     
!DEBUG
      write(nu_debug,*)
      write(nu_debug,'(a, 1x,i3,1x,f7.3)') 
     :   ' indx_permax, per_max ', indx_permax, per_max
!DEBUG

!    Find index for pgv (period = -1.0):
      index_pgv = 0
      find_pgv_index: DO i = 1, nper_gmpe
        IF (per_gmpe(i) == -1.0) THEN
          index_pgv = i
          EXIT find_pgv_index
        END IF
      END DO find_pgv_index
      
!    Find index for pga (period = 0.0):
      index_pga = 0
      find_pga_index: DO i = 1, nper_gmpe
        IF (per_gmpe(i) == 0.0) THEN
          index_pga = i
          EXIT find_pga_index
        END IF
      END DO find_pga_index
      
      IF (index_pga == 0) THEN
        PRINT *,' ERROR; in ngaw2_gm_tmr, '//
     :        'coefficients for pga not found; QUIT'
          STOP
      END IF

      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)      
      f_out = ' '
      read(nu_ctl,'(a)') f_out
      call trim_c(f_out,nc_f_out)
      call get_lun(nu_out)
      open(unit=nu_out,file=f_out(1:nc_f_out),status='unknown')


      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)      
      read(nu_ctl, * ) mech, iregion, z1, irelation
 
      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)      
      loop_order = ' '
      read(nu_ctl,'(a)') loop_order
      call upstr(loop_order)
      call trim_c(loop_order, nc_loop_order)

      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)      
      read(nu_ctl,*) iperflag

      buf = ' '
      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)      
      read(nu_ctl,'(a)') buf
      call trim_c(buf, nc_buf)
      call upstr(buf)
      if (buf(1:1) == 'Y') then
        compute_pgv_pga = .true.
      else
        compute_pgv_pga = .false.
      end if

      
      if (iperflag == 0) then
        log_spaced_periods = .true.
        call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)      
        read(nu_ctl,*) nper, per_start, per_stop
        if (per_start <= 0.0) then
          write(*,*) ' per_start <= 0.0, which is not allowed'//
     :      ' for log-spaced periods. QUITTING'
          close(nu_ctl)
          STOP
        end if        
        if (nper == 1) then
          dlogPer = 0.0
        else
          dlogPer = alog10(per_stop/per_start)/real(nper-1)
        end if
        if (compute_pgv_pga) then
          nper = nper+2
          per(1) = -1.0
          per(2) = 0.0
          iper_offset = 2
        else
          iper_offset = 0
        end if
        do iper=1,nper-iper_offset
          per(iper+iper_offset) = per_start*10**(real(iper-1)*dlogPer)
        end do
      else
        log_spaced_periods = .false.
        call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)      
        read(nu_ctl,*) nper, (per(iper), iper = 1, nper)
      end if        
      
      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)      
      read(nu_ctl,*) imagflag
      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)      
      if (imagflag ==1) then  ! individual magnitudes
        read(nu_ctl,*) nmag, (amag_in(imag), imag=1, nmag)
      else
        read(nu_ctl,*) m_start, dm, m_stop
        nmag = nint((m_stop - m_start)/dm) + 1
        do imag=1, nmag
          amag_in(imag) = m_start + real(imag-1)*dm
        end do
      end if        

      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)      
      read(nu_ctl,*) idistflag
      if (idistflag == 0) then
        log_spaced_distances = .true.
        call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)      
        read(nu_ctl,*) nr, r_start, r_stop
        if (r_start <= 0.0) then
          write(*,*) ' r_start = 0.0, which is not allowed; '//
     :               ' use a small nonzero number instead; QUITTING'
          close(nu_ctl)
          STOP
        end if        
        if (nr == 1) then
          dlogr = 0.0
        else
          dlogr = alog10(r_stop/r_start)/real(nr-1)
        end if
        do ir=1,nr
          r_in(ir) = r_start*10**(real(ir-1)*dlogr)
        end do
      else
        log_spaced_distances = .false.
        call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)      
        read(nu_ctl,*) nr, (r_in(ir), ir = 1, nr)
      end if        
      
      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)      
      read(nu_ctl,*) iv30flag
      if (iv30flag == 0) then
        log_spaced_v30 = .true.
        call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)      
        read(nu_ctl,*) nv30, v30_start, v30_stop
        if (v30_start <= 0.0) then
          write(*,*) ' v30_start = 0.0, which is not allowed; '//
     :               ' use a small nonzero number instead; QUITTING'
          close(nu_ctl)
          STOP
        end if        
        if (nv30 == 1) then
          dlogv30 = 0.0
        else
          dlogv30 = alog10(v30_stop/v30_start)/real(nv30-1)
        end if
        do iv30=1,nv30
          v30_in(iv30) = v30_start*10**(real(iv30-1)*dlogv30)
        end do
      else
        log_spaced_v30 = .false.
        call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)      
        read(nu_ctl,*) nv30, (v30_in(iv30), iv30 = 1, nv30)
      end if        

      write(nu_out,'(a)') 
     : '! f_coeff  = '//f_coeff(1:nc_f_coeff)

      write(nu_out,'(
     :                12x,a, 12x,a, 4x,a, 4x,a, 6x,a, 4x,a, 2x,a,
     :                1x,a, 4x,a,
     :                7x,a, 1x,a, 1x,a, 
     :                5x,a, 1x,a, 1x,a, 
     :                2x,a,
     :                9x,a, 
     :                8x,a, 9x,a,
     :                8x,a, 9x,a,
     :                4x,a,
     :                5x,a,
     :                4x,a, 
     :                3x,a, 4x,a,
     :                3x,a, 4x,a,
     :                2x,a, 1x,a,
     :                4x,a, 5x,a, 2x,a,
     :                5x,a, 5x,a, 3x,a
     :                                  )')
     :           'T', 'F', 'M', 'Rjb', 'R', 'V30', 'mech',
     :           'iregion', 'z1',
     :           'Y(g)', 'Y(g)/e^sig_LnY', 'Y(g)*e^sig_LnY',
     :           'Y(cgs)', 'Y(cgs)/e^sig_LnY', 'Y(cgs)*e^sig_LnY',
     :           'e^sig_lnY',
     :           'fe', 
     :           'fpb', 'fp',
     :           'fsb', 'fs',
     :           'gspread',
     :           'ln(yg)',
     :           'exp(fe)',
     :           'exp(fpb)', 'exp(fp)',
     :           'exp(fsb)', 'exp(fs)',
     :           'pga4nl(g)', 'pga4nl(cgs)',
     :           'amp_lin', 'amp_nl', 'amp_total',
     :           'phi', 'tau', 'sigma'


      
! Loop over values of TMRV:

! Select case:

      SELECT CASE (loop_order(1:nc_loop_order))
        
        CASE ('TMRV')

          DO iper = 1, nper            
            per_desired = per(iper)
            if (per_desired > per_max) per_desired = per_max
            DO imag = 1, nmag      
              m = amag_in(imag)
              DO ir = 1, nr   
                rjb = r_in(ir)
                DO iv30 = 1, nv30 
                  v30 = v30_in(iv30)
                  CALL bssa14_gm_loop_tmrv_core(
     :              nu_out, 
     :              per_desired, indx_permax, index_pga, index_pgv, 
     :              m, rjb, r, v30, mech, iregion, z1, irelation,      
     :              per_gmpe, nper_gmpe, c_gmpe, amref_gmpe, 
     :              rref_gmpe, delta_c3,
     :              e_gmpe, amh_gmpe, h_gmpe,
     :              clin, vclin, vref,
     :              f1, f3, f4, f5,
     :              f6, f7,
     :              phi1, phi2, 
     :              tau1, tau2,
     :              R1, R2, delta_phiR,
     :              v1, v2, delta_phiV    
     :                                                          )
                END DO  
              END DO  
            END DO              
          END DO  


        CASE ('VTMR')
        
          DO iv30 = 1, nv30      
            v30 = v30_in(iv30)
            DO iper = 1, nper            
              per_desired = per(iper)
              if (per_desired > per_max) per_desired = per_max
              DO imag = 1, nmag 
                m = amag_in(imag)
                DO ir = 1, nr      
                  rjb = r_in(ir)
                  CALL bssa14_gm_loop_tmrv_core(
     :              nu_out, 
     :              per_desired, indx_permax, index_pga, index_pgv, 
     :              m, rjb, r, v30, mech, iregion, z1, irelation,      
     :              per_gmpe, nper_gmpe, c_gmpe, amref_gmpe, 
     :              rref_gmpe, delta_c3,
     :              e_gmpe, amh_gmpe, h_gmpe,
     :              clin, vclin, vref,
     :              f1, f3, f4, f5,
     :              f6, f7,
     :              phi1, phi2, 
     :              tau1, tau2,
     :              R1, R2, delta_phiR,
     :              v1, v2, delta_phiV    
     :                                                          )
                END DO  
              END DO              
            END DO  
          END DO  

        CASE ('VTRM')

          DO iv30 = 1, nv30      
            v30 = v30_in(iv30)
            DO iper = 1, nper            
              per_desired = per(iper)
              if (per_desired > per_max) per_desired = per_max
              DO ir = 1, nr      
                rjb = r_in(ir)
                DO imag = 1, nmag 
                  m = amag_in(imag)
                  CALL bssa14_gm_loop_tmrv_core(
     :              nu_out, 
     :              per_desired, indx_permax, index_pga, index_pgv, 
     :              m, rjb, r, v30, mech, iregion, z1, irelation,      
     :              per_gmpe, nper_gmpe, c_gmpe, amref_gmpe, 
     :              rref_gmpe, delta_c3,
     :              e_gmpe, amh_gmpe, h_gmpe,
     :              clin, vclin, vref,
     :              f1, f3, f4, f5,
     :              f6, f7,
     :              phi1, phi2, 
     :              tau1, tau2,
     :              R1, R2, delta_phiR,
     :              v1, v2, delta_phiV    
     :                                                          )
                END DO  
              END DO              
            END DO  
          END DO  

        CASE ('VMRT')

          DO iv30 = 1, nv30      
            v30 = v30_in(iv30)
            DO imag = 1, nmag 
              m = amag_in(imag)
              DO ir = 1, nr      
                rjb = r_in(ir)
                DO iper = 1, nper            
                  per_desired = per(iper)
                  if (per_desired > per_max) per_desired = per_max
                  CALL bssa14_gm_loop_tmrv_core(
     :              nu_out, 
     :              per_desired, indx_permax, index_pga, index_pgv, 
     :              m, rjb, r, v30, mech, iregion, z1, irelation,      
     :              per_gmpe, nper_gmpe, c_gmpe, amref_gmpe, 
     :              rref_gmpe, delta_c3,
     :              e_gmpe, amh_gmpe, h_gmpe,
     :              clin, vclin, vref,
     :              f1, f3, f4, f5,
     :              f6, f7,
     :              phi1, phi2, 
     :              tau1, tau2,
     :              R1, R2, delta_phiR,
     :              v1, v2, delta_phiV    
     :                                                          )
                END DO  
              END DO              
            END DO  
          END DO  

        CASE DEFAULT
        
          PRINT *, ' Loop order string '//loop_order(1:nc_loop_order)//
     :             ' not valid; QUITTING'
          STOP

      END SELECT



!######################################

      close(nu_ctl)
      close(nu_out)
      close(nu_debug)
      stop
       
      end

! ----------------- BEGIN BSSA14_GM_LOOP_TMRV_CORE ---------------------------
      SUBROUTINE bssa14_gm_loop_tmrv_core(
     :        nu_out, 
     :        per_desired, indx_permax, index_pga, index_pgv, 
     :        m, rjb, r, v30, mech, iregion, z1, irelation,      
     :        per_gmpe, nper_gmpe, c_gmpe, amref_gmpe, 
     :        rref_gmpe, delta_c3,
     :        e_gmpe, amh_gmpe, h_gmpe,
     :        clin, vclin, vref,
     :        f1, f3, f4, f5,
     :        f6, f7,
     :        phi1, phi2, 
     :        tau1, tau2,
     :        R1, R2, delta_phiR,
     :        v1, v2, delta_phiV    
     :                                                          )

! 04/20/15 - Written, based on TMRS_RV_DRVR_CORE in TMRS_LOOP_RV_DRVR
! 10/03/17 - Add frequency to output
      
      integer, PARAMETER :: 
     :                      nper_max = 120,
     :                      ncoeff_stage1 = 3,
     :                      ncoeff_stage2 = 7,
     :                      nregions = 3
       
      character :: y_type*10
      
      integer :: indx_permax, nper_gmpe
      
      real freq_out
      
      real m

      real per_gmpe(120), e_gmpe(0:6,120), amh_gmpe(120), c_gmpe(3,120), 
     :     amref_gmpe(120), rref_gmpe(120), h_gmpe(120)
      
      real delta_c3(0:nregions,120),
     :     clin(120), vclin(120), vref(120),
     :     f1(120), f3(120), f4(120), f5(120),
     :     f6(120), f7(120) 
     
      real r1(120), r2(120), delta_phiR(120), 
     :     delta_phiV(120), v1(120), v2(120), 
     :     phi1(120), phi2(120),
     :     tau1(120), tau2(120)

        r_pga4nl = sqrt(rjb**2+h_gmpe(index_pga)**2)
        call y_bssa14_no_site_amps(
     :           m, r_pga4nl, mech, iregion,
     :           c_gmpe(:,index_pga), amref_gmpe(index_pga), 
     :           rref_gmpe(index_pga), delta_c3(:,index_pga),
     :           e_gmpe(:,index_pga), amh_gmpe(index_pga),
     :           pga4nl, fe_pga4nl, fpb_pga4nl, fp_pga4nl, 
     :           gspread_pga4nl)     
 
! Select case:

        y_type = ' '
        IF (per_desired <= 0) THEN
          y_type = 'PGVorPGA'
        ELSE
          y_type = 'PSA'        
        END IF
        
        SELECT CASE (y_type)
        
        
        CASE ('PGVorPGA')
        
            if (per_desired == -1.0) then
              IF (index_pgv == 0) THEN
                print *,' ERROR; coefficients for pgv not found; QUIT'
                stop
              END IF
              indx_per = index_pgv
            else if (per_desired == 0.0) then
              IF (index_pga == 0) THEN
                print *,' ERROR; coefficients for pga not found; QUIT'
                stop
              END IF
              indx_per = index_pga
            else
              print *,' PER = ', per_desired, 
     :            ', not -1.0 (PGV) or 0.0 (PGA); QUITTING!!'
              STOP
            end if

            r = sqrt(rjb**2+h_gmpe(indx_per)**2)
      
            call bssa14_gm_sub4y(
     :        per_desired,
     :        m, rjb, r, v30, mech, iregion, z1, irelation,      
     :        pga4nl,
     :        c_gmpe(:,indx_per), amref_gmpe(indx_per), 
     :        rref_gmpe(indx_per), delta_c3(:,indx_per),
     :        e_gmpe(:,indx_per), amh_gmpe(indx_per),
     :        clin(indx_per), vclin(indx_per), vref(indx_per),
     :        f1(indx_per), f3(indx_per), f4(indx_per), f5(indx_per),
     :        f6(indx_per), f7(indx_per),
     :        phi1(indx_per), phi2(indx_per), 
     :        tau1(indx_per), tau2(indx_per),
     :        R1(indx_per), R2(indx_per), delta_phiR(indx_per),
     :        v1(indx_per), v2(indx_per), delta_phiV(indx_per),
     :        y, fe, fpb, fp, fsb, fs, 
     :        phi, tau, sigma, 
     :        gspread,     
     :        amp_lin, amp_nl, amp_total)     
     
        CASE ('PSA')         

! Trap for per outside the range of the tabulated periods (but only if not pga or pgv):

            if (per_desired < per_gmpe(1) .or. 
     :               per_desired > per_gmpe(indx_permax))   then
              print *,' ERROR: per = ', per_desired,
     :     ' is outside the tabulated range in per_gmpe table'
              print *,' per_gmpe '
              do i = 1, nper_gmpe
                print *, per_gmpe(i)
              end do
              print *, ' QUITTING!!!'
              stop
            end if
       

! Find lower index of the two tabulated periods that span the interval
! containing the input period
! assume that the periods are sorted from small to large, and that period value
! for index nper_gmpe corresponds to the maximum period (and thus 
! the indices 1 and 2 are for pgv and pgv (per_gmpe = -1 and 0)).

            iflag_per = 0
            IF (per_desired == per_gmpe(nper_gmpe)) THEN
              indx1 = nper_gmpe-1
              indx2 = nper_gmpe
              iflag_per = 1
            ELSE
              DO i = 1, nper_gmpe-1  
                IF (per_desired >= per_gmpe(i) .and. 
     :              per_desired < per_gmpe(i+1) ) 
     :                                                       THEN
                  indx1 = i
                  indx2 = i+1
                  iflag_per = 1
                  EXIT
                END IF
              END DO
            END IF
            if (iflag_per == 0) then
              PRINT *,' ERROR: could not find per = ',per_desired,
     :      ' in per_gmpe table'
              PRINT *,' per_gmpe '
              DO i = 1, nper_gmpe
                PRINT *, per_gmpe(i)
              END DO
              PRINT *, ' QUITTING!!!'
              STOP
            end if
  
            per1 = per_gmpe(indx1)
            per2 = per_gmpe(indx2)

! Now evaluate PSA and sigma for the periods on either side of the desired period 
! (the evaluations are so quick that there is no need to treat
! the case of the desired and tabulated periods being the same).

 
! In previous version evaluate sigma first because ba08_gm_sub4y uses it as an input parameter
! (although it is only used in computing expsiglny). I no longer use expsiglny that is returned
! from ba08_gm_sub4y; for the sake of minimal complications, however,
! I am not going to remove sigt from the input arguments of ba08_gm_sub4y now, and I need to call
! ba08_gm_sub4y with a dummy value.  

            r_t1 = sqrt(rjb**2+h_gmpe(indx1)**2)
      
            call bssa14_gm_sub4y(
     :        per_desired,
     :        m, rjb, r_t1, v30, mech, iregion, z1, irelation,      
     :        pga4nl,
     :        c_gmpe(:,indx1), amref_gmpe(indx1), 
     :        rref_gmpe(indx1), delta_c3(:,indx1),
     :        e_gmpe(:,indx1), amh_gmpe(indx1),
     :        clin(indx1), vclin(indx1), vref(indx1),
     :        f1(indx1), f3(indx1), f4(indx1), f5(indx1),
     :        f6(indx1), f7(indx1),
     :        phi1(indx1), phi2(indx1), 
     :        tau1(indx1), tau2(indx1),
     :        R1(indx1), R2(indx1), delta_phiR(indx1),
     :        v1(indx1), v2(indx1), delta_phiV(indx1),
     :        y_t1, fe_t1, fpb_t1, fp_t1, fsb_t1, fs_t1, 
     :        phi_t1, tau_t1, sigma_t1, 
     :        gspread_t1, 
     :        amp_lin_t1, amp_nl_t1, amp_total_t1)    
     
     
            r_t2 = sqrt(rjb**2+h_gmpe(indx2)**2)
      
            call bssa14_gm_sub4y(
     :        per_desired,
     :        m, rjb, r_t2, v30, mech, iregion, z1, irelation,      
     :        pga4nl,
     :        c_gmpe(:,indx2), amref_gmpe(indx2), 
     :        rref_gmpe(indx2), delta_c3(:,indx2),
     :        e_gmpe(:,indx2), amh_gmpe(indx2),
     :        clin(indx2), vclin(indx2), vref(indx2),
     :        f1(indx2), f3(indx2), f4(indx2), f5(indx2),
     :        f6(indx2), f7(indx2),
     :        phi1(indx2), phi2(indx2), 
     :        tau1(indx2), tau2(indx2),
     :        R1(indx2), R2(indx2), delta_phiR(indx2),
     :        v1(indx2), v2(indx2), delta_phiV(indx2),
     :        y_t2, fe_t2, fpb_t2, fp_t2, fsb_t2, fs_t2, 
     :        phi_t2, tau_t2, sigma_t2, 
     :        gspread_t2, 
     :        amp_lin_t2, amp_nl_t2, amp_total_t2)     

     
            slope_logy =  (alog10(y_t2) - alog10(y_t1))/
     :           alog10(per2/per1)
            y = 10**(alog10(y_t1) + 
     :                slope_logy*alog10(per_desired/per1)) 
        
 
            slope_phi =  (phi_t2 - phi_t1)/
     :           alog10(per2/per1)
            phi = phi_t1 + 
     :                slope_phi*alog10(per_desired/per1) 
        
            slope_tau =  (tau_t2 - tau_t1)/
     :           alog10(per2/per1)
            tau = tau_t1 + 
     :                slope_tau*alog10(per_desired/per1) 
        
            slope_sigma =  (sigma_t2 - sigma_t1)/
     :           alog10(per2/per1)
            sigma = sigma_t1 + 
     :                slope_sigma*alog10(per_desired/per1) 
        
            slope_fe =  (fe_t2 - fe_t1)/
     :           alog10(per2/per1)
            fe =  fe_t1 + 
     :                slope_fe*alog10(per_desired/per1) 
     
            slope_fpb =  (fpb_t2 - fpb_t1)/
     :           alog10(per2/per1)
            fpb =  fpb_t1 + 
     :                slope_fpb*alog10(per_desired/per1) 
     
            slope_fp =  (fp_t2 - fp_t1)/
     :           alog10(per2/per1)
            fp =  fp_t1 + 
     :                slope_fp*alog10(per_desired/per1) 
     
            slope_fsb =  (fsb_t2 - fsb_t1)/
     :           alog10(per2/per1)
            fsb =  fsb_t1 + 
     :                slope_fsb*alog10(per_desired/per1) 
     
            slope_fs =  (fs_t2 - fs_t1)/
     :           alog10(per2/per1)
            fs =  fs_t1 + 
     :                slope_fs*alog10(per_desired/per1) 
     
            slope_gspread =  (gspread_t2 - gspread_t1)/
     :           alog10(per2/per1)
            gspread =  gspread_t1 + 
     :                slope_gspread*alog10(per_desired/per1) 
        
            slope_logamp_lin = (alog10(amp_lin_t2)-alog10(amp_lin_t1))/
     :           alog10(per2/per1)
            amp_lin = 10**(alog10(amp_lin_t1) + 
     :                slope_logamp_lin*alog10(per_desired/per1)) 
        
            slope_logamp_nl =  (alog10(amp_nl_t2) - alog10(amp_nl_t1))/
     :           alog10(per2/per1)
            amp_nl = 10**(alog10(amp_nl_t1) + 
     :                slope_logamp_nl*alog10(per_desired/per1)) 
        
            slope_logamp_total =  
     :                  (alog10(amp_total_t2) - alog10(amp_total_t1))/
     :           alog10(per2/per1)
            amp_total = 10**(alog10(amp_total_t1) + 
     :                slope_logamp_total*alog10(per_desired/per1)) 
        
            slope_logr =  
     :                  (alog10(r_t2) - alog10(r_t1))/
     :           alog10(per2/per1)
            r = 10**(alog10(r_t1) + 
     :                slope_logr*alog10(per_desired/per1)) 
     
        CASE DEFAULT
        
          PRINT *, ' Neither PGVorPGA or PSA were selected; '//
     :              'something is wrong; QUITTING'
          STOP

        END SELECT


! Convert from g to cgs units:

        if (per_desired >= 0.0) then
          yg   =y
          ycgs = 981.0 * yg
        else
          yg   = y 
          ycgs = yg
        end if

        expsiglny = exp(sigma)

        if (per_desired < 0.0) then
          freq_out = -1.0
        else if (per_desired == 0.0) then
          freq_out = 999.0
        else
          freq_out = 1.0/per_desired
        end if
        
        write(nu_out,'(
     :      1x,es12.5,1x,es12.5,1x,f4.2,1x,f6.2,1x,f6.2,1x,f6.1,4x,i2,
     :      7x,i1, 1x,f5.2,
     :      1x,es10.3, 2(5x,es10.3),
     :      1x,es10.3, 2(7x,es10.3),
     :      1x,es10.3,
     :      1x,es10.3,
     :      2(1x,es10.3),
     :      2(1x,es10.3),
     :      1x,es10.3,
     :      1x,es10.3,
     :      1x,es10.3,
     :      2(1x,es10.3),
     :      2(1x,es10.3),
     :      1x,es10.3, 2x,es10.3,
     :      3(1x,es10.3),
     :      1x,f7.4, 1x,f7.4, 1x,f7.4
     :                                                 )')
     :      per_desired, freq_out, m, rjb, r, v30, mech, 
     :      iregion, z1, 
     :      yg, yg/expsiglny, yg*expsiglny,
     :      ycgs, ycgs/expsiglny, ycgs*expsiglny,
     :      expsiglny,
     :      fe, 
     :      fpb, fp, 
     :      fsb, fs,
     :      gspread,
     :      alog(yg),
     :      exp(fe), 
     :      exp(fpb), exp(fp), 
     :      exp(fsb), exp(fs),
     :      pga4nl, 981.0*pga4nl,
     :      amp_lin, amp_nl, amp_total,
     :      phi, tau, sigma

 
 
      END
! ----------------- END BSSA14_GM_LOOP_TMRV_CORE ---------------------------


      include 'bssa14_gm_tmr_subroutines.for'  !This file is made by calling 
                                               !make_bssa14gm_tmr_subroutine_file.bat
                                               !from within cl_bssa14_gm_loop_tmrv.bat

      
!      include 'y_bssa14_no_site_amps.for'

!      include 'bssa14_gm_sub4y.for'
      
!      include '\forprogs\skipcmnt.for'
        
!      include '\forprogs\get_lun.for'
!      include '\forprogs\lin_interp.for'
!      include '\forprogs\locate.for'
!      include '\forprogs\skip.for'
!      include '\forprogs\upstr.for'
!      include '\forprogs\trim_c.for'

 
