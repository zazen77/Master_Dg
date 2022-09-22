! Control file for program bssa14_gm_loop_tmrv.for
!Revision of program involving a change in the parameter file on this date:
   04/20/15
!Header to add to output file (no "!" at beginning!)
 [blank]
!file containing regression coefficients:
C:\nga_w2\paper4eqspectra_bssa14\BSSA14_Coefficients_071314_Revisedf4_071514.csv
!header lines to skip:
 1
!name of output file:
  bssa14_gm_loop_tmrv.out
! mech, region, z1, irelation:   
! mech (mech= 0, 1, 2, 3 for unspecified, ss,ns,rs), 
! region(0=unspecified,1=CA+Taiwan+New Zealand;2=China+Turkey;3=Italy,Japan), z1, irelation
!   Sediment depth parameters:  The equations use the difference between the depth, in km,
!     to a shear-wave velocity of 1 km/s and the depth for the specified Vs30.   This
!     differential depth (delta_z1 in Boore et al, 2014) is not specified directly here,
!     but is computed in the program from z1 and a specified relation between 
!     depth and Vs30 (mu_z1 in equation 10 in Boore et al, 2014).
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
  1       0	-1.0    1
!Order of loops over T, M, R, and Vs (the first letter in the 4-character
! string is the outer loop, etc):
! The following are allowed as of 04/04/11: VTMR, VTRM, VMRT, and TMRV
 TMRV
!log-spaced periods (0) or individual periods (1)
 1
!If log-spaced, compute pgv and pga also? (Y/N):
! Note: need an entry (not used) even if individual periods
  Y
!if log-spaced periods: nper, per_start, per_stop:
!if individual periods: nper, list of nper periods:
! 200 0.01 10.0
 5 0.01, 0.1, 0.2, 1.0, 2.0
!linearly spaced magnitudes (0) or individual magnitudes (1)?
 0
!if individual magnitudes: nmag, list of nmag magnitudes:
!if linearly spaced magnitudes: m_start, delta_m, m_stop
! 3 4.0 6.0 7.0
 4.0 4 8.0
!log-spaced distances (0) or individual distances (1)
 0
!if log-spaced distances: nr, r_start, r_stop:
!if individual distances: nr, list of nr distances:
 10 1.0 200.0
! 3 10.0 40.0 80.0
!log-spaced V30 (0) or individual V30 (1)
 0
!if log-spaced V30: nv30, v30_start, v30_stop:
!if individual V30: nv30, list of nv30 distances:
 20 180.0 1100.0
! 2 255.0 760.0

 
 


 