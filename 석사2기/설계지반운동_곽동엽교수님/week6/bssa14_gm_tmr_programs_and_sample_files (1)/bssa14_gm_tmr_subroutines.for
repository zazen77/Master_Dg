
! --------------------------------------- y_bssa14_no_site_amps --------------------------------
      subroutine y_bssa14_no_site_amps(
     :           m, r, mech, iregion,
     :           c, mref, rref, delta_c3,
     :           e, mh,
     :           y, fe, fpb, fp, gspread)     
     
! **** USE THIS VERSION WITH COEFFICIENTS IN PEER REPORT FORMAT ****

! NOTE: y in g, unless pgv!!!!  

! ncoeff_s1, ncoeff_s2 are not included as arguments because it is
! assumed that they are 3 and 7, respectively.

! Assume no site amp

! mech = 0, 1, 2, 3 for unspecified, SS, NS, RS, respectively

! Dates: 02/27/13 - Modified from y_ngaw2_no_site_amps
!        07/01/13 - Renamed from y_bssa13_no_site_amps
 
      IMPLICIT none
      
      real :: m, r, 
     :        c(1:3), mref, rref, delta_c3(0:3), 
     :        e(0:6), mh, 
     :        alny, y, fpb, fp, fe, fs, fault_type_term, gspread     
     
      integer :: mech, iregion

      if ( e(mech) == 0.0) then   ! no coefficient
        print *,' From y_bssa14_no_site_amps, '//
     :       ' e(mech) == 0.0; QUITTING'
        stop
      end if
        
      if (mech /=0 .and. mech /= 1 .and. mech /= 2 
     :             .and. mech /= 3)                      then
        print *,' From y_bssa14_no_site_amps, mech = ', 
     :             mech,
     :       ' not a valid value; QUITTING'
        stop
      else
        fault_type_term = e(mech)
      end if
        

      if (m < mh ) then      
        fe = e(4)*(m-mh)
     :       + e(5)*(m-mh)**2 
      else
        fe = e(6)*(m-mh)
      end if
 
      fe = fault_type_term + fe
 
      gspread = c(1) + c(2)*(m-mref)

      fpb = gspread*alog(r/rref)
     :   + c(3)*(r-rref)

! Note: delta_c3 is initialized to 0, and only values for indices 1,2,3
! are read in.  So specifying iregion = 0 should be OK
      fp = fpb + delta_c3(iregion)*(r-rref)
      
      fs = 0.0  ! No amps
         
      alny = fe + fp + fs 
      
      y = exp(alny)
      
      return
      end
! --------------------------------------- y_bssa14_no_site_amps --------------------------------
! --------------------------------------------------------- bssa14_gm_sub4y
      subroutine bssa14_gm_sub4y(
     :        per,
     :        m, rjb, r, v30, mech, iregion, z1, irelation,     
     :        pga4nl,
     :        c, amref, 
     :        rref, delta_c3,
     :        e, amh,
     :        clin, vclin, vref,
     :        f1, f3, f4, f5,
     :        f6, f7,
     :        phi1, phi2, tau1, tau2,
     :        r1, r2, delta_phiR, 
     :        v1, v2, delta_phiV,
     :        y, fe, fpb, fp, fsb, fs, 
     :        phi, tau, sigma, 
     :        gspread,     
     :        amp_lin, amp_nl, amp_total)     
     
 
!Input arguments:

!          per, m, rjb, 
!          mech, v30, 
!          e_gmpe, amh_gmpe, c_gmpe, amref_gmpe,
!          rref_gmpe, h_gmpe, 
!          v30ref, 
!          sigt_gmpe,
!          e_pga4nl, amh_pga4nl, c_pga4nl, amref_pga4nl,
!          rref_pga4nl, h_pga4nl,


!Output arguments:

!          y, expsiglny,
!          pga4nl, fm_pga4nl, fd_pga4nl, 
!          fs_lin_pga4nl, fs_nonlin_pga4nl,
!          bnl, amp_lin, amp_nl, amp_total,
!          fm_gmpe, fd_gmpe, fs_gmpe,
!          gspread


! Computes NGAW2 NGA motions for given input variables

! Note: For generality, include a magnitude-dependent anelastic
! term, although the coefficients are 0.0 for the 2007 BA NGA GMPEs

 
! Dates: 02/27/13 - Written by D.M. Boore, based on ngaw2_gm_sub4y
!        04/09/13 - Revise in light of the new report and coefficient table.
!                 - Sediment depth term in terms of delta_z1, which requires the relation
!                   between z1 and Vs30.  There are two relations in the report
!                   (equations (4.9a) and (4.9b)), one for Japan and one for California.
!                   I need to add an input argument to specify which relation to use, but for
!                   now, just skip the sediment term altogether (even if it is included, specifying
!                   Z1<0 will turn it off).
!                 - Return a large value of phi if Rjb> 300
!        04/11/13 - Ignore restriction on Rjb for phi computation
!        07/01/13 - Renamed from bssa13_gm_sub4y
!        01/09/14 - Adjust mu1_vs30 output from units of m to km.

      implicit none
      
      real  per, m, rjb, r, v30, z1,      
     :        pga4nl,
     :        c(1:3), amref, 
     :        rref, delta_c3(0:3),
     :        e(0:6), amh,
     :        clin, vclin, vref,
     :        f1, f2, f3, f4, f5,
     :        f6, f7,
     :        r1, r2, delta_phiR, 
     :        delta_phiV, v1, v2,
     :        phi1, phi2, 
     :        tau1, tau2,
     :        y, fe, fpb, fp, fsb, fs, 
     :        phiM, phiMR, phi, tau, sigma, 
     :        gspread,     
     :        amp_lin, amp_nl, amp_total, bnl,
     :        delta_z1, f_delta_z1,
     :        mu1_vs30
     
      real y_xamp
     
      integer :: mech, iregion, irelation

        

!GET Y FOR GMPE, WITHOUT SITE AMPS:

      call y_bssa14_no_site_amps(
     :           m, r, mech, iregion,
     :           c, amref, 
     :           rref, delta_c3,
     :           e, amh,
     :           y_xamp, fe, fpb, fp, 
     :           gspread)     
      
!Compute site amps

      if (v30 <= vclin) then
        amp_lin  =  (v30/vref)**clin
      else
        amp_lin  =  (vclin/vref)**clin
      end if
        
      f2 = f4 * (exp(f5*(amin1(v30,760.0)-360.0)) -
     :             exp(f5*(760.0-360.0)           )   )
      
      bnl = f2
      
      amp_nl = exp(f1+f2*alog( (pga4nl+f3)/f3 ))

      amp_total = amp_lin * amp_nl
      
      fsb  = alog(amp_total)   ! Natural log!!
      
      if (z1 < 0.0) then
        delta_z1 = 0.0
      else
        if (irelation == 1 .or. irelation == 2) then
          delta_z1 = z1 - mu1_vs30(v30, irelation)/1000.0
        else
          delta_z1 = 0.0
        end if
      end if
      
      if (per < 0.65) then
        f_delta_z1 = 0.0
      else
        if (delta_z1 > f7/f6) then
          f_delta_z1 = f7
        else
          f_delta_z1 = f6*delta_z1
        end if
      end if
      
      fs = fsb + f_delta_z1
      
      y = exp(fs) * y_xamp
       
!Compute phi, tau, and sigma   

      if (m <= 4.5) then
        tau = tau1
      else if (m >= 5.5) then
        tau = tau2
      else 
        tau = tau1+(tau2-tau1)*(m-4.5)
      end if
      
      if (m <= 4.5) then
        phiM = phi1
      else if (m >= 5.5) then
        phiM = phi2
      else 
        phiM = phi1+(phi2-phi1)*(m-4.5)
      end if
      
      if (Rjb <= R1) then
        phiMR = phiM
      else if (Rjb > R2) then
        phiMR=phiM+delta_phiR
      else
        phiMR = phiM+delta_phiR*(alog(Rjb/R1)/alog(R2/R1))
      end if

      if (v30 <= v1) then
        phi = phiMR - delta_phiV
      else if (v30 >= v2) then
        phi = phiMR
      else
        phi = phiMR - delta_phiV*(alog(v2/v30)/alog(v2/v1))
      end if
      
      sigma = sqrt(phi**2 + tau**2)
         
  
      return
      end subroutine bssa14_gm_sub4y
! --------------------------------------------------------- bssa14_gm_sub4y

! --------------------------------------------------------- mu1_vs30
      function mu1_vs30(vs30, irelation)     

! irelation   relation
!       1     California
!       2     Japan

! NOTE: the units of mu1_vs30 are m. not km

! Dates: 05/10/13 - Written by D.M. Boore 

      implicit none
      
      real :: ln_mu1, mu1_vs30, vs30 
     
      integer :: irelation
      
      if (irelation == 1) then 
      
        ln_mu1 = 
     :   (-7.15/4.0)*alog((vs30**4+570.94**4)/(1360.0**4+570.94**4))
        mu1_vs30 = exp(ln_mu1)
        
      else if (irelation == 2) then
      
        ln_mu1 = 
     :   (-5.23/2.0)*alog((vs30**2+412.39**2)/(1360.0**2+412.39**2))
        mu1_vs30 = exp(ln_mu1)      
      
      else
      
        ! No relation, return negative value
        mu1_vs30 = -9.9
      
      end if
      

      return
      
      end function mu1_vs30
! --------------------------------------------------------- mu1_vs30
! ------------------------------------------------------------------ skipcmnt
      subroutine skipcmnt(nu, comment, ncomments)

! Skip text comments in the file attached to unit nu, but save skipped 
! comments in character array comment.  Skip at least one line, and more as 
! long as the lines are preceded by "|" or "!".

! Dates: 04/16/01 - Written by D. Boore
!        12/07/01 - Added check for eof
!        11/04/03 - Use trim_c to trim possible leading blank
!        02/03/07 - Initialize comments to blank
!        04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)

      character comment(*)*(*), buf*80

      ncomments = 0
100   buf = ' '
      read (nu,'(a)',end=999) buf
      call trim_c(buf,nc_buf)
      if (buf(1:1) .eq.'!' .or. buf(1:1) .eq.'|' .or. 
     :                     ncomments + 1 .eq. 1) then
        ncomments = ncomments + 1
        comment(ncomments) = ' '
        comment(ncomments) = buf(1:nc_buf)
        goto 100
      else 
        backspace nu
      end if

999   continue
 
      return
      end
! ------------------------------------------------------------------ skipcmnt

! --------------------------- BEGIN GET_LUN ----------------
      subroutine get_lun(lun)

! Finds a logical unit number not in use; returns
! -1 if it cannot find one.

! Dates -- 05/19/98 - Written by D. Boore, following
!                     Larry Baker's suggestion
!        04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)

      logical isopen
      do i = 99,10,-1
        inquire (unit=i, opened=isopen)
        if(.not.isopen) then
          lun = i
          return
        end if
      end do
      lun = -1

      return
      end
! --------------------------- END GET_LUN ----------------
     
      subroutine lin_interp(x, y, n, j, x_intrp, y_intrp)
      
* Computes linearly interpolated value of y

* Values out of range are assigned end values

* Dates: 03/16/05 - Written by D. Boore
*        07/24/05 - Added index j for end cases

      real x(*), y(*), slope, x_intrp, y_intrp
      integer j, n
      
      if (x_intrp .le. x(1)) then
        j = 1
        y_intrp = y(1)
        return
      end if
      
      if (x_intrp .ge. x(n)) then
        j = n
        y_intrp = y(n)
        return
      end if          
      
      call locate(x,n,x_intrp,j)
      
      slope = (y(j+1) - y(j))/(x(j+1)-x(j))
      y_intrp = y(j) + slope*(x_intrp - x(j))
      
      return
      end
      


!* --------------------- BEGIN LOCATE -----------------
      SUBROUTINE locate(xx,n,x,j)
      
! Comments added by D. Boore on 26feb2010:
!  finds j such that xx(j) < x <= xx(j+1)
!  EXCEPT if x = xx(1), then j = 1 (logically it would be 0 from
!  the above relation, but this is the same returned value of j
!  for a point out of range).
!  Also, if x = xx(n), j = n-1, which is OK
!  Note that j = 0 or j = n indicates that x is out of range.
!
! See the program test_locate.for to test this routine.

! Dates: 04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)

      INTEGER j,n
      REAL x,xx(n)
      INTEGER jl,jm,ju
      jl=0
      ju=n+1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      goto 10
      endif
      if(x.eq.xx(1))then
        j=1
      else if(x.eq.xx(n))then
        j=n-1
      else
        j=jl
      endif
      return
      END
!* --------------------- END LOCATE -----------------

! ---------------------- BEGIN SKIP -------------------
      subroutine SKIP(lunit, nlines)
        if (nlines .lt. 1) then
          return
        else
          do i = 1, nlines
             read(lunit, *)
          end do
          return
        end if
      end
! ---------------------- END SKIP -------------------
! --------------------- BEGIN UPSTR ----------------------------------
      Subroutine UPSTR ( text )
! Converts character string in TEXT to uppercase
! Dates: 03/12/96 - Written by Larry Baker
!        04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)

!
      Implicit   None
!
      Character  text*(*)
!
      Integer    j
      Character  ch
!
      Do 1000 j = 1,LEN(text)
         ch = text(j:j)
         If ( LGE(ch,'a') .and. LLE(ch,'z') ) Then
            text(j:j) = CHAR ( ICHAR(ch) - ICHAR('a') + ICHAR('A') )
         End If
 1000    Continue
!
      Return
      End
! --------------------- END UPSTR ----------------------------------
! --------------------------- BEGIN TRIM_C -----------------------
      subroutine trim_c(cstr, nchar)

! strips leading and trailing blanks from cstr, returning the
! result in cstr, which is now nchar characters in length

! Strip off tabs also.

! Here is a sample use in constructing a column header, filled out with 
! periods:

!* Read idtag:
!        idtag = ' '
!        read(nu_in, '(1x,a)') idtag
!        call trim_c(idtag, nc_id)
!* Set up the column headings:
!        colhead = ' '
!        colhead = idtag(1:nc_id)//'......' ! nc_id + 6 > length of colhead

! Dates: 12/23/97 - written by D. Boore
!        12/08/00 - pad with trailing blanks.  Otherwise some comparisons
!                   of the trimmed character string with other strings
!                   can be in error because the original characters are left
!                   behind after shifting.  For example, here is a string
!                   before and after shifting, using the old version:
!                      col:12345
!                           MTWH  before
!                          MTWHH  after (but nc = 4).
!        03/21/01 - Check for a zero length input string
!        11/09/01 - Change check for zero length string to check for all blanks
!        10/19/09 - Strip off tabs
!        04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)

      character cstr*(*)

      if(cstr .eq. ' ') then
        nchar = 0
        return
      end if

      nend = len(cstr)

! Replace tabs with blanks:

      do i = 1, nend
        if(ichar(cstr(i:i)) .eq. 9) then
           cstr(i:i) = ' '
        end if
      end do



!      if(nend .eq. 0) then
!        nchar = 0
!        return
!      end if

      do i = nend, 1, -1
        if (cstr(i:i) .ne. ' ') then
           nchar2 = i
           goto 10
        end if
      end do

10    continue

      do j = 1, nchar2
        if (cstr(j:j) .ne. ' ') then
          nchar1 = j
          goto 20
        end if
      end do

20    continue
   
      nchar = nchar2 - nchar1 + 1
      cstr(1:nchar) = cstr(nchar1: nchar2)
      if (nchar .lt. nend) then
        do i = nchar+1, nend
          cstr(i:i) = ' '
        end do
      end if

      return
      end
! --------------------------- END TRIM_C -----------------------