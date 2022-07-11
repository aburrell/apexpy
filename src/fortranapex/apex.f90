! apex.f90

module apexmodule

  implicit none

  real(8), parameter       :: pi=3.14159265358979323846D0
  real(8), parameter       :: dtor=pi / 180D0, rtod=180D0 / pi, pid2=pi / 2D0, twopi=2D0 * pi
  real(8), parameter       :: Req=6378.1370D0, eps=1.D0 / 298.257223563D0
  real(8), parameter       :: Re=Req * (1 - eps / 3D0), ecc2=eps * (2 - eps)
  real(8), parameter       :: missing=- 9999E0

contains

  real(8) function fint(x1, x2, x3, y1, y2, y3, xfit)

          ! Second degree interpolation used by FNDAPX
          ! INPUTS:
          ! X1   = point 1 ordinate value
          ! X2   = point 2 ordinate value
          ! X3   = point 3 ordinate value
          ! Y1   = point 1 abscissa value
          ! Y2   = point 2 abscissa value
          ! Y3   = point 3 abscissa value
          ! XFIT = ordinate value to fit
          ! RETURNS:
          ! YFIT = abscissa value corresponding to XFIT
          !
          ! MODIFICATIONS:
          ! Apr 2004: Change from subroutine to function, rename variables and
          ! add intermediates which are otherwise calculated twice


    real(8), intent(in)  :: x1, x2, x3, y1, y2, y3, xfit
    real(8)              :: x12, x13, x23, xf1, xf2, xf3

    x12 = x1 - x2
    x13 = x1 - x3
    x23 = x2 - x3
    xf1 = xfit - x1
    xf2 = xfit - x2
    xf3 = xfit - x3

    fint = (y1 * x23 * xf2 * xf3 - y2 * x13 * xf1 * xf3 + y3 * x12 * xf1 * xf2) / (x12 * x13 * x23)
  end function fint

end module apexmodule

module dipole
  implicit none
  real(8)                  :: colat, elon, vp, ctp, stp
end module dipole

module fldcomd
  implicit none
  real(8)                  :: bx, by, bz, bb
end module fldcomd

module apxin
  implicit none
  real(8)                  :: yapx(3, 3)
end module apxin

module itra
  implicit none
  integer(4)               :: nstp
  real(8)                  :: sgn, ds
  real(8)                  :: y(3), yold(3), yp(3, 4)
end module itra



subroutine apex(date, igrffilein, dlat, dlon, alt, a, alat, alon, bmag, xmag, ymag, zmag, v)

  ! Calculate apex radius, latitude, longitude; and magnetic field and
  ! scalar magnetic potential.
  !
  ! INPUTS:
  ! DATE = Year and fraction (1990.0 = 1990 January 1, 0 UT)
  ! DLAT = Geodetic latitude in degrees
  ! DLON = Geodetic longitude in degrees
  ! ALT = Altitude in km
  !
  ! RETURNS:
  ! A    = (Apex height + REQ)/REQ, where REQ = equatorial Earth radius.
  ! A is analogous to the L value in invariant coordinates.
  ! ALAT = Apex latitude in degrees (negative in S. magnetic hemisphere)
  ! ALON = Apex longitude (geomagnetic longitude of apex) in degrees
  ! BMAG = geomagnetic field magnitude (nT)
  ! XMAG = geomagnetic field component (nT): north
  ! YMAG = geomagnetic field component (nT): east
  ! ZMAG = geomagnetic field component (nT): downward
  ! V    = geomagnetic potential (T-m)
  !
  ! COMMON BLOCKS:
  ! COMMON /DIPOLE/ COLAT,ELON,VP,CTP,STP
  !
  ! DIPOLE has IGRF variables obtained from routines in magfld.f:
  ! COLAT = Geocentric colatitude of geomagnetic dipole north pole (deg)
  ! ELON  = East longitude of geomagnetic dipole north pole (deg)
  ! VP    = Magnitude (T-m) of dipole component of magnetic potential at
  ! geomagnetic pole and geocentric radius of 6371.0088 km
  ! CTP   = cosine of COLAT
  ! STP   = sine   of COLAT
  !
  ! ------------------------------------------------------------------------------
  ! HISTORY:
  ! Aug 1994: First version completed on the 22nd by A.D. Richmond.
  ! May 1999: Revise DS calculation in LINAPX to avoid divide by zero.
  ! Apr 2004: - Change definition of earth's equatorial radius (REQ)
  ! from the IAU-1966 spheroid (6378.160 km) to the WGS-1984
  ! spheroid (6378.137 km); see description below.
  ! - Revise comments toward a consistent format so they are
  ! easy to read.
  ! - Replace computed GO TO in ITRACE with IF blocks.
  ! - Refine FNDAPX to insure |Bdown/Btot| < 1.E-6 at apex
  ! Nov 2009: Change definition of earth's mean radius (RE) from 6371.2
  ! to the WGS84 value (6371.0088), by J.T. Emmert, NRL
  ! Feb 2021: Modified by Ashton Reimer to pass IGRF coefficients file
  ! to the COFRM subroutine call.
  ! Jul 2022: Revised to fortran 90 standards by L. Lamarche, SRI International.
  !
  ! ------------------------------------------------------------------------------
  ! Reference Spheroid Change                March 2004
  !
  ! Apex geomagnetic coordinates are based on the International Reference
  ! Geomagnetic Field (IGRF) which involves the earth's shape when converting
  ! geographic latitude and altitude to geocentric coordinates.  For this
  ! purpose, the earth is assumed to be an ellipsoid which is fatter at the
  ! equator than the poles.  The World Geodetic System 1984 spheroid
  ! (WGS-1984) is recommended in the recent release of IGRF-9 because it is
  ! used to position current satellite magnetic data (EOS Volume 84 Number 46
  ! November 18 2003).  This differs from previous IGRF releases which favored
  ! the International Astronomical Union 1966 spheroid (IAU-1966) so the Apex
  ! program conversion from geographic to geocentric coordinates in subroutine
  ! CONVRT of file magfld.f has been revised from the IAU-1966 spheroid to the
  ! WGS-1984 spheroid.
  !
  ! The spheroid used to prepare earlier IGRF releases is not always known but
  ! changing spheroids now produces differences at the earth's surface less
  ! than 1 nT, significantly less than other error sources: viz., 9 nT RMS
  ! error for older measurements reporting 1 nT accuracy, 200 nT in the
  ! vicinity of magnetized rocks, or as much as 1000 nT during and after a
  ! geomagnetic storm (www.ngdc.noaa.gov/IAGA/vmod/index.html).
  !
  ! The elliptical shape is characterized in subroutine CONVRT by eccentricity
  ! (e) which is related to the the earth's equatorial radius (a) and polar
  ! radius (b) by
  !
  ! e**2  = 1 - (b/a)**2     (1)
  !
  ! This term is part of an eighth order Lagrange expansion formula (Astron.
  ! J.  Vol. 66, p. 15-16, 1961) designed to give eight digit conversion
  ! accuracy.  The following table summarizes the relevant spheroids:
  !
  ! a           b           e**2          Source
  ! ----------- -----------  -----------  --------------
  ! -           -        0.006722670  Astron J. 1961
  ! 6378.160 km 6356.775 km  0.006701642  IAU-1966
  ! 6378.137 km 6356.752 km  0.006694478  WGS-1984
  !
  ! The previous formulation in CONVRT used the oblateness factor (epsilon),
  ! a surrogate for eccentricity, where
  !
  ! e**2 = 2*epsilon - epsilon**2
  !
  ! with oblateness revised from the 1961 paper's suggested 1/297 to 1/298.25
  ! in accordance with the IAU-1966 spheroid.  Now CONVRT is reformulated to
  ! use equation 1 with the WGS-1984 spheroid's major axis (a) and minor axis
  ! (b) for which epsilon is approximately 1/298.2528.
  !
  ! In addition to earth's equatorial radius and polar radius, the reference
  ! radius (Re) of 6371.2 km is explicit in the IGRF formula for magnetic
  ! potential and implicit in the derived magnetic field coefficients.  The
  ! reference radius has not changed in IGRF releases.
  !
  ! ------------------------------------------------------------------------------

  use dipole
  use apexmodule

  implicit none


  real(8), intent(in)      :: date
  character(len=1000), intent(in)   :: igrffilein
  real(8), intent(in)      :: dlat, dlon, alt
  real(8)      :: a, alon, bmag, xmag, ymag, zmag, v
  real(8)      :: alat, bx, by, bz, x, y, z
  real(8)      :: clatp, polon, vpol

  call cofrm(date, igrffilein)
  call dypol(clatp, polon, vpol)
  colat = clatp
  ctp = cos(clatp * dtor)
  stp = sqrt(1 - CTP * CTP)
  elon = polon
  vp = vpol
  call linapx(dlat, dlon, alt, a, alat, alon, xmag, ymag, zmag, bmag)
  xmag = xmag * 1.e5
  ymag = ymag * 1.e5
  zmag = zmag * 1.e5
  bmag = bmag * 1.e5
  call gd2cart(dlat, dlon, alt, x, y, z)
  call feldg(3, x / Re, y / Re, z / Re, bx, by, bz, v)
  return
end subroutine apex


subroutine linapx(gdlat, glon, alt, a, alat, alon, xmag, ymag, zmag, f)

! Transform geographic coordinates to Apex coordinates.
!
! INPUTS:
! GDLAT = Latitude  (degrees, positive northward)
! GLON  = Longitude (degrees, positive eastward)
! ALT   = Height of starting point (km above mean sea level)
!
! OUTPUTS:
! A     = (Apex height + REQ)/REQ, where REQ = equatorial Earth radius.
! A is analogous to the L value in invariant coordinates.
! ALAT  = Apex Lat. (deg)
! ALON  = Apex Lon. (deg)
! XMAG  = Geomagnetic field component (gauss): north
! YMAG  = Geomagnetic field component (gauss): east
! ZMAG  = Geomagnetic field component (gauss): down
! F     = Geomagnetic field magnitude (gauss)
!
! Trace the geomagnetic field line from the given location to find the
! apex of the field line.  Before starting iterations to trace along
! the field line: (1) Establish a step size (DS, arc length in km)
! based on the geomagnetic dipole latitude; (2) determine the step
! direction from the sign of the vertical component of the geomagnetic
! field; and (3) convert to geocentric cartesian coordinates.  Each
! iteration increments a step count (NSTP) and calls ITRACE to move
! along the the field line until reaching the iteration count limit
! (MAXS) or passing the apex (IAPX=2) and then calling FNDAPX to
! determine the apex location from the last three step locations
! (YAPX); however, if reaching the iteration limit, apex coordinates
! are calculated by DIPAPX which assumes a simplified dipole field.
!
! COMMON BLOCKS:
! COMMON /APXIN/   YAPX(3,3)
! COMMON /DIPOLE/  COLAT,ELON,VP,CTP,STP
! COMMON /FLDCOMD/ BX, BY, BZ, BB
! COMMON /ITRA/    NSTP, Y(3), YOLD(3), SGN, DS
!
! APXIN has step locations determined in ITRACE:
! YAPX  = Matrix of cartesian coordinates (loaded columnwise) of the
! three points about the apex.  Set in subroutine ITRACE.
!
! DIPOLE has IGRF variables obtained from routines in magfld.f:
! COLAT = Geocentric colatitude of geomagnetic dipole north pole (deg)
! ELON  = East longitude of geomagnetic dipole north pole (deg)
! VP    = Magnitude (T-m) of dipole component of magnetic potential at
! geomagnetic pole and geocentric radius of 6371.0088 km
! CTP   = cosine of COLAT
! STP   = sine   of COLAT
!
! FLDCOMD has geomagnetic field at current trace point:
! BX    = X component (Gauss)
! BY    = Y component (Gauss)
! BZ    = Z component (Gauss)
! BB    = Magnitude   (Gauss)
!
! ITRA has field line tracing variables determined in LINAPX:
! NSTP  = Step count.
! Y     = Array containing current tracing point cartesian coordinates.
! YOLD  = Array containing previous tracing point cartesian coordinates.
! SGN   = Determines direction of trace.
! DS    = Step size (arc length in km).
!
! REFERENCES:
! Stassinopoulos E. G. , Mead Gilbert D., X-841-72-17 (1971) GSFC,
! Greenbelt, Maryland
!
! EXTERNALS:
! GD2CART = Convert geodetic to geocentric cartesian coordinates (in magfld.f)
! CONVRT  = Convert geodetic to geocentric cylindrical or geocentric spherical
! and back (in magfld.f).
! FELDG   = Obtain IGRF magnetic field components (in magfld.f).
! ITRACE  = Follow a geomagnetic field line
! DIPAPX  = Compute apex coordinates assuming a geomagnetic dipole field
! FNDAPX  = Compute apex coordinates from the last three traced field line points
!
! ------------------------------------------------------------------------------
! HISTORY:
! Oct 1973: Initial version completed on the 29th by Wally Clark, NOAA
! ERL Lab.
! Feb 1988: Revised on the 1st by Harsh Anand Passi, NCAR.
! Aug 1994: Revision by A. D. Richmond, NCAR.
! Nov 2009: Change definition of earth's mean radius (RE) from 6371.2
! to the WGS84 value (6371.0088), by J.T. Emmert, NRL.
! Jul 2022: Revised to fortran 90 standards by L. Lamarche, SRI International.



  use fldcomd
  use apxin
  use dipole
  use itra
  use apexmodule

  implicit none

  real(8)               :: gdlat, glon, alt, a, alat, alon, xmag, ymag, zmag, f
  real(8)               :: bnrth, beast, bdown, babs
  real(8)               :: cgml2, gclat, ht, r, rho
  real(8)               :: singml, xlat, xlon
  integer(4)            :: i, j, iapx
  integer(4), parameter :: maxs = 200

  ! Set step size based on the goemagnetic dipole latitude of the starting point
  call convrt(2, gdlat, alt, gclat, r)
  singml = ctp * sin(gclat * dtor) + stp * cos(gclat * dtor) * cos((glon - elon) * dtor)
  ! May 1999: avoid possible divide by zero (when SINGML = 1.): the old version limited DS to its value at 60 deg GM latitude with: DS = .06*R/(1.-SINGML*SINGML) - 370. IF (DS .GT. 1186.) DS = 1186.
  cgml2 = amax1(0.25, 1.- singml * singml)
  ds = 0.06 * r / cgml2 - 370.

  ! Initialize yapx array
  do j = 1, 3
    do i = 1, 3
      yapx(i, j) = 0.
    end do
  end do

  ! Convert from geodetic to earth centered cartesian coordinates
  call gd2cart(gdlat, glon, alt, y(1), y(2), y(3))
  nstp = 0

  ! Get magnetic field components to determine the direction for tracing the field line
  call feldg(1, gdlat, glon, alt, xmag, ymag, zmag, f)
  sgn = sign(1.D0, - zmag)

  ! Use cartesian coordinates to get magnetic field components (from which gradients steer the tracing)
  iapx = 1
  do while (iapx .EQ. 1)
    call feldg(2, y(1) / Re, y(2) / Re, y(3) / Re, bx, by, bz, bb)
    nstp = nstp + 1
    if (nstp .ge. maxs) exit
    call itrace(iapx)
  end do
  if (nstp < maxs) then
    call fndapx(alt, zmag, a, alat, alon)
  else
    rho = sqrt(y(1) ** 2 + y(2) ** 2)
    call convrt(3, xlat, ht, rho, y(3))
    xlon = rtod * atan2(y(2), y(1))
    call feldg(1, xlat, xlon, ht, bnrth, beast, bdown, babs)
    call dipapx(xlat, xlon, ht, bnrth, beast, bdown, a, alon)
    alat = - sgn * rtod * acos(sqrt(1./ a))
  end if

  return

end subroutine linapx


subroutine itrace(iapx)

! Follow a geomagnetic field line until passing its apex
!
! INPUTS:
! (all are in common blocks)
! OUTPUTS:
! IAPX = 2 (when apex passed) or 1 (not)
!
! This uses the 4-point Adams formula after initialization.
! First 7 iterations advance point by 3 steps.
!
! COMMON BLOCKS:
! COMMON /APXIN/   YAPX(3,3)
! COMMON /FLDCOMD/ BX, BY, BZ, BB
! COMMON /ITRA/    NSTP, Y(3), YOLD(3), SGN, DS
!
! APXIN has step locations determined in ITRACE:
! YAPX  = Matrix of cartesian coordinates (loaded columnwise) of the
! three points about the apex.  Set in subroutine ITRACE.
!
! FLDCOMD has geomagnetic field at current trace point:
! BX    = X component (Gauss)
! BY    = Y component (Gauss)
! BZ    = Z component (Gauss)
! BB    = Magnitude   (Gauss)
!
! ITRA has field line tracing variables determined in LINAPX:
! NSTP  = Step count.
! Y     = Array containing current tracing point cartesian coordinates.
! YOLD  = Array containing previous tracing point cartesian coordinates.
! SGN   = Determines direction of trace.
! DS    = Step size (arc length in km).
!
! REFERENCES:
! Stassinopoulos E. G. , Mead Gilbert D., X-841-72-17 (1971) GSFC,
! Greenbelt, Maryland
! ------------------------------------------------------------------------------
! HISTORY:
! Oct 1973: Initial version completed on the 29th by W. Clark, NOAA ERL
! Laboratory.
! Feb 1988: Revised by H. Passi, NCAR.
! Apr 2004: Replace computed GO TO with IF blocks because some compilers
! are threatening to remove this old feature
! Jul 2022: Revised to fortran 90 standards by L. Lamarche, SRI International.



  use apxin
  use fldcomd
  use itra

  implicit none

  integer(4), intent(out)  :: iapx
  real(8)                  :: d12, d2, d24, d6, rc, rp
  integer(4)               :: i, j


  iapx = 1

  ! Cartesian component magnetic field (partial) derivitives steer the trace
  yp(1, 4) = sgn * bx / bb
  yp(2, 4) = sgn * by / bb
  yp(3, 4) = sgn * bz / Bb

  d2 = ds / 2.
  d6 = ds / 6.
  d12 = ds / 12.
  d24 = ds / 24.

  if (nstp .le. 7) then
    do i = 1, 3
      if (nstp == 1) then
        ! d2 = ds/2.
        ! d6 = ds/6.
        ! d12 = ds/12.
        ! d24 = ds/24.
        yp(i, 1) = yp(i, 4)
        yold(i) = y(i)
        yapx(i, 1) = y(i)
        y(i) = yold(i) + ds * yp(i, 1)

      else if (nstp == 2) then
        yp(i, 2) = yp(i, 4)
        y(i) = yold(i) + d2 * (yp(i, 2) + yp(i, 1))

      else if (nstp == 3) then
        y(i) = yold(i) + d6 * (2.* yp(i, 4) + yp(i, 2) + 3.* yp(i, 1))

      else if (nstp == 4) then
        yp(i, 2) = yp(i, 4)
        yapx(i, 2) = y(i)
        yold(i) = y(i)
        y(i) = yold(i) + d2 * (3.* yp(i, 2) - yp(i, 1))

      else if (nstp == 5) then
        y(i) = yold(i) + d12 * (5.* yp(i, 4) + 8.* yp(i, 2) - yp(i, 1))

      else if (nstp == 6) then
        yp(i, 3) = yp(i, 4)
        yold(i) = y(i)
        yapx(i, 3) = y(i)
        y(i) = yold(i) + d12 * (23.* yp(i, 3) - 16.* yp(i, 2) + 5.* yp(i, 1))

      else if (nstp == 7) then
        yapx(i, 1) = yapx(i, 2)
        yapx(i, 2) = yapx(i, 3)
        y(i) = yold(i) + d24 * (9.* yp(i, 4) + 19.* yp(i, 3) - 5.* yp(i, 2) + yp(i, 1))
        yapx(i, 3) = y(i)

      end if
    end do

    if (nstp == 6 .or. nstp == 7) then
      rc = sqrt(yapx(1, 3) ** 2 + yapx(2, 3) ** 2 + yapx(3, 3) ** 2)
      rp = sqrt(yapx(1, 2) ** 2 + yapx(2, 2) ** 2 + yapx(3, 2) ** 2)
      if (rc < rp) iapx = 2
    end if

  else
    do i = 1, 3
      yapx(i, 1) = yapx(i, 2)
      yapx(i, 2) = y(i)
      yold(i) = y(i)
      y(i) = yold(i) + d24 * (55.* yp(i, 4) - 59.* yp(i, 3) + 37.* yp(i, 2) - 9.* yp(i, 1))
      yapx(i, 3) = y(i)
      do j = 1, 3
        yp(i, j) = yp(i, j + 1)
      end do
    end do
    rc = sqrt(y(1) ** 2 + y(2) ** 2 + y(3) ** 2)
    rp = sqrt(yold(1) ** 2 + yold(2) ** 2 + yold(3) ** 2)
    if (rc < rp) iapx = 2

  end if

  return
end subroutine itrace


subroutine fndapx(alt, zmag, a, alat, alon)

! Find apex coordinates once tracing (in subroutine ITRACE) has
! signalled that the apex has been passed.
! INPUTS:
! ALT  = Altitude of starting point
! ZMAG = Downward component of geomagnetic field at starting point
! OUTPUT
! A    = Apex radius, defined as (Apex height + Req)/Req, where
! Req = equatorial Earth radius.
! A is analogous to the L value in invariant coordinates.
! ALAT = Apex Lat. (deg)
! ALON = Apex Lon. (deg)
!
! COMMON BLOCKS:
! COMMON /APXIN/  YAPX(3,3)
! COMMON /DIPOLE/ COLAT,ELON,VP,CTP,STP
!
! APXIN has step locations determined in ITRACE:
! YAPX  = Matrix of cartesian coordinates (loaded columnwise) of the
! three points about the apex.  Set in subroutine ITRACE.
!
! DIPOLE has IGRF variables obtained from routines in magfld.f:
! COLAT = Geocentric colatitude of geomagnetic dipole north pole (deg)
! ELON  = East longitude of geomagnetic dipole north pole (deg)
! VP    = Magnitude (T-m) of dipole component of magnetic potential at
! geomagnetic pole and geocentric radius of 6371.0088 km
! CTP   = cosine of COLAT
! STP   = sine   of COLAT
!
! EXTERNALS:
! FINT = Second degree interpolation routine
! ------------------------------------------------------------------------------
! HISTORY:
! Oct 1973: Initial version completed on the 23rd by Clark, W., NOAA
! Boulder.
! Aug 1994: Revision on the 3rd by A.D. Richmond, NCAR
! Apr 2004: Repair problem noted by Dan Weimer where the apex location
! produced by FINT may still have a non-zero vertical magnetic
! field component.
! Jul 2022: Revised to fortran 90 standards by L. Lamarche, SRI International.


  use apxin
  use dipole
  use apexmodule

  implicit none

  real(8)       :: alt, zmag, a, alat, alon
  real(8)       :: abdob, ang, ba, bda, bea, bna, hta, rho
  real(8)       :: r, rasq, sang, cang, ste, cte, stfcpa, stfspa, xlon
  real(8)       :: gdlt, gdln, ht, bn, be, bmag
  integer(4)    :: i, nitr
  real(8)       :: bd(3), y(3)

  ! Get geodetic height and vertical (downward) compontent of the magnetic field at least three points found by ITRACE
  do i = 1, 3
    rho = sqrt(yapx(1, i) ** 2 + yapx(2, i) ** 2)
    call convrt(3, gdlt, ht, rho, yapx(3, i))
    gdln = rtod * atan(yapx(2, i), yapx(1, i))
    call feldg(1, gdlt, gdln, ht, bn, be, bd(i), bmag)
  end do

  ! Interpolate to where Bdown=0 to find cartesian coordinates at dip equatior
  nitr = 0
  do while (nitr < 4)    ! 4 was chosen because tests rarely required 2 iterations
    y(1) = fint(bd(1), bd(2), bd(3), yapx(1, 1), yapx(1, 2), yapx(1, 3), 0.D0)
    y(2) = fint(bd(1), bd(2), bd(3), yapx(2, 1), yapx(2, 2), yapx(2, 3), 0.D0)
    y(3) = fint(bd(1), bd(2), bd(3), yapx(3, 1), yapx(3, 2), yapx(3, 3), 0.D0)

          ! Insure negligible Bdown or
          !
          ! |Bdown/Btot| < 2.E-6
          !
          ! For instance, Bdown must be less than 0.1 nT at low altitudes where
          ! Btot ~ 50000 nT.  This ratio can be exceeded when interpolation is
          ! not accurate; i.e., when the middle of the three points interpolated
          ! is too far from the dip equator.  The three points were initially
          ! defined with equal spacing by ITRACE, so replacing point 2 with the
          ! most recently fit location will reduce the interpolation span.

    rho = sqrt(y(1) ** 2 + y(2) ** 2)
    gdln = rtod * atan2(y(2), y(1))
    call convrt(3, gdlt, hta, rho, y(3))
    call feldg(1, gdlt, gdln, hta, bna, bea, bda, ba)
    abdob = abs(bda / ba)

    if (abdob > 2.e-6) then
      nitr = nitr + 1
      yapx(1, 2) = y(1)
      yapx(2, 2) = y(2)
      yapx(3, 2) = y(3)
      bd(2) = bda
    else
      exit
    end if
  end do

  if (abdob > 2.e-6) then
    write(0, '("APEX: Imprecise fit of apex: |Bdown/B| "(1PE7.1))') abdob
  end if


  ! Ensure altitude of the Apex is at least the initial altitude when
  ! defining the Apex radius then use it to define the Apex latitude whose
  ! hemisphere (sign) is inferred from the sign of the dip angle at the
  ! starting point

  a = (Req + amax1(alt, hta)) / Req
  if (a < 1.) then
    ! write(0, '("APEX: A cannot be less than 1; A, REQ, HTA "(1P3))') a, Req, hta
    write(0, *) "APEX: A cannot be less than 1; A, REQ, HTA ", a, Req, hta
    call exit(1)
  end if
  rasq = acos(sqrt(1./ a)) * rtod
  alat = sign(rasq, zmag)

  ! ALON is the dipole longitude of the apex and is defined using
  ! spherical coordinates where
  ! GP   = geographic pole.
  ! GM   = geomagnetic pole (colatitude COLAT, east longitude ELON).
  ! XLON = longitude of apex.
  ! TE   = colatitude of apex.
  ! ANG  = longitude angle from GM to apex.
  ! TP   = colatitude of GM.
  ! TF   = arc length between GM and apex.
  ! PA   = ALON be geomagnetic longitude, i.e., Pi minus angle measured
  ! counterclockwise from arc GM-apex to arc GM-GP.
  ! then, spherical-trigonometry formulas for the functions of the angles
  ! are as shown below.  Notation uses C=cos, S=sin and STFCPA = sin(TF) * cos(PA),
  ! STFSPA = sin(TF) * sin(PA)

  xlon = atan2(y(2), y(1))
  ang = xlon - elon * dtor
  cang = cos(ang)
  sang = sin(ang)
  r = sqrt(y(1) ** 2 + y(2) ** 2 + y(3) ** 2)
  cte = y(3) / r
  ste = sqrt(1.- cte * cte)
  stfcpa = ste * ctp * cang - cte * stp
  stfspa = sang * ste
  alon = atan2(stfspa, stfcpa) * rtod

  return
end subroutine fndapx


subroutine dipapx(gdlat, gdlon, alt, bnorth, beast, bdown, a, alon)

! Compute A, ALON from local magnetic field using dipole and spherical
! approximation.
!
! INPUTS:
! GDLAT  = geodetic latitude, degrees
! GDLON  = geodetic longitude, degrees
! ALT    = altitude, km
! BNORTH = geodetic northward magnetic field component (any units)
! BEAST  = eastward magnetic field component
! BDOWN  = geodetic downward magnetic field component
! OUTPUTS:
! A      = apex radius, 1 + h_A/R_eq
! ALON   = apex longitude, degrees
!
! Use spherical coordinates and define:
! GP    = geographic pole.
! GM    = geomagnetic pole (colatitude COLAT, east longitude ELON).
! G     = point at GDLAT,GDLON.
! E     = point on sphere below apex of dipolar field line passing
! through G.
! TD    = dipole colatitude of point G, found by applying dipole
! formula for dip angle to actual dip angle.
! B     = Pi plus local declination angle.  B is in the direction
! from G to E.
! TG    = colatitude of G.
! ANG   = longitude angle from GM to G.
! TE    = colatitude of E.
! TP    = colatitude of GM.
! A     = longitude angle from G to E.
! APANG = A + ANG
! PA    = geomagnetic longitude, i.e., Pi minus angle measured
! counterclockwise from arc GM-E to arc GM-GP.
! TF    = arc length between GM and E.
! Then, using notation C=cos, S=sin, COT=cot, spherical-trigonometry
! formulas for the functions of the angles are as shown below.  Note:
! STFCPA = sin(TF) * cos(PA)
! STFSPA = sin(TF) * sin(PA)
!
! COMMON BLOCKS:
! COMMON /DIPOLE/ COLAT,ELON,VP,CTP,STP
!
! DIPOLE has IGRF variables obtained from routines in magfld.f:
! COLAT = Geocentric colatitude of geomagnetic dipole north pole (deg)
! ELON  = East longitude of geomagnetic dipole north pole (deg)
! VP    = Magnitude (T-m) of dipole component of magnetic potential at
! geomagnetic pole and geocentric radius of 6371.0088 km
! CTP   = cosine of COLAT
! STP   = sine   of COLAT
! ------------------------------------------------------------------------------
! HISTORY:
! May 1994:  Completed on the 1st by A. D. Richmond
! Nov 2009: Change definition of earth's mean radius (RE) from 6371.2
! to the WGS84 value (6371.0088), by J.T. Emmert, NRL.
! Jul 2022: Revised to fortran 90 standards by L. Lamarche, SRI International.


  use dipole
  use apexmodule

  implicit none

  real(8)        :: gdlat, gdlon, alt, bnorth, beast, bdown, a, alon
  real(8)        :: ang, bhor, ca, cang, capang, cb, cottd, ctd, cte, ctg
  real(8)        :: ha, r, sa, sang, sb, std, ste, stfcpa, stg, sapang, stfspa


  bhor = sqrt(bnorth * bnorth + beast * beast)
  if (bhor == 0.) then
    alon = 0.
    a = 1.E34
    return
  end if

  cottd = bdown * 0.5 / bhor
  std = 1./ sqrt(1.+ cottd * cottd)
  ctd = cottd * std
  sb = - beast / bhor
  cb = - bnorth / bhor
  ctg = sin(gdlat * dtor)
  stg = cos(gdlat * dtor)
  ang = (gdlon - elon) * dtor
  sang = sin(ang)
  cang = cos(ang)
  cte = ctg * std + stg * ctd * cb
  ste = sqrt(1.- cte * cte)
  sa = sb * ctd / ste
  ca = (std * stg - ctd * ctg * cb) / ste
  capang = ca * cang - sa * sang
  sapang = ca * sang + sa * cang
  stfcpa = ste * ctp * capang - cte * stp
  stfspa = sapang * ste
  alon = atan2(stfspa, stfcpa) * rtod
  r = alt + Re
  ha = alt + r * cottd * cottd
  a = 1. + ha / Req

  return
end subroutine dipapx


! real(8) function fint(x1, x2, x3, y1, y2, y3, xfit)
!
! implicit none
!
! real(8), intent(in)  :: x1, x2, x3, y1, y2, y3, xfit
! real(8)              :: x12, x13, x23, xf1, xf2, xf3
!
! x12 = x1-x2
! x13 = x1-x3
! x23 = x2-x3
! xf1 = xfit-x1
! xf2 = xfit-x2
! xf3 = xfit-x3
!
! fint = (y1*x23*xf2*xf3 - y2*x13*xf1*xf3 + y3*x12*xf1*xf2)/(x12*x13*x23)
! return
! end function fint
