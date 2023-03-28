! **********************************************
!
! File Name: magfld.f90
! Authors:
! Translated to Fortran 90 by Leslie Lamarche
! Date: 01/27/2022
!
! ***********************************************


module magcof

  implicit none

  integer(4)       :: nmax1
  real(8)          :: ichg
  real(8)          :: gb(225), gv(225)

end module magcof

module coeffmodule
  ! Req = Equatorial radius of Earth in km (WGS84 value)
  ! eps = flatness of ellipsoidal Earth (WGS84 value)
  ! Re = Mean radius of Earth in km
  ! ecc2 = squared eccentricity of ellipsoidal Earth
  real(8), parameter       :: pi=3.14159265358979323846D0
  real(8), parameter       :: dtor=pi / 180D0, rtod=180D0 / pi, pid2=pi / 2D0, twopi=2D0 * pi
  real(8), parameter       :: Req=6378.1370D0, Rep = 6356.7520D0
  real(8), parameter       :: eps=1.D0 / 298.257223563D0, Re=Req * (1 - eps / 3D0)
  ! real(8), parameter       :: ecc2 = 1D0 - (Rep / Req)
  real(8), parameter       :: ecc2 = eps * (2D0 - eps)
  real(8), parameter       :: missing=- 9999E0
end module coeffmodule


subroutine cofrm(date, filename)

  ! Define the International Geomagnetic Reference Field (IGRF) as a
  ! scalar potential field using a truncated series expansion with
  ! Schmidt semi-normalized associated Legendre functions of degree n and
  ! order m.  The polynomial coefficients are a function of time and are
  ! interpolated between five year epochs or extrapolated at a constant
  ! rate after the last epoch.
  !
  ! INPUTS:
  ! date = yyyy.fraction (UT)
  ! fielname = filename for IGRF coefficient file
  ! OUTPUTS (variables in imported module magcof):
  ! nmax1 = Maximum order of spherical harmonic coefficients used
  ! gb   = Coefficients for magnetic field calculation
  ! gv   = Coefficients for magnetic potential calculation
  ! ichg = Flag indicating when gb,gv have been changed in cofrm
  !
  ! It is fatal to supply a date before the first epoch.  A warning is
  ! issued to Fortran unit 0 (stderr) if date is later than the
  ! recommended limit, five years after the last epoch.
  !
  ! HISTORY (blame):
  ! Apr 1983:  Written by Vincent B. Wickwar (Utah State Univ.) including
  ! secular variation acceleration rate set to zero in case the IGRF
  ! definition includes such second time derivitives.  The maximum degree
  ! (n) defined was 10.
  !
  ! Jun 1986:  Updated coefficients adding Definitive Geomagnetic Reference
  ! Field (DGRF) for 1980 and IGRF for 1985 (EOS Volume 7 Number 24).  The
  ! designation DGRF means coefficients will not change in the future
  ! whereas IGRF coefficients are interim pending incorporation of new
  ! magnetometer data.  Common block MAG was replaced by MAGCOF, thus
  ! removing variables not used in subroutine FELDG.  (Roy Barnes)
  !
  ! Apr 1992 (Barnes):  Added DGRF 1985 and IGRF 1990 as given in EOS
  ! Volume 73 Number 16 April 21 1992.  Other changes were made so future
  ! updates should:  (1) Increment NDGY; (2) Append to EPOCH the next IGRF
  ! year; (3) Append the next DGRF coefficients to G1DIM and H1DIM; and (4)
  ! replace the IGRF initial values (G0, GT) and rates of change indices
  ! (H0, HT).
  !
  ! Apr 1994 (Art Richmond): Computation of GV added, for finding magnetic
  ! potential.
  !
  ! Aug 1995 (Barnes):  Added DGRF for 1990 and IGRF for 1995, which were
  ! obtained by anonymous ftp to geomag.gsfc.nasa.gov (cd pub, mget table*)
  ! as per instructions from Bob Langel (langel@geomag.gsfc.nasa.gov) with
  ! problems reported to baldwin@geomag.gsfc.nasa.gov.
  !
  ! Oct 1995 (Barnes):  Correct error in IGRF-95 G 7 6 and H 8 7 (see email
  ! in folder).  Also found bug whereby coefficients were not being updated
  ! in FELDG when IENTY did not change so ICHG was added to flag date
  ! changes.  Also, a vestigial switch (IS) was removed from COFRM; it was
  ! always zero and involved 3 branch if statements in the main polynomial
  ! construction loop now numbered 200.
  !
  ! Feb 1999 (Barnes):  Explicitly initialize GV(1) in COFRM to avoid the
  ! possibility of compiler or loader options initializing memory to
  ! something else (e.g., indefinite).  Also simplify the algebra in COFRM
  ! with no effect on results.
  !
  ! Mar 1999 (Barnes):  Removed three branch if's from FELDG and changed
  ! statement labels to ascending order.
  !
  ! Jun 1999 (Barnes):  Corrected RTOD definition in GD2CART.
  !
  ! May 2000 (Barnes):  Replace IGRF 1995, add IGRF 2000, and extend the
  ! earlier DGRF's back to 1900.  The coefficients came from an NGDC web
  ! page.  Related documentation is in $APXROOT/docs/igrf.2000.*  where
  ! $APXROOT, defined by 'source envapex', is traditionally ~bozo/apex).
  !
  ! Mar 2004 (Barnes):  Replace 1995 and 2000 coefficients; now both are
  ! DGRF.  Coefficients for 2000 are degree 13 with precision increased to
  ! tenths nT and accommodating this instigated changes:  (1) degree (NMAX)
  ! is now a function of epoch (NMXE) to curtail irrelevant looping over
  ! unused high order terms (n > 10 in epochs before 2000) when calculating
  ! GB; (2) expand coefficients data statement layout for G1D and H1D,
  ! formerly G1DIM and H1DIM; (3) omit secular variation acceleration terms
  ! which were always zero; (4) increase array dimensions in common block
  ! MAGCOF and associated arrays G and H in FELDG; (5) change earth's shape
  ! in CONVRT from the IAU-1966 to the WGS-1984 spheroid; (6) eliminate
  ! reference to 'definitive' in variables in COFRM which were not always
  ! definitive; (7) change G to GB in COFRM s.t. arrays GB and GV in common
  ! block MAGCOF are consistently named in all subroutines; (8) remove
  ! unused constants in all five subroutines.  See EOS Volume 84 Number 46
  ! November 18 2003, www.ngdc.noaa.gov/IAGA/vmod/igrf.html or local files
  ! $APXROOT/docs/igrf.2004.*
  !
  ! Sept. 2005 (Maute): update with IGRF10 from
  ! http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html use script
  ! ~maute/apex.d/apex_update/igrf2f Note that the downloaded file the start
  ! column of the year in the first line has to be before the start of each
  ! number in the same column
  !
  ! Jan. 2010 (Maute) update with IGRF11 (same instructions as Sep. 2005
  ! comment
  !
  ! May 2020 (Achim Morschhauser): Update with routine to read
  ! IGRF coefficients file directly.
  !
  ! Jul 2022 (Lamarche): Revise to fortran 90 standards


    use magcof
    use igrf

    implicit none

    real(8), intent(in)              :: date
    character(len=1000), intent(in)  :: filename
    real(8)                          :: f, f0, t, to5
    integer(4)                       :: i, i1, iy, iy1, m, n, mm, nn, ngh


    ! Do not need to load new coefficients if date has not changed
    ichg = 0
    if (date == datel) then
      return
    else
      datel = date
      ichg = 1
    end if

    ! Read in the IGRF file, if it has not already been loaded
    if (.not. allocated(gyr)) then
      call read_igrf(filename)
    end if
    ngh = nght * nepo

    ! Test the desired output epoch
    if (date < epoch(1)) then
      write(0, '("COFRM:  DATE "(F9.3)" preceeds earliest available "(F6.1))') date, epoch(1)
      call exit(1)
    elseif (date > epoch(nepo) + 5.) then
      write(0, '("COFRM:  DATE "(F9.3)" is after the last recommended for extrapolation "(F6.1))') date, epoch(1) + 5.
      call exit(1)
    end if

    ! Set the date and time index
    do iy = 1, nepo
      if (date < epoch(iy)) exit
    end do
    iy = iy - 1

    ngh = nght * nepo
    nmax1 = int(nmxe(iy))
    ! time = date
    t = date - epoch(iy)
    to5 = t / 5.
    iy1 = iy + 1
    gb(1) = 0.0
    gv(1) = 0.0
    i = 2
    f0 = - 1.0d-5

    do n = 1, nmax1
      f0 = f0 * real(n) / 2.
      f = f0 / sqrt(2.0)
      nn = n + 1
      mm = 1

      if (iy < nepo) then  ! interoplate (m=0 terms)
        gb(i) = (gyr(nn, mm, iy) + (gyr(nn, mm, iy1) - gyr(nn, mm, iy)) * to5) * f0
      end if
      if (iy == nepo) then  ! extrapolate (m=0 terms)
        gb(i) = (gyr(nn, mm, iy) + gt(nn, mm) * t) * f0
      end if
      gv(i) = gb(i) / real(nn)
      i = i + 1
      do m = 1, n
        f = f / sqrt(real(n - m + 1) / real(n + m))
        nn = n + 1
        mm = m + 1
        i1 = i + 1
        if (iy < nepo) then  ! interpolate (m>0 terms)
          gb(i) = (gyr(nn, mm, iy) + (gyr(nn, mm, iy1) - gyr(nn, mm, iy)) * to5) * f
          gb(i1) = (hyr(nn, mm, iy) + (hyr(nn, mm, iy1) - hyr(nn, mm, iy)) * to5) * f
        else                    ! extrapolate (m>0 terms)
          gb(i) = (gyr(nn, mm, iy) + gt(nn, mm) * t) * f
          gb(i1) = (hyr(nn, mm, iy) + ht(nn, mm) * t) * f
        end if
        gv(i) = gb(i) / real(nn)
        gv(i1) = gb(i1) / real(nn)
        i = i + 2
      end do
    end do

    return
end subroutine cofrm


subroutine dypol(colat, elon, vp)

  ! Computes parameters for dipole component of geomagnetic field.
  ! cofrm must be called before calling dypol!
  ! 940504 A. D. Richmond
  !
  ! INPUTS (variables in imported module magcof):
  ! nmax1 = Maximum order of spherical harmonic coefficients used
  ! gb   = Coefficients for magnetic field calculation
  ! gv   = Coefficients for magnetic potential calculation
  ! ichg = Flag indicating when gb,gv have been changed
  !
  ! RETURNS:
  ! colat = Geocentric colatitude of geomagnetic dipole north pole
  ! (deg)
  ! elon  = East longitude of geomagnetic dipole north pole (deg)
  ! vp    = Magnitude, in T.m, of dipole component of magnetic
  ! potential at geomagnetic pole and geocentric radius
  ! of 6371.2 km
  !
  ! IMPORTS:
  ! magcof
  ! coeffmodule


  use magcof
  use coeffmodule

  implicit none

  real(8), intent(out)      :: colat, elon, vp
  real(8)                   :: ctp, gpl, stp

  ! Compute geographic colatitude and logitude of the north pole of earth centered dipole
  gpl = sqrt(gb(2) ** 2 + gb(3) ** 2 + gb(4) ** 2)
  ctp = gb(2) / gpl
  stp = sqrt(1 - ctp * ctp)
  colat = acos(ctp) * rtod
  elon = atan2(gb(4), gb(3)) * rtod

  ! Compute magnitude of magnetic potential at pole, radius Re
  vp = 0.2 * gpl * Re
  ! .2 = 2*(10**-4 T/gauss)*(1000 m/km) (2 comes through f0 in cofrm).

  return
end subroutine dypol


subroutine feldg(ienty, glat, glon, alt, bnrth, beast, bdown, babs)

  ! Compute the DGRF/IGRF field components at the point glat,glon,alt.
  ! cofrm must be called to establish coefficients for correct date
  ! prior to calling feldg.
  !
  ! ienty is an input flag controlling the meaning and direction of the
  ! remaining formal arguments:
  ! ienty = 1
  ! INPUTS:
  ! glat = Latitude of point (deg)
  ! glon = Longitude (east=+) of point (deg)
  ! alt  = Ht of point (km)
  ! RETURNS:
  ! bnrth = north component of field vector (Gauss)
  ! beast = east component of field vector  (Gauss)
  ! bdown = downward component of field vector (Gauss)
  ! babs  = magnitude of field vector (Gauss)
  !
  ! ienty = 2
  ! INPUTS:
  ! glat = X coordinate (in units of earth radii 6371.2 km )
  ! glon = Y coordinate (in units of earth radii 6371.2 km )
  ! alt  = Z coordinate (in units of earth radii 6371.2 km )
  ! RETURNS:
  ! bnrth = X component of field vector (Gauss)
  ! beast = Y component of field vector (Gauss)
  ! bdown = Z component of field vector (Gauss)
  ! babs  = Magnitude of field vector (Gauss)
  !
  ! ienty = 3
  ! INPUTS:
  ! glat = X coordinate (in units of earth radii 6371.2 km )
  ! glon = Y coordinate (in units of earth radii 6371.2 km )
  ! alt  = Z coordinate (in units of earth radii 6371.2 km )
  ! RETURNS:
  ! bnrth = Dummy variable
  ! beast = Dummy variable
  ! bdown = Dummy variable
  ! babs  = Magnetic potential (T.m)
  !
  ! INPUTS (variables in imported module magcof):
  ! nmax1 = Maximum order of spherical harmonic coefficients used
  ! gb    = Coefficients for magnetic field calculation
  ! gv    = Coefficients for magnetic potential calculation
  ! ichg  = Flag indicating when GB,GV have been changed
  !
  ! HISTORY:
  ! Apr 1983: written by Vincent B. Wickwar (Utah State Univ.).
  !
  ! May 1994 (A.D. Richmond): Added magnetic potential calculation
  !
  ! Oct 1995 (Barnes): Added ICHG
  !
  ! Jul 2022 (Lamarche): Revise to fortran 90 standards


  use magcof
  use coeffmodule

  implicit none

  integer(4), intent(in) :: ienty
  real(8), intent(in)    :: glat, glon, alt
  real(8), intent(out)   :: bnrth, beast, bdown, babs
  integer(4)             :: ientyp, ih, ihm, ihmax, il, ilm, imax
  integer(4)             :: is, k, last, m, mk, i
  real(8)                :: brho, bxxx, byyy, bzzz, cp, sp, ct, st, f
  real(8)                :: rlat, rlon, rq, s, t, x, xxx, y, yyy, z, zzz
  real(8)                :: g(255), h(255), xi(3)

  ientyp = - 10000

  if (ienty == 1) then
    is = 1
    rlat = glat * dtor
    ct = sin(rlat)
    st = cos(rlat)
    rlon = glon * dtor
    cp = cos(rlon)
    sp = sin(rlon)
    call gd2cart(glat, glon, alt, xxx, yyy, zzz)
    xxx = xxx / Re
    yyy = yyy / Re
    zzz = zzz / Re
  else
    is = 2
    xxx = glat
    yyy = glon
    zzz = alt
  end if
  rq = 1./ (xxx ** 2 + yyy ** 2 + zzz ** 2)
  xi(1) = xxx * rq
  xi(2) = yyy * rq
  xi(3) = zzz * rq
  ihmax = nmax1 * nmax1 + 1
  last = ihmax + nmax1 + nmax1
  imax = nmax1 + nmax1 - 1

  ! ientyp is not saved between calls, so this will always be run
  ! This may be incorrect?
  if (ienty /= ientyp .or. ichg == 1) then
    ientyp = ienty
    ichg = 0
    if (ienty /= 3) then
      do i = 1, last
        g(i) = gb(i)
      end do
    else
      do i = 1, last
        g(i) = gv(i)
      end do
    end if
  end if

  do i = ihmax, last
    h(i) = g(i)
  end do

  mk = 3
  if (imax == 1) mk = 1

  do k = 1, mk, 2
    i = imax
    ih = ihmax
    do while (i .ge. k)
      il = ih - i
      f = 2./ FLOAT(i - k + 2)
      x = xi(1) * f
      y = xi(2) * f
      z = xi(3) * (f + f)

      i = i - 2
      if (i > 1) then
        do m = 3, i, 2
          ihm = ih + m
          ilm = il + m
          h(ilm + 1) = g(ilm + 1) + z * h(ihm + 1) + x * (h(ihm + 3) - h(ihm - 1)) - y * (h(ihm + 2) + h(ihm - 2))
          h(ilm)   = g(ilm)  + z * h(ihm)   + x * (h(ihm + 2) - h(ihm - 2)) + y * (h(ihm + 3) + h(ihm - 1))
        end do
        h(il + 2) = g(il + 2) + z * h(ih + 2) + x * h(ih + 4) - y * (h(ih + 3) + h(ih))
        h(il + 1) = g(il + 1) + z * h(ih + 1) + y * h(ih + 4) + x * (h(ih + 3) - h(ih))
      elseif (i == 1) then
        h(il + 2) = g(il + 2) + z * h(ih + 2) + x * h(ih + 4) - y * (h(ih + 3) + h(ih))
        h(il + 1) = g(il + 1) + z * h(ih + 1) + y * h(ih + 4) + x * (h(ih + 3) - h(ih))
      end if
      h(il) = g(il) + z * h(ih) + 2.* (x * h(ih + 1) + y * h(ih + 2))
      ih = il
    end do
  end do

  s = .5 * h(1) + 2.* (h(2) * xi(3) + h(3) * xi(1) + h(4) * xi(2))
  t = (rq + rq) * sqrt(rq)
  bxxx = t * (h(3) - s * xxx)
  byyy = t * (h(4) - s * yyy)
  bzzz = t * (h(2) - s * zzz)
  babs = sqrt(bxxx ** 2 + byyy ** 2 + bzzz ** 2)
  if (is == 1) then            ! (convert back to geodetic)
    beast = byyy * cp - bxxx * sp
    brho  = byyy * sp + bxxx * cp
    bnrth = bzzz * st - brho * ct
    bdown = - bzzz * ct - brho * st
  elseif (is == 2) then        ! (leave in earth centered cartesian)
    bnrth = bxxx
    beast = byyy
    bdown = bzzz
  end if

  ! Magnetic potential computation makes use of the fact that the
  ! calculation of V is identical to that for r*Br, if coefficients
  ! in the latter calculation have been divided by (n+1) (coefficients
  ! GV).  Factor .1 converts km to m and gauss to tesla.
  if (ienty == 3) babs = (bxxx * xxx + byyy * yyy + bzzz * zzz) * Re *.1

  return
end subroutine feldg


subroutine gd2cart(gdlat, glon, alt, x, y, z)

  ! Convert geodetic to cartesian coordinates by calling CONVRT
  ! 940503 A. D. Richmond

  use magcof
  use coeffmodule

  implicit none

  real(8), intent(in)     :: gdlat, glon, alt
  real(8), intent(out)    :: x, y, z
  real(8)                 :: ang, rho

  call convrt(1, gdlat, alt, rho, z)
  ang = glon * dtor
  x = rho * cos(ang)
  y = rho * sin(ang)

  return
end subroutine gd2cart


subroutine convrt(i, gdlat, alt, x1, x2)

  ! Convert space point from geodetic to geocentric or vice versa.
  !
  ! i is an input flag controlling the meaning and direction of the
  ! remaining formal arguments:
  !
  ! i = 1  (convert from geodetic to cylindrical geocentric)
  ! INPUTS:
  ! gdlat = Geodetic latitude (deg)
  ! alt   = Altitude above reference ellipsoid (km)
  ! RETURNS:
  ! x1    = Distance from Earth's rotation axis (km)
  ! x2    = Distance above (north of) Earth's equatorial plane (km)
  !
  ! i = 2  (convert from geodetic to spherical geocentric)
  ! INPUTS:
  ! gdlat = Geodetic latitude (deg)
  ! alt   = Altitude above reference ellipsoid (km)
  ! RETURNS:
  ! x1    = Geocentric latitude (deg)
  ! x2    = Geocentric distance (km)
  !
  ! i = 3  (convert from cylindrical geocentric to geodetic)
  ! INPUTS:
  ! x1    = Distance from Earth's rotation axis (km)
  ! x2    = Distance from Earth's equatorial plane (km)
  ! RETURNS:
  ! gdlat = Geodetic latitude (deg)
  ! alt   = Altitude above reference ellipsoid (km)
  !
  ! i = 4  (convert from spherical geocentric to geodetic)
  ! INPUTS:
  ! x1    = Geocentric latitude (deg)
  ! x2    = Geocentric distance (km)
  ! RETURNS:
  ! gdlat = Geodetic latitude (deg)
  ! alt   = Altitude above reference ellipsoid (km)
  !
  !
  ! HISTORY:
  ! 940503 (A. D. Richmond):  Based on a routine originally written
  ! by V. B. Wickwar.
  !
  ! Mar 2004: (Barnes) Revise spheroid definition to WGS-1984 to conform
  ! with IGRF-9 release (EOS Volume 84 Number 46 November 18 2003).
  !
  ! Jul 2022 (Lamarche): Revise to fortran 90 standards
  !
  ! REFERENCE: ASTRON. J. VOL. 66, p. 15-16, 1961


  use magcof
  use coeffmodule

  implicit none

  integer(4), intent(in) :: i
  real(8)                :: gdlat, alt, x1, x2
  real(8)                :: a2, a4, a6, a8, c2cl, c4cl, ccl, coslat, sinlat, d, dltcl
  real(8)                :: gclat, rho, ri, rkm, s2cl, s4cl, s6cl, s8cl, scl, sgl, z

  ! A lot of these are the related to Req and esp - remove reduncancy
  real(8), parameter :: e2 = ecc2, e4 = e2 * e2, e6 = e4 * e2, e8 = e4 * e4,           &
                     ome2req = (1.- e2) * Req,                                     &
                     a21 = (512.* e2 + 128.* e4 + 60.* e6 + 35.* e8) / 1024.,         &
                     a22 = (e6 + e8) / 32., a23 = (e6 + e8) / 32.,                    &
                     a41 = - (64.* e4 + 48.* e6 + 35.* e8) / 1024.,                   &
                     a42 = (4.* e4 + 2.* e6 + e8) / 16., a43 = 15.* e8 / 256.,         &
                     a44 = - e8 / 16., a61 = 3.* (4.* e6 + 5.* e8) / 1024.,             &
                     a62 = - 3.* (e6 + e8) / 32., a63 = 35.* (4.* e6 + 3.* e8) / 768.,   &
                     a81 = - 5.* e8 / 2048., a82 = 64.* e8 / 2048.,                    &
                     a83 = - 252.* e8 / 2048., a84 = 320.* e8 / 2048.

  if (i < 3) then
    ! Geodetic to Geocentric
    ! Compute rho, z
    sinlat = sin(gdlat * dtor)
    coslat = sqrt(1.- sinlat * sinlat)
    d = sqrt(1.- e2 * sinlat * sinlat)
    z = (alt + ome2req / d) * sinlat
    rho = (alt + req / d) * coslat
    x1 = rho
    x2 = z
    if (i == 1) return

    ! Compute gclat, rkm
    rkm = sqrt(z * z + rho * rho)
    gclat = rtod * atan2(z, rho)
    x1 = gclat
    x2 = rkm
    return
  end if

  if (i == 3) then
    ! Geocentric to geodetic
    rho = x1
    z = x2
    rkm = sqrt(z * z + rho * rho)
    scl = z / rkm
    gclat = asin(scl) * rtod
  elseif (i == 4) then
    gclat = x1
    rkm = x2
    scl = sin(gclat * dtor)
  else
    return
  end if

  ri = req / rkm
  a2 = ri * (a21 + ri * (a22 + ri * a23))
  a4 = ri * (a41 + ri * (a42 + ri * (a43 + ri * a44)))
  a6 = ri * (a61 + ri * (a62 + ri * a63))
  a8 = ri * (a81 + ri * (a82 + ri * (a83 + ri * a84)))
  ccl = sqrt(1.- scl * scl)
  s2cl = 2.* scl * ccl
  c2cl = 2.* ccl * ccl - 1.
  s4cl = 2.* s2cl * c2cl
  c4cl = 2.* c2cl * c2cl - 1.
  s8cl = 2.* s4cl * c4cl
  s6cl = s2cl * c4cl + c2cl * s4cl
  dltcl = s2cl * a2 + s4cl * a4 + s6cl * a6 + s8cl * a8
  gdlat = dltcl * rtod + gclat
  sgl = sin(gdlat * dtor)
  alt = rkm * cos(dltcl) - req * sqrt(1.- e2 * sgl * sgl)
  return

end subroutine convrt
