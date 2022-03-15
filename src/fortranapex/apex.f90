! apex.f90

module apexmodule

  implicit none

! ! check all these to make sure they actually need to be shared variables
!   real(8)                  :: colat, elon, vp, ctp, stp
!   real(8)                  :: bx, by, bz, bb
!   real(8)                  :: yapx(3,3)
!   integer(4)               :: nstp
!   real(8)                  :: sgn, ds
!   real(8)                  :: y(3)

  real(8), parameter       :: pi=3.14159265358979323846D0
  real(8), parameter       :: dtor=pi/180D0, rtod=180D0/pi, pid2=pi/2D0, twopi=2D0*pi
  real(8), parameter       :: Req=6378.1370D0, eps=1.D0/298.257223563D0
  real(8), parameter       :: Re=Req*(1-eps/3D0), ecc2=eps*(2-eps)
  real(4), parameter       :: missing=-9999E0

contains

  real(8) function fint(x1, x2, x3, y1, y2, y3, xfit)

    real(8), intent(in)  :: x1, x2, x3, y1, y2, y3, xfit
    real(8)              :: x12, x13, x23, xf1, xf2, xf3

    x12 = x1-x2
    x13 = x1-x3
    x23 = x2-x3
    xf1 = xfit-x1
    xf2 = xfit-x2
    xf3 = xfit-x3

    fint = (y1*x23*xf2*xf3 - y2*x13*xf1*xf3 + y3*x12*xf1*xf2)/(x12*x13*x23)
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
  real(8)                  :: yapx(3,3)
end module apxin

module itra
  implicit none
  integer(4)               :: nstp
  real(8)                  :: sgn, ds
  real(8)                  :: y(3), yold(3), yp(3,4)
end module itra
! COMMON /FLDCOMD/ BX, BY, BZ, BB
! COMMON /APXIN/   YAPX(3,3)
! COMMON /DIPOLE/  COLAT,ELON,VP,CTP,STP
! COMMON /ITRA/    NSTP, Y(3), YP(3), SGN, DS

subroutine apex(date, igrffilein, dlat, dlon, alt, a, alat, alon, bmag, xmag, ymag, zmag, v)

  use dipole
  use apexmodule

  implicit none

  ! real, parameter :: re = 6371.0088, dtor = .01745329251994330
  ! COMMON /DIPOLE/ COLAT,ELON,VP,CTP,STP

  real(4), intent(in)      :: date
  character(len=1000), intent(in)   :: igrffilein
  real(8), intent(in)      :: dlat, dlon, alt
  real(8)      :: a, alon, bmag, xmag, ymag, zmag, v
  real(8)      :: alat, bx, by, bz, x, y, z
  real(8)      :: clatp, polon, vpol

  call cofrm(date, igrffilein)
  call dypol(clatp, polon, vpol)
  colat = clatp
  ctp = cos(clatp*dtor)
  stp = sqrt(1 - CTP*CTP)
  elon = polon
  vp = vpol
  ! print *, 'LINE 17', alt
  call linapx(dlat, dlon, alt, a, alat, alon, xmag, ymag, zmag, bmag)
  ! print *, 'LINE 19'
  xmag = xmag*1.e5
  ymag = ymag*1.e5
  zmag = zmag*1.e5
  bmag = bmag*1.e5
  ! print *, 'LINE 24'
  call gd2cart(dlat, dlon, alt, x, y, z)
  ! print *, 'LINE 26'
  call feldg(3, x/Re, y/Re, z/Re, bx, by, bz, v)
  ! print *, 'LINE 25'
  return
end subroutine apex


subroutine linapx(gdlat, glon, alt, a, alat, alon, xmag, ymag, zmag, f)

  use fldcomd
  use apxin
  use dipole
  use itra
  use apexmodule

  implicit none

  real(8)               :: gdlat, glon, alt, a, alat, alon, xmag, ymag, zmag, f
  real(8)               :: bnrth, beast, bdown, babs
  real(8)               :: cgml2, gclat, gcml2, ht, r, rho
  real(8)               :: singml, xlat, xlon
  integer(4)            :: i, j, iapx
  integer(4), parameter :: maxs = 200
  ! rtod = 57.2957795130823, re = 6371.0088, &
                     ! dtor = .01745329251994330, req = 6378.137
  ! COMMON /FLDCOMD/ BX, BY, BZ, BB
  ! COMMON /APXIN/   YAPX(3,3)
  ! COMMON /DIPOLE/  COLAT,ELON,VP,CTP,STP
  ! COMMON /ITRA/    NSTP, Y(3), YP(3), SGN, DS
  ! print *, 'LINE 41'
  call convrt(2, gdlat, alt, gclat, r)
  singml = ctp*sin(gclat*dtor) + stp*cos(gclat*dtor)*cos((glon-elon)*dtor)
  cgml2 = amax1(0.25, 1.-singml*singml)
  ds = 0.06*r/cgml2 - 370.
  ! print *, 'LINE 46'
  do j=1,3
    do i=1,3
      yapx(i,j) = 0.
    enddo
  enddo
  ! print *, 'LINE 52'
  call gd2cart(gdlat, glon, alt, y(1), y(2), y(3))
  ! print *, gdlat, glon, alt, y
  nstp = 0
  ! print *, 'LINE 55'
  call feldg(1, gdlat, glon, alt, xmag, ymag, zmag, f)
  ! print *, xmag, ymag, zmag, f
  ! print *, 'LINE 57'
  sgn = sign(1.D0, -zmag)
  ! print *, 'LINE 58'
  ! call feldg(2, y(1)/Re, y(2)/Re, y(3)/Re, bx, by, bz, bb)
  ! print *, bx, by, bz, bb
  ! nstp = nstp+1
  ! print *, 'LINE 61'
  iapx = 1
  do while (iapx .EQ. 1)
    call feldg(2, y(1)/Re, y(2)/Re, y(3)/Re, bx, by, bz, bb)
    ! print *, 'LINE 158', y
    nstp = nstp+1
    if (nstp .ge. maxs) exit
    call itrace(iapx)
    ! print *, 'LINE 162', y
  enddo
  if (nstp .lt. maxs) then
    call fndapx(alt, zmag, a, alat, alon)
  else
    rho = sqrt(y(1)**2 + y(2)**2)
    call convrt(3, xlat, ht, rho, y(3))
    xlon = rtod*atan2(y(2), y(1))
    call feldg(1, xlat, xlon, ht, bnrth, beast, bdown, babs)
    call dipapx(xlat, xlon, ht, bnrth, beast, bdown, a, alon)
    alat = -sgn*rtod*acos(sqrt(1./a))
  endif
  ! print *, 'LINE 79'
  return
end subroutine linapx


subroutine itrace(iapx)

  use apxin
  use fldcomd
  use itra
  ! use apexmodule

  implicit none

  integer(4), intent(out)  :: iapx
  real(8)                  :: d12, d2, d24, d6, rc, rp
  integer(4)               :: i, j
  ! real(8)                  :: yold(3), yp(3,4)

  ! COMMON /APXIN/   YAPX(3,3)
  ! COMMON /FLDCOMD/ BX, BY, BZ, BB
  ! COMMON /ITRA/    NSTP, Y(3), YOLD(3), SGN, DS
  ! DIMENSION YP(3,4)
  ! SAVE

  ! rdus(d,e,f) = sqrt(d**2 + e**2 + f**2)

  iapx = 1

  yp(1,4) = sgn*bx/bb
  yp(2,4) = sgn*by/bb
  yp(3,4) = sgn*bz/Bb

  d2 = ds/2.
  d6 = ds/6.
  d12 = ds/12.
  d24 = ds/24.

  if (nstp .le. 7) then
    do i=1,3
      if (nstp .eq. 1) then
        ! d2 = ds/2.
        ! d6 = ds/6.
        ! d12 = ds/12.
        ! d24 = ds/24.
        yp(i,1) = yp(i,4)
        yold(i) = y(i)
        yapx(i,1) = y(i)
        y(i) = yold(i) + ds*yp(i,1)

      else if (nstp .eq. 2) then
        yp(i,2) = yp(i,4)
        y(i) = yold(i) + d2*(yp(i,2)+yp(i,1))

      else if (nstp .eq. 3) then
        y(i) = yold(i) + d6*(2.*yp(i,4)+yp(i,2)+3.*yp(i,1))

      else if (nstp .eq. 4) then
        yp(i,2) = yp(i,4)
        yapx(i,2) = y(i)
        yold(i) = y(i)
        y(i) = yold(i) + d2*(3.*yp(i,2)-yp(i,1))

      else if (nstp .eq. 5) then
        y(i) = yold(i) + d12*(5.*yp(i,4)+8.*yp(i,2)-yp(i,1))

      else if (nstp .eq. 6) then
        yp(i,3) = yp(i,4)
        yold(i) = y(i)
        yapx(i,3) = y(i)
        y(i) = yold(i) + d12*(23.*yp(i,3)-16.*yp(i,2)+5.*yp(i,1))

      else if (nstp .eq. 7) then
        yapx(i,1) = yapx(i,2)
        yapx(i,2) = yapx(i,3)
        y(i) = yold(i) + d24*(9.*yp(i,4)+19.*yp(i,3)-5.*yp(i,2)+yp(i,1))
        yapx(i,3) = y(i)

      endif
    enddo

    if (nstp .eq. 6 .or. nstp .eq. 7) then
      ! rc = rdus(YAPX(1,3), YAPX(2,3), YAPX(3,3))
      ! rp = rdus(YAPX(1,2), YAPX(2,2), YAPX(3,2))
      rc = sqrt(yapx(1,3)**2 + yapx(2,3)**2 + yapx(3,3)**2)
      rp = sqrt(yapx(1,2)**2 + yapx(2,2)**2 + yapx(3,2)**2)
      ! print *, yapx
      if (rc .lt. rp) iapx=2
    endif

  else
    do i=1,3
      yapx(i,1) = yapx(i,2)
      yapx(i,2) = y(i)
      yold(i) = y(i)
      y(i) = yold(i) + d24*(55.*yp(i,4)-59.*yp(i,3)+37.*yp(i,2)-9.*yp(i,1))
      yapx(i,3) = y(i)
      do j=1,3
        yp(i,j) = yp(i,j+1)
      enddo
    enddo
    ! rc = rdus(Y(1), Y(2), Y(3))
    ! rp = rdus(YOLD(1), YOLD(2), YOLD(3))
    rc = sqrt(y(1)**2 + y(2)**2 + y(3)**2)
    rp = sqrt(yold(1)**2 + yold(2)**2 +yold(3)**2)
    if (rc .lt. rp) iapx=2

  endif

  return
end subroutine itrace


subroutine fndapx(alt, zmag, a, alat, alon)

  use apxin
  use dipole
  use apexmodule

  implicit none

  ! real, parameter :: rtod = 57.2957795130823, dtor = .01745329251994330, req = 6378.137
  ! COMMON /APXIN/  YAPX(3,3)
  ! COMMON /DIPOLE/ COLAT,ELON,VP,CTP,STP
  ! DIMENSION BD(3), Y(3)
  real(8)       :: alt, zmag, a, alat, alon
  real(8)       :: abdob, ang, ba, bda, bea, bna, hta, rho
  real(8)       :: r, rasq, sang, cang, ste, cte, stfcpa, stfspa, xlon
  real(8)       :: gdlt, gdln, ht, bn, be, bmag
  integer(4)    :: i, nitr
  real(8)       :: bd(3), y(3)

  do i=1,3
    rho = sqrt(yapx(1,i)**2 + yapx(2,i)**2)
    call convrt(3, gdlt, ht, rho, yapx(3,i))
    gdln = rtod*atan(yapx(2,i),yapx(1,i))
    call feldg(1, gdlt, gdln, ht, bn, be, bd(i), bmag)
  enddo
  ! print *, bd, yapx

  nitr = 0
  do while (nitr .lt. 4)
    y(1) = fint(bd(1), bd(2), bd(3), yapx(1,1), yapx(1,2), yapx(1,3), 0.D0)
    y(2) = fint(bd(1), bd(2), bd(3), yapx(2,1), yapx(2,2), yapx(2,3), 0.D0)
    y(3) = fint(bd(1), bd(2), bd(3), yapx(3,1), yapx(3,2), yapx(3,3), 0.D0)

    rho = sqrt(y(1)**2 + y(2)**2)
    gdln = rtod*atan2(y(2), y(1))
    call convrt(3, gdlt, hta, rho, y(3))
    call feldg(1, gdlt, gdln, hta, bna, bea, bda, ba)
    abdob = abs(bda/ba)

    if (abdob .gt. 2.e-6) then
      nitr = nitr + 1
      yapx(1,2) = y(1)
      yapx(2,2) = y(2)
      yapx(3,2) = y(3)
      bd(2) = bda
    else
      exit
    endif
  enddo

  if (abdob .gt. 2.e-6) then
    write(0, '("APEX: Imprecise fit of apex: |Bdown/B| "(1PE7.1))') abdob
  endif

  a = (Req + amax1(alt, hta))/Req
  if (a .lt. 1.) then
    ! write(0, '("APEX: A cannot be less than 1; A, REQ, HTA "(1P3))') a, Req, hta
    write(0, *) "APEX: A cannot be less than 1; A, REQ, HTA ", a, Req, hta
    call exit(1)
  endif
  rasq = acos(sqrt(1./a))*rtod
  alat = sign(rasq, zmag)

  xlon = atan2(y(2), y(1))
  ang = xlon-elon*dtor
  cang = cos(ang)
  sang = sin(ang)
  r = sqrt(y(1)**2 + y(2)**2 + y(3)**2)
  cte = y(3)/r
  ste = sqrt(1.-cte*cte)
  stfcpa = ste*ctp*cang - cte*stp
  stfspa = sang*ste
  alon = atan2(stfspa,stfcpa)*rtod

  return
end subroutine fndapx


subroutine dipapx(gdlat, gdlon, alt, bnorth, beast, bdown, a, alon)

  use dipole
  use apexmodule

  implicit none

  real(8)        :: gdlat, gdlon, alt, bnorth, beast, bdown, a, alon
  real(8)        :: ang, bhor, ca, cang, capang, cb, cottd, ctd, cte, ctg
  real(8)        :: ha, r, sa, sang, sapabg, sb, std, ste, stfcpa, stg, sapang, stfspa

  ! real, parameter :: rtod = 57.2957795130823, re = 6371.0088, &
  !                    dtor = .01745329251994330, req = 6378.137
  ! COMMON /DIPOLE/ COLAT,ELON,VP,CTP,STP

  bhor = sqrt(bnorth*bnorth + beast*beast)
  if (bhor .eq. 0.) then
    alon = 0.
    a = 1.E34
    return
  endif

  cottd = bdown*0.5/bhor
  std = 1./sqrt(1.+cottd*cottd)
  ctd = cottd*std
  sb = -beast/bhor
  cb = -bnorth/bhor
  ctg = sin(gdlat*dtor)
  stg = cos(gdlat*dtor)
  ang = (gdlon - elon)*dtor
  sang = sin(ang)
  cang = cos(ang)
  cte = ctg*std + stg*ctd*cb
  ste = sqrt(1.-cte*cte)
  sa = sb*ctd/ste
  ca = (std*stg - ctd*ctg*cb)/ste
  capang = ca*cang - sa*sang
  sapang = ca*sang + sa*cang
  stfcpa = ste*ctp*capang - cte*stp
  stfspa = sapang*ste
  alon = atan2(stfspa,stfcpa)*rtod
  r = alt + Re
  ha = alt + r*cottd*cottd
  a = 1. + ha/Req

  return
end subroutine dipapx


! real(8) function fint(x1, x2, x3, y1, y2, y3, xfit)
!
!   implicit none
!
!   real(8), intent(in)  :: x1, x2, x3, y1, y2, y3, xfit
!   real(8)              :: x12, x13, x23, xf1, xf2, xf3
!
!   x12 = x1-x2
!   x13 = x1-x3
!   x23 = x2-x3
!   xf1 = xfit-x1
!   xf2 = xfit-x2
!   xf3 = xfit-x3
!
!   fint = (y1*x23*xf2*xf3 - y2*x13*xf1*xf3 + y3*x12*xf1*xf2)/(x12*x13*x23)
!   return
! end function fint
