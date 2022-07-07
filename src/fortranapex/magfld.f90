!**********************************************
!
!  File Name: magfld.f90
!  Authors:
!  Translated to Fortran 90 by Leslie Lamarche
!  Date: 01/27/2022
!
!***********************************************

module magfldmodule

  implicit none

  integer(4)       :: nmax1
  real(8)          :: ichg
  real(8)          :: gb(225), gv(225)

end module magfldmodule

module coeffmodule
  real(8), parameter       :: pi=3.14159265358979323846D0
  real(8), parameter       :: dtor=pi/180D0, rtod=180D0/pi, pid2=pi/2D0, twopi=2D0*pi
  real(8), parameter       :: Req=6378.1370D0, eps=1.D0/298.257223563D0
  real(8), parameter       :: Re=Req*(1-eps/3D0), ecc2=eps*(2-eps)
  ! real(4), parameter       :: missing=-9999E0
end module coeffmodule

! module igrfparammodule
!
!   real(4)             :: datel
!   integer(4)                  :: nepo, nght
!   real(4), allocatable     :: epoch(:), nmxe(:)
!   real(8), allocatable                     :: gyr(:,:,:), hyr(:,:,:)
!   real(8), allocatable                     :: gt(:,:), ht(:,:)
!
! end module igrfparammodule

subroutine cofrm(date,filename)

    use magfldmodule
    use igrf
    ! use igrfparammodule

    implicit none

    ! COMMON /MAGCOF/ nmax1, GB, GV, ICHG

    real(4)                     :: date
    real(4)                     :: f, f0, t, to5
    integer   :: i, i1, iy, iy1, m, n, mm, nn, ngh

    character(len=1000)       :: filename

    ! integer(4)                  :: nmax1
    ! real(4)                     :: GB(1:255), GV(1:225)
    ! real(4)                     :: ICHG

    ! igrf parameters are saved so they only need to be loaded once
    ! SAVE DATEL, GYR, HYR, GT, HT, NEPO, EPOCH, NGHT, NMXE
    ! Commenting these out might break something...
    ! DATA datel /-999./
    ! DATA ichg /-99999/
    ichg = 0
    if (date .eq. datel) then
      return
    else
      datel = date
      ichg = 1
    endif

    if (.not. allocated(gyr)) then
      ! call read_igrf(filename, gyr, hyr, gt, ht, nepo, nght, epoch, nmxe)
      call read_igrf(filename)
    endif
    ngh = nght*nepo

    if (date .lt. epoch(1)) then
      write(0, '("COFRM:  DATE "(F9.3)" preceeds earliest available "(F6.1))') date, epoch(1)
      call exit(1)
    elseif (date .gt. epoch(nepo)+5.) then
      write(0, '("COFRM:  DATE "(F9.3)" is after the last recommended for extrapolation "(F6.1))') date, epoch(1)+5.
      call exit(1)
    endif

    iy = 1
    do while (date .gt. epoch(iy))
      iy = iy + 1
    enddo

    ngh = nght*nepo
    nmax1 = nmxe(iy)
    ! time = date
    t = date-epoch(iy)
    to5 = t/5.
    iy1 = iy + 1
    gb(1) = 0.0
    gv(1) = 0.0
    i = 2
    f0 = -1.0d-5

    do n=1,nmax1
      f0 = f0*real(n)/2.
      f = f0/sqrt(2.0)
      nn = n+1
      mm =1
      if (iy .lt. nepo) then
        gb(i) = (gyr(nn,mm,iy) + (gyr(nn,mm,iy1) - gyr(nn,mm,iy))*to5) * f0
      endif
      if (iy .eq. nepo) then
        gb(i) = (gyr(nn,mm,iy) + gt(nn,mm)*t) * f0
      endif
      gv(i) = gb(i)/real(nn)
      i = i+1
      do m=1,n
        f = f/sqrt(real(n-m+1)/real(n+m))
        nn = n+1
        mm = m+1
        i1 = i+1
        if (iy .lt. nepo) then
          gb(i) = (gyr(nn,mm,iy) + (gyr(nn,mm,iy1) - gyr(nn,mm,iy))*to5) * f
          gb(i1) = (hyr(nn,mm,iy) + (hyr(nn,mm,iy1) - hyr(nn,mm,iy))*to5) * f
        else
          gb(i) = (gyr(nn,mm,iy) + gt(nn,mm)*t) * f
          gb(i1) = (hyr(nn,mm,iy) + ht(nn,mm)*t) * f
        endif
        gv(i) = gb(i)/real(nn)
        gv(i1) = gb(i1)/real(nn)
        i = i+2
      enddo
    enddo

    return
end subroutine cofrm


subroutine dypol(colat, elon, vp)

  use magfldmodule
  use coeffmodule

  implicit none

  real(8)       :: colat, elon, vp
  real(8)       :: ctp, gpl, stp

  ! real, parameter :: rtod = 57.2957795130823
  ! real, parameter :: re = 6371.2
  ! COMMON /MAGCOF/ nmax1, GB(255), GV(225), ICHG

  gpl = sqrt(gb(2)**2 + gb(3)**2 + gb(4)**2)
  ctp = gb(2)/gpl
  stp = sqrt(1 - ctp*ctp)
  colat = acos(ctp)*rtod
  elon = atan2(gb(4),gb(3))*rtod

  ! Compute magnitude of magnetic potential at pole, radius Re
  vp = 0.2*gpl*Re

  return
end subroutine dypol


subroutine feldg(ienty, glat, glon, alt, bnrth, beast, bdown, babs)

  use magfldmodule
  use coeffmodule

  implicit none

  integer(4)         :: ienty, ientyp, ih, ihm, ihmax, il, ilm, imax
  integer(4)         :: is, k, last, m, mk, i
  real(8), intent(in)    :: glat, glon, alt
  real(8), intent(out)   :: bnrth, beast, bdown, babs
  real(8)            :: brho, bxxx, byyy, bzzz, cp, sp, ct, st, f
  real(8)            :: rlat, rlon, rq, s, t, x, xxx, y, yyy, z, zzz

  ! real, parameter :: rtod = 57.2957795130823
  ! real, parameter :: re = 6371.2
  ! COMMON /MAGCOF/ nmax1,GB(255),GV(225),ICHG
  ! DIMENSION G(255), H(255), XI(3)

  real(8)            :: g(255), h(255), xi(3)
  ! Not sure what's going on with these two
  ! SAVE IENTYP, g
  ! DATA IENTYP/-10000/
  ! print *, 'LINE 129'
  if (ienty .eq. 1) then
    is = 1
    rlat = glat*dtor
    ct = sin(rlat)
    st = cos(rlat)
    rlon = glon*dtor
    cp = cos(rlon)
    sp = sin(rlon)
    call gd2cart(glat, glon, alt, xxx, yyy, zzz)
    xxx = xxx/Re
    yyy = yyy/Re
    zzz = zzz/Re
  else
    is = 2
    xxx = glat
    yyy = glon
    zzz = alt
  endif
  rq = 1./(xxx**2 + yyy**2 + zzz**2)
  xi(1) = xxx*rq
  xi(2) = yyy*rq
  xi(3) = zzz*rq
  ihmax = nmax1*nmax1+1
  last = ihmax+nmax1+nmax1
  imax = nmax1+nmax1-1
  ! print *, 'LINE 155'
  if (ienty .ne. ientyp .or. ichg .eq. 1) then
    ientyp = ienty
    ichg = 0
    if (ienty .ne. 3) then
      do i=1,last
        g(i) = gb(i)
      enddo
    else
      do i=1,last
        g(i) = gv(i)
      enddo
    endif
  endif
  ! print *, 'LINE 169'
  do i=ihmax,last
    h(i) = g(i)
  enddo
  ! print *, 'LINE 173'
  mk = 3
  if (imax .eq. 1) mk=1
  ! print *, 'LINE 176'
  do k=1,mk,2
    i = imax
    ih = ihmax
    ! print *, 'i=', i, 'k=', k
    do while (i .ge. k)
      il = ih-i
      f = 2./FLOAT(i-k+2)
      x = xi(1)*f
      y = xi(2)*f
      z = xi(3)*(f+f)
      ! print *, 'LINE 187'
      i = i-2
      ! IF (I .LT. 1) GO TO 90
      ! IF (I .EQ. 1) GO TO 80
      ! print *, 'LINE 191'
      if (i .gt. 1) then
        do m=3,i,2
          ihm = ih+m
          ilm = il+m
          ! print *, 'm=', m, 'ihm=', ihm, 'ilm=', ilm
          ! print *, 'LINE 197'
          h(ilm+1) = g(ilm+1)+ z*h(ihm+1) + x*(h(ihm+3)-h(ihm-1)) - y*(h(ihm+2)+h(ihm-2))
          ! print *, 'LINE 199'
          h(ilm)   = g(ilm)  + z*h(ihm)   + x*(h(ihm+2)-h(ihm-2)) + y*(h(ihm+3)+h(ihm-1))
          ! print *, 'LINE 201'
        enddo
        h(il+2) = g(il+2) + z*h(ih+2) + x*h(ih+4) - y*(h(ih+3)+h(ih))
        h(il+1) = g(il+1) + z*h(ih+1) + y*h(ih+4) + x*(h(ih+3)-h(ih))

      elseif (i .eq. 1) then
        h(il+2) = g(il+2) + z*h(ih+2) + x*h(ih+4) - y*(h(ih+3)+h(ih))
        h(il+1) = g(il+1) + z*h(ih+1) + y*h(ih+4) + x*(h(ih+3)-h(ih))

      endif
      ! print *, 'LINE 207'
      h(il) = g(il) + z*h(ih) + 2.*(x*h(ih+1)+y*h(ih+2))
      ih = il
    enddo
    ! IF (I .GE. K) GO TO 60
  enddo
  ! print *, 'LINE 213'
  s = .5*h(1)+2.*(h(2)*xi(3)+h(3)*xi(1)+h(4)*xi(2))
  t = (rq+rq)*sqrt(rq)
  bxxx = t*(h(3)-s*xxx)
  byyy = t*(h(4)-s*yyy)
  bzzz = t*(h(2)-s*zzz)
  babs = sqrt(bxxx**2+byyy**2+bzzz**2)
  if (is .eq. 1) then            ! (convert back to geodetic)
    beast = byyy*cp-bxxx*sp
    brho  = byyy*sp+bxxx*cp
    bnrth = bzzz*st-brho*ct
    bdown = -bzzz*ct-brho*st
  elseif (is .eq. 2) then        ! (leave in earth centered cartesian)
    bnrth = bxxx
    beast = byyy
    bdown = bzzz
  endif
  ! print *, 'LINE 230'
  if (ienty .eq. 3) babs = (bxxx*xxx + byyy*yyy + bzzz*zzz)*Re*.1
  ! print *, 'LINE 232'
  return
end subroutine feldg


subroutine gd2cart(gdlat, glon, alt, x, y, z)

  use magfldmodule
  use coeffmodule

  implicit none

  real(8)           :: gdlat, glon, alt, x, y, z
  real(8)           :: ang, rho
  ! real, parameter :: dtor = 0.01745329251994330

  call convrt(1, gdlat, alt, rho, z)
  ang = glon*dtor
  x = rho*cos(ang)
  y = rho*sin(ang)

  return
end subroutine gd2cart


subroutine convrt(i, gdlat, alt, x1, x2)

  use magfldmodule
  use coeffmodule

  implicit none

  integer(4)        :: i
  real(8)           :: gdlat, alt, x1, x2
  real(8)           :: a2, a4, a6, a8, c2cl, c4cl, ccl, coslat, sinlat, d, dltcl
  real(8)           :: gclat, rho, ri, rkm, s2cl, s4cl, s6cl, s8cl, scl, sgl, z

  ! A lot of these are the related to Req and esp - remove reduncancy
  real, parameter :: e2 = ecc2, e4 = e2*e2, e6 = e4*e2, e8 = e4*e4,           &
                     ome2req = (1.-e2)*Req,                                     &
                     a21 = (512.*e2 + 128.*e4 + 60.*e6 + 35.*e8)/1024.,         &
                     a22 = (e6 + e8)/32., a23 = (e6+e8)/32.,                    &
                     a41 = -(64.*e4 + 48.*e6 + 35.*e8)/1024.,                   &
                     a42 = (4.*e4 + 2.*e6 + e8)/16., a43 = 15.*e8/256.,         &
                     a44 = -e8/16., a61 = 3.*(4.*e6 + 5.*e8)/1024.,             &
                     a62 = -3.*(e6 + e8)/32., a63 = 35.*(4.*e6 + 3.*e8)/768.,   &
                     a81 = -5.*e8/2048., a82 = 64.*e8/2048.,                    &
                     a83 = -252.*e8/2048., a84 = 320.*e8 /2048.

  if (i .lt. 3) then
    sinlat = sin(gdlat*dtor)
    coslat = sqrt(1.-sinlat*sinlat)
    d = sqrt(1.-e2*sinlat*sinlat)
    z = (alt+ome2req/d)*sinlat
    rho = (alt+req/d)*coslat
    x1 = rho
    x2 = z
    if (i .eq. 1) return

    rkm = sqrt(z*z + rho*rho)
    gclat = rtod*atan2(z,rho)
    x1 = gclat
    x2 = rkm
    return
  endif

  if (i .eq. 3) then
    rho = x1
    z = x2
    rkm = sqrt(z*z+rho*rho)
    scl = z/rkm
    gclat = asin(scl)*rtod
  elseif (i .eq. 4) then
    gclat = x1
    rkm = x2
    scl = sin(gclat*dtor)
  else
    return
  endif

  ri = req/rkm
  a2 = ri*(a21+ri*(a22+ri*a23))
  a4 = ri*(a41+ri*(a42+ri*(a43+ri*a44)))
  a6 = ri*(a61+ri*(a62+ri*a63))
  a8 = ri*(a81+ri*(a82+ri*(a83+ri*a84)))
  ccl = sqrt(1.-scl*scl)
  s2cl = 2.*scl*ccl
  c2cl = 2.*ccl*ccl-1.
  s4cl = 2.*s2cl*c2cl
  c4cl = 2.*c2cl*c2cl-1.
  s8cl = 2.*s4cl*c4cl
  s6cl = s2cl*c4cl+c2cl*s4cl
  dltcl = s2cl*a2+s4cl*a4+s6cl*a6+s8cl*a8
  gdlat = dltcl*rtod+gclat
  sgl = sin(gdlat*dtor)
  alt = rkm*cos(dltcl)-req*sqrt(1.-e2*sgl*sgl)
  return

end subroutine convrt
