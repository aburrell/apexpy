! ******************************************************************************
!
! File Name: makeapexsh.f90
! Authors: John Emmert, Art Richmond
! Date: 11/13/2009
! Version: 1.0
! Description: Creates and saves spherical harmonic expansion coefficients for
! QD coordinate conversion.
! References: Richmond, A. D., Ionospheric Electrodynamics Using Magnetic Apex
! Coordinates, J. Geomag. Geoelectr., 47, 191-212, 1995.
! Emmert, J. T., A. D. Richmond, and D. P. Drob, A computationally
! compact representation of Magnetic-Apex and Quasi-Dipole
! coordinates with smooth base vectors, J. Geophys. Res., 115,
! A08322, doi:10.1029/2010JA015326, 2010.
!
! ******************************************************************************
!
! HISTORY (blame):
!
! 25 Feb 2021: Modified by Ashton Reimer to pass IGRF coefficients file to the
! cofrm and apex subroutine calls.
!
! ******************************************************************************
!
! makeapxsh
! Computes and saves harmonic coefficients for coordinate conversions.
!
! call makeapxsh(datafilein, igrffilein, epochgridin, nepochin, lmaxin, mmaxin, nmaxin)
!
! INPUT:
! datafilein  = Name of output data file that will contain the conversion
! coefficients
! igrffilein  = Name of the input IGRF coefficients file
! epochgridin = Array of ordered epochs (decimal years, yyyy.y) for which
! coefficients are to be computed.
! nepochin    = Number of elements in EPOCHS.
! lmaxin      = Maximum order of vertical polynomial expansion.
! mmaxin      = Maximum order of spherical harmonic expansion.
! nmaxin      = Maximum degree of spherical harmonic expansion. N must be equal
! or greater than M.
!
! DEPENDENCIES
! apex.f90, magfld.f90, apexsh.f90
!
! ******************************************************************************

subroutine makeapxsh(datafilein, igrffilein, epochgridin, nepochin, lmaxin, mmaxin, nmaxin)

    use apxshmodule
    use magcof

    implicit none

    character(128), intent(in)      :: datafilein
    character(len=1000), intent(in) :: igrffilein
    real(8), intent(in)             :: epochgridin(0:31)
    integer(4), intent(in)          :: nmaxin, mmaxin, lmaxin, nepochin

    integer(4), parameter       :: sampfact=4, iun=12
    integer(4)                  :: nterm1, nlon, nlat, nalt
    integer(4)                  :: iepoch, iterm, ilon, ilat, ialt
    integer(4)                  :: l, ish, i, j, Rpt(0:2)
    real(8)                     :: glon, glat, alt, lonspace, latspace, lat0
    real(8)                     :: rhospace
    real(8)                     :: norm1, norm2, cosmplat, sinmplat, cosmplon
    real(8)                     :: sinmplon
    real(8)                     :: thetag, phig, thetaq, phiq, latwgt
    real(8)                     :: xq0, yq0, zq0, xg, yg, zg
    real(8), allocatable        :: altgrid(:), altfn(:,:), shg(:)
    real(8), allocatable        :: Gg(:,:), Gq(:,:), Fg(:), Fq(:), cfdiag(:)
    real(8), allocatable        :: coefftemp(:)
    real(8), allocatable        :: Dxq(:), Dyq(:), Dzq(:), Dxg(:), Dyg(:)
    real(8), allocatable        :: Dzg(:)
    real(8)                     :: a, alat, alon, dum1, dum2, dum3, dum4, dum5

    external cofrm, apex

    ! PARSE INPUT AND ALLOCATE ARRAYS
    nepoch = nepochin
    nmax = nmaxin
    mmax = mmaxin
    lmax = lmaxin
    if (mmax > nmax) mmax = nmax
    ntermsh = mmax * (2 * nmax - mmax + 1) + nmax + 1
    nterm = (lmax + 1) * ntermsh
    nterm1 = lmax * ntermsh
    Rpt = (/ nmax + 1, nmax + 2, 1 /)
    nlon = sampfact * mmax
    nlat = sampfact * nmax
    nalt = sampfact * lmax
    call allocatearrays
    coeff0 = 0
    do iepoch = 0, nepoch - 1
      epochgrid(iepoch) = epochgridin(iepoch)
    end do
    if (allocated(altgrid)) then
      deallocate(altgrid, Gg, Gq, Fg, Fq, Dxq, Dyq, Dzq, Dxg, Dyg, Dzg, shg, cfdiag, coefftemp)
    end if
    allocate(altgrid(1:nalt), altfn(0:lmax, 1:nalt), shg(0:ntermsh - 1))
    allocate(Gg(0:nterm1 - 1, 0:nterm1 - 1), Gq(0:nterm1 - 1, 0:nterm1 - 1))
    allocate(Fg(0:nterm1 - 1), Fq(0:nterm1 - 1), cfdiag(0:nterm1 - 1), coefftemp(0:nterm1 - 1))
    allocate(Dxq(0:nterm1 - 1), Dyq(0:nterm1 - 1), Dzq(0:nterm1 - 1))
    allocate(Dxg(0:nterm1 - 1), Dyg(0:nterm1 - 1), Dzg(0:nterm1 - 1))
    lonspace = 360E0 / real(nlon)
    latspace = 180E0 / real(nlat)
    lat0 = - 90E0 + latspace / 2D0
    rhospace = 1D0 / dble(nalt)
    norm1 = dsqrt(2D0 / 3D0)
    norm2 = dsqrt(4D0 / 3D0)

    ! COMPUTE ALTITUDE BASIS FUNCTIONS
    do ialt = 1, nalt
      rho = rhospace * dble(ialt)
      altgrid(ialt) = Re / rho - Re
      altfn(0, ialt) = 1D0
      do l = 1, lmax
        altfn(l, ialt) = altfn(l - 1, ialt) * rho
      end do
    end do

    ! LOOP THROUGH EPOCHS, COMPUTE QD COORDINATE DATA AND EXPANSION COEFFICIENTS
    do iepoch = 0, nepoch - 1
      print '("COMPUTING COEFFICIENTS FOR EPOCH ",f6.1)', epochgrid(iepoch)
      Gg = 0; Gq = 0
      Dxq = 0; Dyq = 0; Dzq = 0
      Dxg = 0; Dyg = 0; Dzg = 0

      ! Retrieve IGRF
      call cofrm(epochgrid(iepoch), igrffilein)

      ! Compute the dipole rotation parameters
      sinmplat = dble(gb(2) / sqrt(gb(2) * gb(2) + gb(3) * gb(3) + gb(4) * gb(4)))
      cosmplat = dsqrt(1 - sinmplat * sinmplat)
      sinmplon = dble(gb(4) / sqrt(gb(3) * gb(3) + gb(4) * gb(4)))
      cosmplon = dble(gb(3) / sqrt(gb(3) * gb(3) + gb(4) * gb(4)))

      ! Set the expansion coeficients
      coeff0(Rpt, iepoch, 0) = (/ norm2 * sinmplat * cosmplon, norm2 * sinmplat * sinmplon, - norm1 * cosmplat /)
      coeff0(Rpt, iepoch, 1) = (/- norm2 * sinmplon,          norm2 * cosmplon,                      0D0 /)
      coeff0(Rpt, iepoch, 2) = (/ norm2 * cosmplat * cosmplon, norm2 * cosmplat * sinmplon,  norm1 * sinmplat /)
      coeff0(Rpt, iepoch, 3) = (/ norm2 * sinmplat * cosmplon, - norm2 * sinmplon, norm1 * cosmplat * cosmplon /)
      coeff0(Rpt, iepoch, 4) = (/ norm2 * sinmplat * sinmplon,  norm2 * cosmplon, norm1 * cosmplat * sinmplon /)
      coeff0(Rpt, iepoch, 5) = (/- norm2 * cosmplat,                      0D0, norm1 * sinmplat         /)

      ! LOOP THROUGH GEODETIC LATITUDE AND LONGITUDE GRIDPOINTS
      do ilat = 0, nlat - 1
        glat = latspace * real(ilat) + lat0
        thetag = (90D0 - dble(glat)) * dtor
        latwgt = dsin(thetag)
        ! print *, 'LINE 156'!, latwgt
        do ilon = 0, nlon - 1

          ! COMPUTE SPHERICAL HARMONICS OF CURRENT GEODETIC LONGITUDE AND
          ! LATITUDE
          glon = lonspace * real(ilon) - 180E0
          phig = dble(glon) * dtor
          call shcalc(thetag, phig)
          shg = sh
          ! COMPUTE REFERENCE (DIPOLE) MAGNETIC LATITUDE AND LONGITUDE OF
          ! CURRENT LOCATION
          xq0 = dot_product(coeff0(Rpt, iepoch, 0), shg(Rpt))
          yq0 = dot_product(coeff0(Rpt, iepoch, 1), shg(Rpt))
          zq0 = dot_product(coeff0(Rpt, iepoch, 2), shg(Rpt))
          ! LOOP THROUGH ALTITUDE GRIDPOINTS
          do ialt = 1, nalt
            ! COMPUTE QD LATITUDE AND LONGITUDE OF CURRENT LOCATION
            alt = sngl(altgrid(ialt))
            call apex(epochgrid(iepoch), igrffilein, glat, glon, alt, a, alat, alon, dum1, dum2, dum3, dum4, dum5)
            cosqlat = dsqrt( (Re + altgrid(ialt)) / (Re + Req * (dble(a) - 1D0)) )
            if (cosqlat > 1D0) cosqlat = 1D0
            phiq = dble(alon) * dtor
            thetaq = dasin(cosqlat)
            if (alat < 0) thetaq = pi - thetaq
            ! COMPUTE RESIDUAL QD COORDINATES (FITTING DATA FOR GEODETIC TO QD
            ! TRANSFORMATION)
            xq = cosqlat * dcos(phiq) - xq0
            yq = cosqlat * dsin(phiq) - yq0
            zq = dcos(thetaq)       - zq0
            ! COMPUTE SPHERICAL HARMONICS OF CURRENT QD LATITUDE AND LONGITUDE
            call shcalc(thetaq, phiq)
            ! COMPUTE RESIDUAL GD COORDINATES (FITTING DATA FOR QD TO GEODETIC
            ! TRANSFORMATION)
            xg = latwgt * dcos(phig) - dot_product(coeff0(Rpt, iepoch, 3), sh(Rpt))
            yg = latwgt * dsin(phig) - dot_product(coeff0(Rpt, iepoch, 4), sh(Rpt))
            zg = dcos(thetag)      - dot_product(coeff0(Rpt, iepoch, 5), sh(Rpt))
            ! COMPUTE BASIS FUNCTIONS FOR GD-TO-QD AND QD-TO-GD TRANSFORMATIONS
            do l = 1, lmax
              do ish = 0, ntermsh - 1
                iterm = ntermsh * (l - 1) + ish
                Fg(iterm) = altfn(l, ialt) * shg(ish)
                Fq(iterm) = altfn(l, ialt) * sh(ish)
              end do
            end do
            ! UPDATE LEAST-SQUARES ARRAYS
            do i = 0, nterm1 - 1
              Dxq(i) = Dxq(i) + Fg(i) * latwgt * xq
              Dyq(i) = Dyq(i) + Fg(i) * latwgt * yq
              Dzq(i) = Dzq(i) + Fg(i) * latwgt * zq
              Dxg(i) = Dxg(i) + Fq(i) * latwgt * xg
              Dyg(i) = Dyg(i) + Fq(i) * latwgt * yg
              Dzg(i) = Dzg(i) + Fq(i) * latwgt * zg
              do j = i, nterm1 - 1
                Gg(i, j) = Gg(i, j) + Fg(j) * latwgt * Fg(i)
                Gq(i, j) = Gq(i, j) + Fq(j) * latwgt * Fq(i)
              end do
            end do
          end do
        end do
      end do
      ! COMPUTE LEAST-SQUARES SOLUTIONS USING CHOLESKY DECOMPOSITION
      call choldc(Gg, nterm1, nterm1, cfdiag)
      call cholsl(Gg, nterm1, nterm1, cfdiag, Dxq, coefftemp)

      coeff0(ntermsh:nterm - 1, iepoch, 0) = coefftemp
      call cholsl(Gg, nterm1, nterm1, cfdiag, Dyq, coefftemp)
      coeff0(ntermsh:nterm - 1, iepoch, 1) = coefftemp
      call cholsl(Gg, nterm1, nterm1, cfdiag, Dzq, coefftemp)
      coeff0(ntermsh:nterm - 1, iepoch, 2) = coefftemp

      call choldc(Gq, nterm1, nterm1, cfdiag)
      call cholsl(Gq, nterm1, nterm1, cfdiag, Dxg, coefftemp)
      coeff0(ntermsh:nterm - 1, iepoch, 3) = coefftemp
      call cholsl(Gq, nterm1, nterm1, cfdiag, Dyg, coefftemp)
      coeff0(ntermsh:nterm - 1, iepoch, 4) = coefftemp
      call cholsl(Gq, nterm1, nterm1, cfdiag, Dzg, coefftemp)
      coeff0(ntermsh:nterm - 1, iepoch, 5) = coefftemp
    end do

    ! WRITE COEFFICIENTS TO OUTPUT FILE
    open(unit=iun, file=trim(datafilein), form='unformatted')
    write(iun) nepoch, nmax, mmax, lmax, nterm
    write(iun) epochgrid, coeff0
    close(iun)

    return

end subroutine makeapxsh

! ******************************************************************************

subroutine choldc(a, n, np, p)

    integer(4), intent(in)      :: n, np
    real(8), intent(inout)      :: a(1:np, 1:np)
    real(8), intent(out)        :: p(1:n)
    integer(4)                  :: i, j, k
    real(8)                     :: sum

    do i = 1, n
      do j = i, n
        sum = a(i, j)
        do k = i - 1, 1, - 1
          sum = sum - a(i, k) * a(j, k)
        end do
        if (i == j) then
          p(i) = dsqrt(sum)
        else
          a(j, i) = sum / p(i)
        end if
      end do
    end do
    return

end subroutine choldc

! ******************************************************************************

subroutine cholsl(a, n, np, p, b, x)

    integer(4), intent(in)      :: n, np
    real(8), intent(in)         :: a(1:np, 1:np), p(1:n), b(1:n)
    real(8), intent(out)        :: x(1:n)
    integer(4)                  :: i, k
    real(8)                     :: sum

    do i = 1, n
      sum = b(i)
      do k = i - 1, 1, - 1
        sum = sum - a(i, k) * x(k)
      end do
      x(i) = sum / p(i)
    end do
    do i = n, 1, - 1
      sum = x(i)
      do k = i + 1, n
        sum = sum - a(k, i) * x(k)
      end do
      x(i) = sum / p(i)
    end do
    return

end subroutine cholsl

! ******************************************************************************
