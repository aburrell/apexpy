! ******************************************************************************
!
! File Name: checkapexsh.f90
! Authors: John Emmert, Art Richmond
! Date: 11/13/2009
! Version: 1.0
! Description: Test driver for Quasi-Dipole coordinate conversion code package.
! References: Richmond, A. D., Ionospheric Electrodynamics Using Magnetic Apex
! Coordinates, J. Geomag. Geoelectr., 47, 191-212, 1995.
! Emmert, J. T., A. D. Richmond, and D. P. Drob, A computationally compact
! representation of Magnetic-Apex and Quasi-Dipole coordinates with smooth base
! vectors, J. Geophys. Res., 115, A08322, doi:10.1029/2010JA015326, 2010.
!
! ******************************************************************************
!
! HISTORY (blame):
!
! 25 Feb 2021: Modified by Ashton Reimer to pass IGRF coefficients file to the
! makeapxsh subroutine call.
!
! ******************************************************************************

program checkapexsh

  implicit none

  integer(4), parameter :: nepochgrid=26
  ! integer(4), parameter :: nepochgrid=3       ! For testing/debugging
  integer(4)            :: lmax=3, nmmax=6
  character(1000)       :: apexshfile='../apexpy/apexsh.dat'
  character(len=1000)   :: igrffilein='../apexpy/igrf13coeffs.txt'
  real(8)               :: epochgrid(0:nepochgrid - 1)
  real(8)               :: epoch
  real(8)               :: glat, glon, alt, hr, prec, error
  integer(4)            :: vecflag
  real(8)               :: qlat, qlon, mlat, mlon, rho
  real(8)               :: f1(1:2), f2(1:2), f
  real(8)               :: d1(1:3), d2(1:3), d3(1:3), d
  real(8)               :: e1(1:3), e2(1:3), e3(1:3)
  integer(4)            :: irho, ilat, ilon
  real(8), parameter    :: Re=6371.0088

  ! GENERATE COEFFICIENT FILE
  epochgrid = (/ 1900.0, 1905.0, 1910.0, 1915.0, 1920.0, 1925.0, 1930.0, &
                1935.0, 1940.0, 1945.0, 1950.0, 1955.0, 1960.0, 1965.0, &
                1970.0, 1975.0, 1980.0, 1985.0, 1990.0, 1995.0, 2000.0, &
                2005.0, 2010.0, 2015.0, 2020.0, 2025.0 /)
  ! epochgrid = (/1995.0,2000.0,2005.0/)      ! For testing/debugging
  call makeapxsh(apexshfile, igrffilein, epochgrid, nepochgrid, lmax, nmmax, nmmax)

  ! HEIGHT PROFILE OF QD COORDINATES
  epoch = 2005.0
  call loadapxsh(apexshfile, epoch)
  glat = - 30.0
  glon = 180.0
  vecflag = 0
  print *
  print *, 'HEIGHT PROFILE OF QUASI-DIPOLE COORDINATES'
  print '(a8,f6.1, a7,f5.1, a7,f6.1)', 'EPOCH=', epoch, ', GLAT=', glat, ', GLON=', glon
  print '(3a9)', 'ALT', 'QLAT', 'QLON'
  do irho = 0, 19
    rho = 1 - float(irho) / 20.0
    alt = Re / rho - Re
    call apxg2q(glat, glon, alt, vecflag, qlat, qlon, f1, f2, f)
    print '(f9.1,2f9.2)', alt, qlat, qlon
  end do

  ! LONGITUDE PROFILE OF QD COORDINATES AND BASE VECTORS
  epoch = 1998.3
  call loadapxsh(apexshfile, epoch)
  glat = 15.0
  alt = 400.0
  vecflag = 1
  print *
  print *, 'LONGITUDE PROFILE OF QUASI-DIPOLE COORDINATES AND BASE VECTORS'
  print '(a8,f6.1, a7,f5.1, a6,f5.1)', 'EPOCH=', epoch, ', GLAT=', glat, ', ALT=', alt
  print '(8a9)', 'GLON', 'QLAT', 'QLON', 'F1(1)', 'F1(2)', 'F2(1)', 'F2(2)', 'F'
  do ilon = - 180, 180, 20
    glon = float(ilon)
    call apxg2q(glat, glon, alt, vecflag, qlat, qlon, f1, f2, f)
    print '(f9.1,2f9.2,5f9.4)', glon, qlat, qlon, f1, f2, f
  end do

  ! LATITUDE PROFILE OF MA COORDINATES AND BASE VECTORS
  epoch = 1998.3
  call loadapxsh(apexshfile, epoch)
  glon = - 60.0
  alt = 250.0
  hr = 110.0
  vecflag = 1
  print *
  print *, 'LATITUDE PROFILE OF MODIFIED APEX COORDINATES AND BASE VECTORS'
  print '(a8,f6.1,a7,f5.1,a6,f5.1,a5,f5.1)', 'EPOCH=', epoch,', GLON=', glon,', ALT=', alt,', HR=', hr
  print '(13a9)', 'GLAT', 'QLAT', 'MLAT', 'MLON', 'D1(1)', 'D1(2)', 'D1(3)', &
                 'D2(1)', 'D2(2)', 'D2(3)', 'D3(1)', 'D3(2)', 'D3(3)'
  do ilat = - 90, 90, 10
    glat = float(ilat)
    call apxg2all(glat, glon, alt, hr, vecflag, qlat, qlon, mlat, mlon, f1, f2, f, d1, d2, d3, d, e1, e2, e3)
    print '(f9.1,3f9.2,9f9.4)', glat, qlat, mlat, mlon, d1, d2, d3
  end do

  ! LATITUDE PROFILE OF QD TO GEODETIC CONVERSION
  epoch = 2002.0
  call loadapxsh(apexshfile, epoch)
  qlon = 90.0
  alt = 1000.0
  prec = 1E-6
  print *
  print *, 'LATITUDE PROFILE OF QUASI-DIPOLE TO GEODETIC CONVERSION'
  print '(a8,f6.1, a7,f5.1, a6,f6.1)', 'EPOCH=', epoch, ', QLON=', qlon, ', ALT=', alt
  print '(3a9,a11)', 'QLAT', 'GLAT', 'GLON', 'ERROR'
  do ilat = - 90, 90, 10
    qlat = float(ilat)
    call apxq2g(qlat, qlon, alt, prec, glat, glon, error)
    print '(f9.1,2f9.2,e11.3)', qlat, glat, glon, error
  end do

end program checkapexsh


! ******************************************************************************
! TEST OUTPUT
! ******************************************************************************

! COMPUTING COEFFICIENTS FOR EPOCH 1995.0
! COMPUTING COEFFICIENTS FOR EPOCH 2000.0
! COMPUTING COEFFICIENTS FOR EPOCH 2005.0
!
! HEIGHT PROFILE OF QUASI-DIPOLE COORDINATES
! EPOCH=2005.0, GLAT=-30.0, GLON= 180.0
! ALT     QLAT     QLON
! 0.0   -34.40  -101.14
! 335.3   -34.50  -101.16
! 707.9   -34.58  -101.18
! 1124.3   -34.64  -101.21
! 1592.8   -34.69  -101.23
! 2123.7   -34.73  -101.26
! 2730.4   -34.74  -101.30
! 3430.5   -34.74  -101.33
! 4247.3   -34.72  -101.37
! 5212.6   -34.68  -101.41
! 6371.0   -34.62  -101.46
! 7786.8   -34.53  -101.51
! 9556.5   -34.43  -101.56
! 11831.9   -34.30  -101.62
! 14865.7   -34.15  -101.68
! 19113.0   -33.97  -101.75
! 25484.0   -33.77  -101.82
! 36102.4   -33.54  -101.90
! 57339.1   -33.29  -101.98
! 121049.2   -33.01  -102.07
!
! LONGITUDE PROFILE OF QUASI-DIPOLE COORDINATES AND BASE VECTORS
! EPOCH=1998.3, GLAT= 15.0, ALT=400.0
! GLON     QLAT     QLON    F1(1)    F1(2)    F2(1)    F2(2)        F
! -180.0    11.79  -110.62   0.9435  -0.1748   0.1410   0.9806   0.9498
! -160.0    14.72   -91.07   0.9509  -0.1321   0.1776   0.9808   0.9561
! -140.0    17.30   -71.48   0.9868  -0.1426   0.1966   0.9677   0.9829
! -120.0    20.30   -51.76   1.0280  -0.1689   0.1765   0.9640   1.0208
! -100.0    23.94   -31.21   1.0350  -0.2011   0.0932   1.0248   1.0794
! -80.0    26.75    -7.84   0.9787  -0.0306  -0.0457   1.1445   1.1188
! -60.0    23.42    16.40   0.9518   0.3743  -0.1901   1.0682   1.0879
! -40.0    14.30    36.50   1.0290   0.4896  -0.2621   0.9166   1.0715
! -20.0     6.95    54.65   1.1146   0.2421  -0.2115   0.9565   1.1173
! 0.0     4.67    73.95   1.1395   0.0195  -0.0822   1.0134   1.1563
! 20.0     5.40    93.41   1.1364  -0.0829   0.0104   0.9917   1.1277
! 40.0     7.33   112.65   1.1104  -0.0949   0.0199   0.9937   1.1052
! 60.0     8.21   132.26   1.0892   0.0133  -0.0156   1.0110   1.1014
! 80.0     7.37   152.06   1.0976   0.0442  -0.0247   1.0148   1.1149
! 100.0     7.19   171.76   1.0892  -0.0205  -0.0120   1.0059   1.0954
! 120.0     7.44  -168.61   1.0564   0.0091  -0.0149   1.0084   1.0654
! 140.0     7.13  -148.98   1.0277  -0.0038   0.0134   0.9980   1.0257
! 160.0     8.48  -129.79   0.9835  -0.1394   0.0825   0.9729   0.9683
! 180.0    11.79  -110.62   0.9435  -0.1748   0.1410   0.9806   0.9498
!
! LATITUDE PROFILE OF MODIFIED APEX COORDINATES AND BASE VECTORS
! EPOCH=1998.3, GLON=-60.0, ALT=250.0, HR=110.0
! GLAT     QLAT     MLAT     MLON    D1(1)    D1(2)    D1(3)    D2(1)    D2(2)    D2(3)    D3(1)    D3(2)    D3(3)
! -90.0   -74.17   -74.35    17.73   0.9266  -0.1168  -0.1031   0.0952   0.9352  -0.2560   0.1507   0.2714   1.0474
! -80.0   -64.46   -64.75    14.51   0.8894  -0.1047  -0.0993   0.0060   0.9026  -0.3419   0.1665   0.4029   1.0666
! -70.0   -54.90   -55.32    12.28   0.8433  -0.1115  -0.0745  -0.0581   0.8599  -0.4309   0.1688   0.5535   1.0819
! -60.0   -45.50   -46.09    10.47   0.8005  -0.1080  -0.0360  -0.1035   0.8034  -0.5203   0.1459   0.7206   1.0836
! -50.0   -36.25   -37.07     9.15   0.7723  -0.0790   0.0051  -0.1309   0.7277  -0.6109   0.0843   0.8919   1.0442
! -40.0   -27.11   -28.28     8.51   0.7676  -0.0225   0.0368  -0.1403   0.6227  -0.7073  -0.0136   1.0449   0.9225
! -30.0   -18.01   -19.80     8.65   0.7902   0.0511   0.0502  -0.1270   0.4712  -0.8096  -0.1185   1.1539   0.6902
! -20.0    -8.91   -12.20     9.59   0.8373   0.1248   0.0423  -0.0811   0.2577  -0.8994  -0.1961   1.1934   0.3597
! -10.0     0.21     8.36    11.19   0.8999   0.1811   0.0159   0.0023  -0.0064  -0.9392  -0.2287   1.1370  -0.0084
! 0.0     9.38    12.54    13.23   0.9645   0.2079  -0.0219   0.1072  -0.2730  -0.9059  -0.2211   0.9917  -0.3251
! 10.0    18.64    20.37    15.43   1.0171   0.2021  -0.0615   0.2036  -0.4977  -0.8165  -0.1943   0.8123  -0.5436
! 20.0    28.10    29.22    17.54   1.0477   0.1713  -0.0932   0.2706  -0.6717  -0.7009  -0.1662   0.6454  -0.6827
! 30.0    37.84    38.62    19.39   1.0532   0.1314  -0.1102   0.3041  -0.8011  -0.5739  -0.1444   0.5036  -0.7795
! 40.0    47.85    48.40    21.00   1.0374   0.1017  -0.1107   0.3136  -0.8870  -0.4415  -0.1293   0.3827  -0.8608
! 50.0    57.94    58.32    22.67   1.0076   0.0979  -0.0981   0.3171  -0.9268  -0.3128  -0.1184   0.2767  -0.9399
! 60.0    67.87    68.12    25.25   0.9700   0.1316  -0.0784   0.3397  -0.9249  -0.2020  -0.1070   0.1829  -1.0176
! 70.0    77.40    77.54    31.59   0.9120   0.2329  -0.0543   0.4312  -0.8861  -0.1211  -0.0910   0.1038  -1.0831
! 80.0    85.75    85.80    66.16   0.5051   0.6823  -0.0049   0.8681  -0.5844  -0.0748  -0.0681   0.0423  -1.1210
! 90.0    82.74    82.82   170.63  -0.9750   0.4002   0.0307   0.2765   0.7999  -0.0143  -0.0381  -0.0069  -1.1215
!
! LATITUDE PROFILE OF QUASI-DIPOLE TO GEODETIC CONVERSION
! EPOCH=2002.0, QLON= 90.0, ALT=1000.0
! QLAT     GLAT     GLON      ERROR
! -90.0   -74.99   123.89  0.382E-05
! -80.0   -74.05    85.32  0.547E-05
! -70.0   -67.48    58.23  0.512E-05
! -60.0   -57.66    42.60  0.233E-04
! -50.0   -45.78    32.34  0.518E-04
! -40.0   -33.10    25.28  0.538E-04
! -30.0   -21.03    20.83  0.132E-04
! -20.0   -10.04    18.32  0.627E-05
! -10.0     0.04    17.07  0.121E-05
! 0.0     9.52    16.57  0.331E-05
! 10.0    18.60    16.41  0.645E-05
! 20.0    27.45    16.21  0.000E+00
! 30.0    36.19    15.65  0.121E-05
! 40.0    44.93    14.34  0.209E-05
! 50.0    53.73    11.84  0.270E-05
! 60.0    62.57     7.25  0.427E-05
! 70.0    71.29    -1.71  0.400E-05
! 80.0    79.07   -23.82  0.592E-05
! 90.0    81.91   -81.31  0.226E-05
