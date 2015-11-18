C  FILE NAME: magfld.f

      SUBROUTINE COFRM (DATE)
C          Define the International Geomagnetic Reference Field (IGRF) as a
C          scalar potential field using a truncated series expansion with
C          Schmidt semi-normalized associated Legendre functions of degree n and
C          order m.  The polynomial coefficients are a function of time and are
C          interpolated between five year epochs or extrapolated at a constant
C          rate after the last epoch.
C
C          INPUTS:
C            DATE = yyyy.fraction (UT)
C          OUTPUTS (in comnon block MAGCOF):
C            NMAX = Maximum order of spherical harmonic coefficients used
C            GB   = Coefficients for magnetic field calculation
C            GV   = Coefficients for magnetic potential calculation
C            ICHG = Flag indicating when GB,GV have been changed in COFRM
C
C          It is fatal to supply a DATE before the first epoch.  A warning is
C          issued to Fortran unit 0 (stderr) if DATE is later than the
C          recommended limit, five years after the last epoch.
C
C          HISTORY (blame):
C          Apr 1983:  Written by Vincent B. Wickwar (Utah State Univ.) including
C          secular variation acceleration rate set to zero in case the IGRF
C          definition includes such second time derivitives.  The maximum degree
C          (n) defined was 10.
C
C          Jun 1986:  Updated coefficients adding Definitive Geomagnetic Reference
C          Field (DGRF) for 1980 and IGRF for 1985 (EOS Volume 7 Number 24).  The
C          designation DGRF means coefficients will not change in the future
C          whereas IGRF coefficients are interim pending incorporation of new
C          magnetometer data.  Common block MAG was replaced by MAGCOF, thus
C          removing variables not used in subroutine FELDG.  (Roy Barnes)
C
C          Apr 1992 (Barnes):  Added DGRF 1985 and IGRF 1990 as given in EOS
C          Volume 73 Number 16 April 21 1992.  Other changes were made so future
C          updates should:  (1) Increment NDGY; (2) Append to EPOCH the next IGRF
C          year; (3) Append the next DGRF coefficients to G1DIM and H1DIM; and (4)
C          replace the IGRF initial values (G0, GT) and rates of change indices
C          (H0, HT).
C
C          Apr 1994 (Art Richmond): Computation of GV added, for finding magnetic
C          potential.
C
C          Aug 1995 (Barnes):  Added DGRF for 1990 and IGRF for 1995, which were
C          obtained by anonymous ftp to geomag.gsfc.nasa.gov (cd pub, mget table*)
C          as per instructions from Bob Langel (langel@geomag.gsfc.nasa.gov) with
C          problems reported to baldwin@geomag.gsfc.nasa.gov.
C
C          Oct 1995 (Barnes):  Correct error in IGRF-95 G 7 6 and H 8 7 (see email
C          in folder).  Also found bug whereby coefficients were not being updated
C          in FELDG when IENTY did not change so ICHG was added to flag date
C          changes.  Also, a vestigial switch (IS) was removed from COFRM; it was
C          always zero and involved 3 branch if statements in the main polynomial
C          construction loop now numbered 200.
C
C          Feb 1999 (Barnes):  Explicitly initialize GV(1) in COFRM to avoid the
C          possibility of compiler or loader options initializing memory to
C          something else (e.g., indefinite).  Also simplify the algebra in COFRM
C          with no effect on results.
C
C          Mar 1999 (Barnes):  Removed three branch if's from FELDG and changed
C          statement labels to ascending order.
C
C          Jun 1999 (Barnes):  Corrected RTOD definition in GD2CART.
C
C          May 2000 (Barnes):  Replace IGRF 1995, add IGRF 2000, and extend the
C          earlier DGRF's back to 1900.  The coefficients came from an NGDC web
C          page.  Related documentation is in $APXROOT/docs/igrf.2000.*  where
C          $APXROOT, defined by 'source envapex', is traditionally ~bozo/apex).
C
C          Mar 2004 (Barnes):  Replace 1995 and 2000 coefficients; now both are
C          DGRF.  Coefficients for 2000 are degree 13 with precision increased to
C          tenths nT and accommodating this instigated changes:  (1) degree (NMAX)
C          is now a function of epoch (NMXE) to curtail irrelevant looping over
C          unused high order terms (n > 10 in epochs before 2000) when calculating
C          GB; (2) expand coefficients data statement layout for G1D and H1D,
C          formerly G1DIM and H1DIM; (3) omit secular variation acceleration terms
C          which were always zero; (4) increase array dimensions in common block
C          MAGCOF and associated arrays G and H in FELDG; (5) change earth's shape
C          in CONVRT from the IAU-1966 to the WGS-1984 spheroid; (6) eliminate
C          reference to 'definitive' in variables in COFRM which were not always
C          definitive; (7) change G to GB in COFRM s.t. arrays GB and GV in common
C          block MAGCOF are consistently named in all subroutines; (8) remove
C          unused constants in all five subroutines.  See EOS Volume 84 Number 46
C          November 18 2003, www.ngdc.noaa.gov/IAGA/vmod/igrf.html or local files
C          $APXROOT/docs/igrf.2004.*
C
C          Sept. 2005 (Maute): update with IGRF10 from 
C          http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html use script 
C          ~maute/apex.d/apex_update/igrf2f Note that the downloaded file the start 
C          column of the year in the first line has to be before the start of each 
C          number in the same column
C   
C          Jan. 2010 (Maute) update with IGRF11 (same instructions as Sep. 2005
C          comment

      DOUBLE PRECISION F,F0
      COMMON /MAGCOF/ NMAX,GB(255),GV(225),ICHG
      DATA ICHG /-99999/
 
C          NEPO = Number of epochs
C          NGH  = Single dimensioned array size of 2D version (GYR or HYR)
C          NGHT = Single dimensioned array size of 2D version (GT  or HT)
      PARAMETER (NEPO = 24, NGH = 225*NEPO, NGHT = 225)
      DIMENSION GYR(15,15,NEPO), HYR(15,15,NEPO), EPOCH(NEPO),
     +          GT (15,15),      HT (15,15),       NMXE(NEPO),
     +          GY1D(NGH),       HY1D(NGH),
     +          GT1D(NGHT),      HT1D(NGHT)
      EQUIVALENCE (GYR(1,1,1),GY1D(1)), (HYR(1,1,1),HY1D(1)),
     +            (GT (1,1),  GT1D(1)), (HT (1,1),  HT1D(1))

      SAVE DATEL, EPOCH, NMXE, GYR, HYR, GT, HT, GY1D, HY1D, GT1D, HT1D
      DATA DATEL /-999./,
     +     EPOCH / 1900, 1905, 1910, 1915, 1920, 1925, 1930, 1935, 1940,
     +             1945, 1950, 1955, 1960, 1965, 1970, 1975, 1980, 1985,
     +             1990, 1995, 2000, 2005, 2010, 2015/,
     +     NMXE  /   10,   10,   10,   10,   10,   10,   10,   10,   10,
     +               10,   10,   10,   10,   10,   10,   10,   10,   10,
     +               10,   10,   13,   13,   13,   13/

C          g(n,m) for 1900
C          Fields across a line are (degree) n=1,13; lines are (order) m=0,13 as indicated
C          in column 6; e.g., for 1965 g(n=3,m=0) = 1297 or g(n=6,m=6) = -111
C
C           1       2       3      4      5      6      7      8      9
C                                        10     11     12     13          (n)
      DATA (GY1D(I),I=1,145) /0,
     O  -31543,   -677,  1022,   876,  -184,    63,    70,    11,     8,
     +                                   -3,     0,     0,     0,   2*0,
     1   -2298,   2905, -1469,   628,   328,    61,   -55,     8,    10,
     +                                   -4,     0,     0,     0,   3*0,
     2             924,  1256,   660,   264,   -11,     0,    -4,     1,
     +                                    2,     0,     0,     0,   4*0,
     3                    572,  -361,     5,  -217,    34,    -9,   -11,
     +                                   -5,     0,     0,     0,   5*0,
     4                           134,   -86,   -58,   -41,     1,    12,
     +                                   -2,     0,     0,     0,   6*0,
     5                                  -16,    59,   -21,     2,     1,
     +                                    6,     0,     0,     0,   7*0,
     6                                         -90,    18,    -9,    -2,
     +                                    4,     0,     0,     0,   8*0,
     7                                                  6,     5,     2,
     +                                    0,     0,     0,     0,   9*0,
     8                                                         8,    -1,
     +                                    2,     0,     0,     0,  10*0,
     9                                                               -1/
      DATA (GY1D(I),I=146,225) /
     +                                    2,     0,     0,     0,  11*0,
     O                                    0,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          g(n,m) for 1905
      DATA (GY1D(I),I=226,370) /0,
     O  -31464,   -728,  1037,   880,  -192,    62,    70,    11,     8,
     +                                   -3,     0,     0,     0,   2*0,
     1   -2298,   2928, -1494,   643,   328,    60,   -54,     8,    10,
     +                                   -4,     0,     0,     0,   3*0,
     2            1041,  1239,   653,   259,   -11,     0,    -4,     1,
     +                                    2,     0,     0,     0,   4*0,
     3                    635,  -380,    -1,  -221,    33,    -9,   -11,
     +                                   -5,     0,     0,     0,   5*0,
     4                           146,   -93,   -57,   -41,     1,    12,
     +                                   -2,     0,     0,     0,   6*0,
     5                                  -26,    57,   -20,     2,     1,
     +                                    6,     0,     0,     0,   7*0,
     6                                         -92,    18,    -8,    -2,
     +                                    4,     0,     0,     0,   8*0,
     7                                                  6,     5,     2,
     +                                    0,     0,     0,     0,   9*0,
     8                                                         8,     0,
     +                                    2,     0,     0,     0,  10*0,
     9                                                               -1/
      DATA (GY1D(I),I=371,450) /
     +                                    2,     0,     0,     0,  11*0,
     O                                    0,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          g(n,m) for 1910
      DATA (GY1D(I),I=451,595) /0,
     O  -31354,   -769,  1058,   884,  -201,    62,    71,    11,     8,
     +                                   -3,     0,     0,     0,   2*0,
     1   -2297,   2948, -1524,   660,   327,    58,   -54,     8,    10,
     +                                   -4,     0,     0,     0,   3*0,
     2            1176,  1223,   644,   253,   -11,     1,    -4,     1,
     +                                    2,     0,     0,     0,   4*0,
     3                    705,  -400,    -9,  -224,    32,    -9,   -11,
     +                                   -5,     0,     0,     0,   5*0,
     4                           160,  -102,   -54,   -40,     1,    12,
     +                                   -2,     0,     0,     0,   6*0,
     5                                  -38,    54,   -19,     2,     1,
     +                                    6,     0,     0,     0,   7*0,
     6                                         -95,    18,    -8,    -2,
     +                                    4,     0,     0,     0,   8*0,
     7                                                  6,     5,     2,
     +                                    0,     0,     0,     0,   9*0,
     8                                                         8,     0,
     +                                    2,     0,     0,     0,  10*0,
     9                                                               -1/
      DATA (GY1D(I),I=596,675) /
     +                                    2,     0,     0,     0,  11*0,
     O                                    0,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          g(n,m) for 1915
      DATA (GY1D(I),I=676,820) /0,
     O  -31212,   -802,  1084,   887,  -211,    61,    72,    11,     8,
     +                                   -3,     0,     0,     0,   2*0,
     1   -2306,   2956, -1559,   678,   327,    57,   -54,     8,    10,
     +                                   -4,     0,     0,     0,   3*0,
     2            1309,  1212,   631,   245,   -10,     2,    -4,     1,
     +                                    2,     0,     0,     0,   4*0,
     3                    778,  -416,   -16,  -228,    31,    -9,   -11,
     +                                   -5,     0,     0,     0,   5*0,
     4                           178,  -111,   -51,   -38,     2,    12,
     +                                   -2,     0,     0,     0,   6*0,
     5                                  -51,    49,   -18,     3,     1,
     +                                    6,     0,     0,     0,   7*0,
     6                                         -98,    19,    -8,    -2,
     +                                    4,     0,     0,     0,   8*0,
     7                                                  6,     6,     2,
     +                                    0,     0,     0,     0,   9*0,
     8                                                         8,     0,
     +                                    1,     0,     0,     0,  10*0,
     9                                                               -1/
      DATA (GY1D(I),I=821,900) /
     +                                    2,     0,     0,     0,  11*0,
     O                                    0,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          g(n,m) for 1920
      DATA (GY1D(I),I=901,1045) /0,
     O  -31060,   -839,  1111,   889,  -221,    61,    73,    11,     8,
     +                                   -3,     0,     0,     0,   2*0,
     1   -2317,   2959, -1600,   695,   326,    55,   -54,     7,    10,
     +                                   -4,     0,     0,     0,   3*0,
     2            1407,  1205,   616,   236,   -10,     2,    -3,     1,
     +                                    2,     0,     0,     0,   4*0,
     3                    839,  -424,   -23,  -233,    29,    -9,   -11,
     +                                   -5,     0,     0,     0,   5*0,
     4                           199,  -119,   -46,   -37,     2,    12,
     +                                   -2,     0,     0,     0,   6*0,
     5                                  -62,    44,   -16,     4,     1,
     +                                    6,     0,     0,     0,   7*0,
     6                                        -101,    19,    -7,    -2,
     +                                    4,     0,     0,     0,   8*0,
     7                                                  6,     6,     2,
     +                                    0,     0,     0,     0,   9*0,
     8                                                         8,     0,
     +                                    1,     0,     0,     0,  10*0,
     9                                                               -1/
      DATA (GY1D(I),I=1046,1125) /
     +                                    3,     0,     0,     0,  11*0,
     O                                    0,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          g(n,m) for 1925
      DATA (GY1D(I),I=1126,1270) /0,
     O  -30926,   -893,  1140,   891,  -230,    61,    73,    11,     8,
     +                                   -3,     0,     0,     0,   2*0,
     1   -2318,   2969, -1645,   711,   326,    54,   -54,     7,    10,
     +                                   -4,     0,     0,     0,   3*0,
     2            1471,  1202,   601,   226,    -9,     3,    -3,     1,
     +                                    2,     0,     0,     0,   4*0,
     3                    881,  -426,   -28,  -238,    27,    -9,   -11,
     +                                   -5,     0,     0,     0,   5*0,
     4                           217,  -125,   -40,   -35,     2,    12,
     +                                   -2,     0,     0,     0,   6*0,
     5                                  -69,    39,   -14,     4,     1,
     +                                    6,     0,     0,     0,   7*0,
     6                                        -103,    19,    -7,    -2,
     +                                    4,     0,     0,     0,   8*0,
     7                                                  6,     7,     2,
     +                                    0,     0,     0,     0,   9*0,
     8                                                         8,     0,
     +                                    1,     0,     0,     0,  10*0,
     9                                                               -1/
      DATA (GY1D(I),I=1271,1350) /
     +                                    3,     0,     0,     0,  11*0,
     O                                    0,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          g(n,m) for 1930
      DATA (GY1D(I),I=1351,1495) /0,
     O  -30805,   -951,  1172,   896,  -237,    60,    74,    11,     8,
     +                                   -3,     0,     0,     0,   2*0,
     1   -2316,   2980, -1692,   727,   327,    53,   -54,     7,    10,
     +                                   -4,     0,     0,     0,   3*0,
     2            1517,  1205,   584,   218,    -9,     4,    -3,     1,
     +                                    2,     0,     0,     0,   4*0,
     3                    907,  -422,   -32,  -242,    25,    -9,   -12,
     +                                   -5,     0,     0,     0,   5*0,
     4                           234,  -131,   -32,   -34,     2,    12,
     +                                   -2,     0,     0,     0,   6*0,
     5                                  -74,    32,   -12,     5,     1,
     +                                    6,     0,     0,     0,   7*0,
     6                                        -104,    18,    -6,    -2,
     +                                    4,     0,     0,     0,   8*0,
     7                                                  6,     8,     3,
     +                                    0,     0,     0,     0,   9*0,
     8                                                         8,     0,
     +                                    1,     0,     0,     0,  10*0,
     9                                                               -2/
      DATA (GY1D(I),I=1496,1575) /
     +                                    3,     0,     0,     0,  11*0,
     O                                    0,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          g(n,m) for 1935
      DATA (GY1D(I),I=1576,1720) /0,
     O  -30715,  -1018,  1206,   903,  -241,    59,    74,    11,     8,
     +                                   -3,     0,     0,     0,   2*0,
     1   -2306,   2984, -1740,   744,   329,    53,   -53,     7,    10,
     +                                   -4,     0,     0,     0,   3*0,
     2            1550,  1215,   565,   211,    -8,     4,    -3,     1,
     +                                    2,     0,     0,     0,   4*0,
     3                    918,  -415,   -33,  -246,    23,    -9,   -12,
     +                                   -5,     0,     0,     0,   5*0,
     4                           249,  -136,   -25,   -33,     1,    11,
     +                                   -2,     0,     0,     0,   6*0,
     5                                  -76,    25,   -11,     6,     1,
     +                                    6,     0,     0,     0,   7*0,
     6                                        -106,    18,    -6,    -2,
     +                                    4,     0,     0,     0,   8*0,
     7                                                  6,     8,     3,
     +                                    0,     0,     0,     0,   9*0,
     8                                                         7,     0,
     +                                    2,     0,     0,     0,  10*0,
     9                                                               -2/
      DATA (GY1D(I),I=1721,1800) /
     +                                    3,     0,     0,     0,  11*0,
     O                                    0,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          g(n,m) for 1940
      DATA (GY1D(I),I=1801,1945) /0,
     O  -30654,  -1106,  1240,   914,  -241,    57,    74,    11,     8,
     +                                   -3,     0,     0,     0,   2*0,
     1   -2292,   2981, -1790,   762,   334,    54,   -53,     7,    10,
     +                                   -4,     0,     0,     0,   3*0,
     2            1566,  1232,   550,   208,    -7,     4,    -3,     1,
     +                                    2,     0,     0,     0,   4*0,
     3                    916,  -405,   -33,  -249,    20,   -10,   -12,
     +                                   -5,     0,     0,     0,   5*0,
     4                           265,  -141,   -18,   -31,     1,    11,
     +                                   -2,     0,     0,     0,   6*0,
     5                                  -76,    18,    -9,     6,     1,
     +                                    6,     0,     0,     0,   7*0,
     6                                        -107,    17,    -5,    -2,
     +                                    4,     0,     0,     0,   8*0,
     7                                                  5,     9,     3,
     +                                    0,     0,     0,     0,   9*0,
     8                                                         7,     1,
     +                                    2,     0,     0,     0,  10*0,
     9                                                               -2/
      DATA (GY1D(I),I=1946,2025) /
     +                                    3,     0,     0,     0,  11*0,
     O                                    0,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          g(n,m) for 1945
      DATA (GY1D(I),I=2026,2170) /0,
     O  -30594,  -1244,  1282,   944,  -253,    59,    70,    13,     5,
     +                                   -3,     0,     0,     0,   2*0,
     1   -2285,   2990, -1834,   776,   346,    57,   -40,     7,   -21,
     +                                   11,     0,     0,     0,   3*0,
     2            1578,  1255,   544,   194,     6,     0,    -8,     1,
     +                                    1,     0,     0,     0,   4*0,
     3                    913,  -421,   -20,  -246,     0,    -5,   -11,
     +                                    2,     0,     0,     0,   5*0,
     4                           304,  -142,   -25,   -29,     9,     3,
     +                                   -5,     0,     0,     0,   6*0,
     5                                  -82,    21,   -10,     7,    16,
     +                                   -1,     0,     0,     0,   7*0,
     6                                        -104,    15,   -10,    -3,
     +                                    8,     0,     0,     0,   8*0,
     7                                                 29,     7,    -4,
     +                                   -1,     0,     0,     0,   9*0,
     8                                                         2,    -3,
     +                                   -3,     0,     0,     0,  10*0,
     9                                                               -4/
      DATA (GY1D(I),I=2171,2250) /
     +                                    5,     0,     0,     0,  11*0,
     O                                   -2,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          g(n,m) for 1950
      DATA (GY1D(I),I=2251,2395) /0,
     O  -30554,  -1341,  1297,   954,  -240,    54,    65,    22,     3,
     +                                   -8,     0,     0,     0,   2*0,
     1   -2250,   2998, -1889,   792,   349,    57,   -55,    15,    -7,
     +                                    4,     0,     0,     0,   3*0,
     2            1576,  1274,   528,   211,     4,     2,    -4,    -1,
     +                                   -1,     0,     0,     0,   4*0,
     3                    896,  -408,   -20,  -247,     1,    -1,   -25,
     +                                   13,     0,     0,     0,   5*0,
     4                           303,  -147,   -16,   -40,    11,    10,
     +                                   -4,     0,     0,     0,   6*0,
     5                                  -76,    12,    -7,    15,     5,
     +                                    4,     0,     0,     0,   7*0,
     6                                        -105,     5,   -13,    -5,
     +                                   12,     0,     0,     0,   8*0,
     7                                                 19,     5,    -2,
     +                                    3,     0,     0,     0,   9*0,
     8                                                        -1,     3,
     +                                    2,     0,     0,     0,  10*0,
     9                                                                8/
      DATA (GY1D(I),I=2396,2475) /
     +                                   10,     0,     0,     0,  11*0,
     O                                    3,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          g(n,m) for 1955
      DATA (GY1D(I),I=2476,2620) /0,
     O  -30500,  -1440,  1302,   958,  -229,    47,    65,    11,     4,
     +                                   -3,     0,     0,     0,   2*0,
     1   -2215,   3003, -1944,   796,   360,    57,   -56,     9,     9,
     +                                   -5,     0,     0,     0,   3*0,
     2            1581,  1288,   510,   230,     3,     2,    -6,    -4,
     +                                   -1,     0,     0,     0,   4*0,
     3                    882,  -397,   -23,  -247,    10,   -14,    -5,
     +                                    2,     0,     0,     0,   5*0,
     4                           290,  -152,    -8,   -32,     6,     2,
     +                                   -3,     0,     0,     0,   6*0,
     5                                  -69,     7,   -11,    10,     4,
     +                                    7,     0,     0,     0,   7*0,
     6                                        -107,     9,    -7,     1,
     +                                    4,     0,     0,     0,   8*0,
     7                                                 18,     6,     2,
     +                                   -2,     0,     0,     0,   9*0,
     8                                                         9,     2,
     +                                    6,     0,     0,     0,  10*0,
     9                                                                5/
      DATA (GY1D(I),I=2621,2700) /
     +                                   -2,     0,     0,     0,  11*0,
     O                                    0,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          g(n,m) for 1960
      DATA (GY1D(I),I=2701,2845) /0,
     O  -30421,  -1555,  1302,   957,  -222,    46,    67,    15,     4,
     +                                    1,     0,     0,     0,   2*0,
     1   -2169,   3002, -1992,   800,   362,    58,   -56,     6,     6,
     +                                   -3,     0,     0,     0,   3*0,
     2            1590,  1289,   504,   242,     1,     5,    -4,     0,
     +                                    4,     0,     0,     0,   4*0,
     3                    878,  -394,   -26,  -237,    15,   -11,    -9,
     +                                    0,     0,     0,     0,   5*0,
     4                           269,  -156,    -1,   -32,     2,     1,
     +                                   -1,     0,     0,     0,   6*0,
     5                                  -63,    -2,    -7,    10,     4,
     +                                    4,     0,     0,     0,   7*0,
     6                                        -113,    17,    -5,    -1,
     +                                    6,     0,     0,     0,   8*0,
     7                                                  8,    10,    -2,
     +                                    1,     0,     0,     0,   9*0,
     8                                                         8,     3,
     +                                   -1,     0,     0,     0,  10*0,
     9                                                               -1/
      DATA (GY1D(I),I=2846,2925) /
     +                                    2,     0,     0,     0,  11*0,
     O                                    0,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          g(n,m) for 1965
      DATA (GY1D(I),I=2926,3070) /0,
     O  -30334,  -1662,  1297,   957,  -219,    45,    75,    13,     8,
     +                                   -2,     0,     0,     0,   2*0,
     1   -2119,   2997, -2038,   804,   358,    61,   -57,     5,    10,
     +                                   -3,     0,     0,     0,   3*0,
     2            1594,  1292,   479,   254,     8,     4,    -4,     2,
     +                                    2,     0,     0,     0,   4*0,
     3                    856,  -390,   -31,  -228,    13,   -14,   -13,
     +                                   -5,     0,     0,     0,   5*0,
     4                           252,  -157,     4,   -26,     0,    10,
     +                                   -2,     0,     0,     0,   6*0,
     5                                  -62,     1,    -6,     8,    -1,
     +                                    4,     0,     0,     0,   7*0,
     6                                        -111,    13,    -1,    -1,
     +                                    4,     0,     0,     0,   8*0,
     7                                                  1,    11,     5,
     +                                    0,     0,     0,     0,   9*0,
     8                                                         4,     1,
     +                                    2,     0,     0,     0,  10*0,
     9                                                               -2/
      DATA (GY1D(I),I=3071,3150) /
     +                                    2,     0,     0,     0,  11*0,
     O                                    0,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          g(n,m) for 1970
      DATA (GY1D(I),I=3151,3295) /0,
     O  -30220,  -1781,  1287,   952,  -216,    43,    72,    14,     8,
     +                                   -3,     0,     0,     0,   2*0,
     1   -2068,   3000, -2091,   800,   359,    64,   -57,     6,    10,
     +                                   -3,     0,     0,     0,   3*0,
     2            1611,  1278,   461,   262,    15,     1,    -2,     2,
     +                                    2,     0,     0,     0,   4*0,
     3                    838,  -395,   -42,  -212,    14,   -13,   -12,
     +                                   -5,     0,     0,     0,   5*0,
     4                           234,  -160,     2,   -22,    -3,    10,
     +                                   -1,     0,     0,     0,   6*0,
     5                                  -56,     3,    -2,     5,    -1,
     +                                    6,     0,     0,     0,   7*0,
     6                                        -112,    13,     0,     0,
     +                                    4,     0,     0,     0,   8*0,
     7                                                 -2,    11,     3,
     +                                    1,     0,     0,     0,   9*0,
     8                                                         3,     1,
     +                                    0,     0,     0,     0,  10*0,
     9                                                               -1/
      DATA (GY1D(I),I=3296,3375) /
     +                                    3,     0,     0,     0,  11*0,
     O                                   -1,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          g(n,m) for 1975
      DATA (GY1D(I),I=3376,3520) /0,
     O  -30100,  -1902,  1276,   946,  -218,    45,    71,    14,     7,
     +                                   -3,     0,     0,     0,   2*0,
     1   -2013,   3010, -2144,   791,   356,    66,   -56,     6,    10,
     +                                   -3,     0,     0,     0,   3*0,
     2            1632,  1260,   438,   264,    28,     1,    -1,     2,
     +                                    2,     0,     0,     0,   4*0,
     3                    830,  -405,   -59,  -198,    16,   -12,   -12,
     +                                   -5,     0,     0,     0,   5*0,
     4                           216,  -159,     1,   -14,    -8,    10,
     +                                   -2,     0,     0,     0,   6*0,
     5                                  -49,     6,     0,     4,    -1,
     +                                    5,     0,     0,     0,   7*0,
     6                                        -111,    12,     0,    -1,
     +                                    4,     0,     0,     0,   8*0,
     7                                                 -5,    10,     4,
     +                                    1,     0,     0,     0,   9*0,
     8                                                         1,     1,
     +                                    0,     0,     0,     0,  10*0,
     9                                                               -2/
      DATA (GY1D(I),I=3521,3600) /
     +                                    3,     0,     0,     0,  11*0,
     O                                   -1,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          g(n,m) for 1980
      DATA (GY1D(I),I=3601,3745) /0,
     O  -29992,  -1997,  1281,   938,  -218,    48,    72,    18,     5,
     +                                   -4,     0,     0,     0,   2*0,
     1   -1956,   3027, -2180,   782,   357,    66,   -59,     6,    10,
     +                                   -4,     0,     0,     0,   3*0,
     2            1663,  1251,   398,   261,    42,     2,     0,     1,
     +                                    2,     0,     0,     0,   4*0,
     3                    833,  -419,   -74,  -192,    21,   -11,   -12,
     +                                   -5,     0,     0,     0,   5*0,
     4                           199,  -162,     4,   -12,    -7,     9,
     +                                   -2,     0,     0,     0,   6*0,
     5                                  -48,    14,     1,     4,    -3,
     +                                    5,     0,     0,     0,   7*0,
     6                                        -108,    11,     3,    -1,
     +                                    3,     0,     0,     0,   8*0,
     7                                                 -2,     6,     7,
     +                                    1,     0,     0,     0,   9*0,
     8                                                        -1,     2,
     +                                    2,     0,     0,     0,  10*0,
     9                                                               -5/
      DATA (GY1D(I),I=3746,3825) /
     +                                    3,     0,     0,     0,  11*0,
     O                                    0,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          g(n,m) for 1985
      DATA (GY1D(I),I=3826,3970) /0,
     O  -29873,  -2072,  1296,   936,  -214,    53,    74,    21,     5,
     +                                   -4,     0,     0,     0,   2*0,
     1   -1905,   3044, -2208,   780,   355,    65,   -62,     6,    10,
     +                                   -4,     0,     0,     0,   3*0,
     2            1687,  1247,   361,   253,    51,     3,     0,     1,
     +                                    3,     0,     0,     0,   4*0,
     3                    829,  -424,   -93,  -185,    24,   -11,   -12,
     +                                   -5,     0,     0,     0,   5*0,
     4                           170,  -164,     4,    -6,    -9,     9,
     +                                   -2,     0,     0,     0,   6*0,
     5                                  -46,    16,     4,     4,    -3,
     +                                    5,     0,     0,     0,   7*0,
     6                                        -102,    10,     4,    -1,
     +                                    3,     0,     0,     0,   8*0,
     7                                                  0,     4,     7,
     +                                    1,     0,     0,     0,   9*0,
     8                                                        -4,     1,
     +                                    2,     0,     0,     0,  10*0,
     9                                                               -5/
      DATA (GY1D(I),I=3971,4050) /
     +                                    3,     0,     0,     0,  11*0,
     O                                    0,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          g(n,m) for 1990
      DATA (GY1D(I),I=4051,4195) /0,
     O  -29775,  -2131,  1314,   939,  -214,    61,    77,    23,     4,
     +                                   -3,     0,     0,     0,   2*0,
     1   -1848,   3059, -2239,   780,   353,    65,   -64,     5,     9,
     +                                   -4,     0,     0,     0,   3*0,
     2            1686,  1248,   325,   245,    59,     2,    -1,     1,
     +                                    2,     0,     0,     0,   4*0,
     3                    802,  -423,  -109,  -178,    26,   -10,   -12,
     +                                   -5,     0,     0,     0,   5*0,
     4                           141,  -165,     3,    -1,   -12,     9,
     +                                   -2,     0,     0,     0,   6*0,
     5                                  -36,    18,     5,     3,    -4,
     +                                    4,     0,     0,     0,   7*0,
     6                                         -96,     9,     4,    -2,
     +                                    3,     0,     0,     0,   8*0,
     7                                                  0,     2,     7,
     +                                    1,     0,     0,     0,   9*0,
     8                                                        -6,     1,
     +                                    3,     0,     0,     0,  10*0,
     9                                                               -6/
      DATA (GY1D(I),I=4196,4275) /
     +                                    3,     0,     0,     0,  11*0,
     O                                    0,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          g(n,m) for 1995
      DATA (GY1D(I),I=4276,4420) /0,
     O  -29692,  -2200,  1335,   940,  -214,    68,    77,    25,     4,
     +                                   -3,     0,     0,     0,   2*0,
     1   -1784,   3070, -2267,   780,   352,    67,   -72,     6,     9,
     +                                   -6,     0,     0,     0,   3*0,
     2            1681,  1249,   290,   235,    68,     1,    -6,     3,
     +                                    2,     0,     0,     0,   4*0,
     3                    759,  -418,  -118,  -170,    28,    -9,   -10,
     +                                   -4,     0,     0,     0,   5*0,
     4                           122,  -166,    -1,     5,   -14,     8,
     +                                   -1,     0,     0,     0,   6*0,
     5                                  -17,    19,     4,     9,    -8,
     +                                    4,     0,     0,     0,   7*0,
     6                                         -93,     8,     6,    -1,
     +                                    2,     0,     0,     0,   8*0,
     7                                                 -2,    -5,    10,
     +                                    2,     0,     0,     0,   9*0,
     8                                                        -7,    -2,
     +                                    5,     0,     0,     0,  10*0,
     9                                                               -8/
      DATA (GY1D(I),I=4421,4500) /
     +                                    1,     0,     0,     0,  11*0,
     O                                    0,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          g(n,m) for 2000
      DATA (GY1D(I),I=4501,4645) /0,
     O-29619.4,-2267.7,1339.6, 932.3,-218.8,  72.3,  79.0,  24.4,   5.0,
     +                                 -2.6,   2.7,  -2.2,  -0.2,   2*0,
     1 -1728.2, 3068.4,-2288.0, 786.8, 351.4,  68.2, -74.0,   6.6,  9.4,
     +                                 -6.0,  -1.7,  -0.3,  -0.9,   3*0,
     2          1670.9,1252.1, 250.0, 222.3,  74.2,   0.0,  -9.2,   3.0,
     +                                  1.7,  -1.9,   0.2,   0.3,   4*0,
     3                  714.5,-403.0,-130.4,-160.9,  33.3,  -7.9,  -8.4,
     +                                 -3.1,   1.5,   0.9,   0.1,   5*0,
     4                         111.3,-168.6,  -5.9,   9.1, -16.6,   6.3,
     +                                 -0.5,  -0.1,  -0.2,  -0.4,   6*0,
     5                                -12.9,  16.9,   6.9,   9.1,  -8.9,
     +                                  3.7,   0.1,   0.9,   1.3,   7*0,
     6                                       -90.4,   7.3,   7.0,  -1.5,
     +                                  1.0,  -0.7,  -0.5,  -0.4,   8*0,
     7                                               -1.2,  -7.9,   9.3,
     +                                  2.0,   0.7,   0.3,   0.7,   9*0,
     8                                                      -7.0,  -4.3,
     +                                  4.2,   1.7,  -0.3,  -0.4,  10*0,
     9                                                            -8.2/
      DATA (GY1D(I),I=4646,4725) /
     +                                  0.3,   0.1,  -0.4,   0.3,  11*0,
     O                                 -1.1,   1.2,  -0.1,  -0.1,  12*0,
     1                                         4.0,  -0.2,   0.4,  13*0,
     2                                               -0.4,   0.0,  14*0,
     3                                                       0.1,  16*0/
C          g(n,m) for 2005
      DATA (GY1D(I),I=4726,4870) /0,
     O-29554.63,-2337.24,1336.30,920.55,-227.00, 73.60, 79.88, 24.80,
     +                         5.58, -2.17,  2.95, -2.15, -0.16,   2*0,
     1-1669.05,3047.69,-2305.83,797.96,354.41, 69.56,-74.46,  7.62,
     +                          9.76, -6.12, -1.60, -0.29, -0.88,   3*0,
     2         1657.76,1246.39,210.65,208.95, 76.74, -1.65,-11.73,
     +                          3.58,   1.42, -1.88,  0.21,  0.30,  4*0,
     3                 672.51,-379.86,-136.54,-151.34, 38.73, -6.88,
     +                       -6.94,   -2.35,  1.44,  0.89,  0.28,   5*0,
     4                      100.00,-168.05,-14.58, 12.30,-18.11,  5.01,
     +                     	    -0.15, -0.31, -0.38, -0.43,   6*0,
     5                     	   -13.55, 14.58,  9.37, 10.17,-10.76,
     +                     	     3.06,  0.29,  0.96,  1.18,   7*0,
     6                     		  -86.36,  5.42,  9.36, -1.25,
     +                     	     0.29, -0.79, -0.30, -0.37,   8*0,
     7                     			   1.94,-11.25,  8.76,
     +                     	     2.06,  0.53,  0.46,  0.75,   9*0,
     8                     				 -4.87, -6.66,
     +                     	     3.77,  1.80, -0.35, -0.26,  10*0,
     9                                                          -9.22/
      DATA (GY1D(I),I=4871,4950) /
     +                                -0.21,  0.16, -0.36,  0.35,  11*0,
     O                                -2.09,  0.96,  0.08, -0.05,  12*0,
     1                                        3.99, -0.49,  0.41,  13*0,
     2                                              -0.08, -0.10,  14*0,
     3                                                     -0.18,  16*0/
C          g(n,m) for 2010
      DATA (GY1D(I),I=4951,5095) /0,
     O -29496.57,-2396.06,1339.85,912.66, -230.87,72.78,80.44,24.41,5.5,
     +                                 -1.94, 3.05, -2.12, -0.09,   2*0,
     1 -1586.42, 3026.34, -2326.54,808.97,357.29,68.69, -75.0,8.21,9.45,
     +                                -6.24, -1.48, -0.21, -0.89,   3*0,
     2       1668.17, 1232.1, 166.58, 200.26, 75.92, -4.55, -14.5, 3.45,
     +                                    0.89, -2.03, 0.3, 0.31,   4*0,
     3            633.73, -356.83, -141.05, -141.4, 45.24, -5.59, -5.27,
     +                                   -1.07, 1.65, 1.04, 0.42,   5*0,
     4                        89.4, -163.17, -22.83, 14.0, -19.34, 3.13,
     +                                -0.16, -0.51, -0.63, -0.45,   6*0,
     5                                -8.03, 13.1, 10.46, 11.61, -12.38,
     +                                    2.45, 0.54, 0.95, 1.08,   7*0,
     6                                       -78.09, 1.64, 10.85, -0.76,
     +                                -0.33, -0.79, -0.11, -0.31,   8*0,
     7                                               4.92, -14.05, 8.43,
     +                                    2.13, 0.37, 0.52, 0.78,   9*0,
     8                                                     -3.54, -8.42,
     +                                  3.09, 1.79, -0.39, -0.18,  10*0,
     9                                                           -10.08/
      DATA (GY1D(I),I=5096,5175) /
     +                                  -1.03, 0.12, -0.37, 0.38,  11*0,
     O                                    -2.8, 0.75, 0.21, 0.02,  12*0,
     1                                         3.75, -0.77, 0.42,  13*0,
     2                                               0.04, -0.26,  14*0,
     3                                                     -0.26,  16*0/
C          g(n,m) for 2015
      DATA (GY1D(I),I=5176,5320) /0,
     O  -29442.0, -2445.1, 1350.7, 907.6, -232.6, 70.0, 81.6, 24.2, 5.4,
     +                                      -1.9, 3.1, -1.9, 0.0,   2*0,
     1    -1501.0, 3012.9, -2352.3, 813.7, 360.1, 67.7, -76.1, 8.8, 8.8,
     +                                    -6.3, -1.5, -0.2, -0.9,   3*0,
     2             1676.7, 1225.6, 120.4, 192.4, 72.7, -6.8, -16.9, 3.1,
     +                                       0.1, -2.3, 0.4, 0.4,   4*0,
     3                  582.0, -334.9, -140.9, -129.9, 51.8, -3.2, -3.3,
     +                                        0.5, 2.0, 1.2, 0.5,   5*0,
     4                            70.4, -157.5, -28.9, 15.0, -20.6, 0.7,
     +                                   -0.5, -0.8, -0.8,  -0.5,   6*0,
     5                                      4.1, 13.2, 9.4, 13.4, -13.3,
     +                                        1.8, 0.6, 0.9, 1.0,   7*0,
     6                                          -70.9, -2.8, 11.7, -0.1,
     +                                     -0.7, -0.7, 0.1, -0.2,   8*0,
     7                                                  6.8, -15.9, 8.7,
     +                                        2.1, 0.2, 0.5, 0.8,   9*0,
     8                                                       -2.0, -9.1,
     +                                      2.4, 1.7, -0.3, -0.1,  10*0,
     9                                                            -10.5/
      DATA (GY1D(I),I=5321,5400) /
     +                                     -1.8, -0.2, -0.4, 0.3,  11*0,
     O                                       -3.6, 0.4, 0.2, 0.1,  12*0,
     1                                            3.5, -0.9, 0.5,  13*0,
     2                                                 0.0, -0.4,  14*0,
     3                                                      -0.3,  16*0/
C          h(n,m) for 1900
      DATA (HY1D(I),I=1,145) /16*0,
     1    5922,  -1061,  -330,   195,  -210,    -9,   -45,     8,   -20,
     +                                    2,     0,     0,     0,   3*0,
     2            1121,     3,   -69,    53,    83,   -13,   -14,    14,
     +                                    1,     0,     0,     0,   4*0,
     3                    523,  -210,   -33,     2,   -10,     7,     5,
     +                                    2,     0,     0,     0,   5*0,
     4                           -75,  -124,   -35,    -1,   -13,    -3,
     +                                    6,     0,     0,     0,   6*0,
     5                                    3,    36,    28,     5,    -2,
     +                                   -4,     0,     0,     0,   7*0,
     6                                         -69,   -12,    16,     8,
     +                                    0,     0,     0,     0,   8*0,
     7                                                -22,    -5,    10,
     +                                   -2,     0,     0,     0,   9*0,
     8                                                       -18,    -2,
     +                                    4,     0,     0,     0,  10*0,
     9                                                                2/
      DATA (HY1D(I),I=146,225) /
     +                                    0,     0,     0,     0,  11*0,
     O                                   -6,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          h(n,m) for 1905
      DATA (HY1D(I),I=226,370) /16*0,
     1    5909,  -1086,  -357,   203,  -193,    -7,   -46,     8,   -20,
     +                                    2,     0,     0,     0,   3*0,
     2            1065,    34,   -77,    56,    86,   -14,   -15,    14,
     +                                    1,     0,     0,     0,   4*0,
     3                    480,  -201,   -32,     4,   -11,     7,     5,
     +                                    2,     0,     0,     0,   5*0,
     4                           -65,  -125,   -32,     0,   -13,    -3,
     +                                    6,     0,     0,     0,   6*0,
     5                                   11,    32,    28,     5,    -2,
     +                                   -4,     0,     0,     0,   7*0,
     6                                         -67,   -12,    16,     8,
     +                                    0,     0,     0,     0,   8*0,
     7                                                -22,    -5,    10,
     +                                   -2,     0,     0,     0,   9*0,
     8                                                       -18,    -2,
     +                                    4,     0,     0,     0,  10*0,
     9                                                                2/
      DATA (HY1D(I),I=371,450) /
     +                                    0,     0,     0,     0,  11*0,
     O                                   -6,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          h(n,m) for 1910
      DATA (HY1D(I),I=451,595) /16*0,
     1    5898,  -1128,  -389,   211,  -172,    -5,   -47,     8,   -20,
     +                                    2,     0,     0,     0,   3*0,
     2            1000,    62,   -90,    57,    89,   -14,   -15,    14,
     +                                    1,     0,     0,     0,   4*0,
     3                    425,  -189,   -33,     5,   -12,     6,     5,
     +                                    2,     0,     0,     0,   5*0,
     4                           -55,  -126,   -29,     1,   -13,    -3,
     +                                    6,     0,     0,     0,   6*0,
     5                                   21,    28,    28,     5,    -2,
     +                                   -4,     0,     0,     0,   7*0,
     6                                         -65,   -13,    16,     8,
     +                                    0,     0,     0,     0,   8*0,
     7                                                -22,    -5,    10,
     +                                   -2,     0,     0,     0,   9*0,
     8                                                       -18,    -2,
     +                                    4,     0,     0,     0,  10*0,
     9                                                                2/
      DATA (HY1D(I),I=596,675) /
     +                                    0,     0,     0,     0,  11*0,
     O                                   -6,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          h(n,m) for 1915
      DATA (HY1D(I),I=676,820) /16*0,
     1    5875,  -1191,  -421,   218,  -148,    -2,   -48,     8,   -20,
     +                                    2,     0,     0,     0,   3*0,
     2             917,    84,  -109,    58,    93,   -14,   -15,    14,
     +                                    1,     0,     0,     0,   4*0,
     3                    360,  -173,   -34,     8,   -12,     6,     5,
     +                                    2,     0,     0,     0,   5*0,
     4                           -51,  -126,   -26,     2,   -13,    -3,
     +                                    6,     0,     0,     0,   6*0,
     5                                   32,    23,    28,     5,    -2,
     +                                   -4,     0,     0,     0,   7*0,
     6                                         -62,   -15,    16,     8,
     +                                    0,     0,     0,     0,   8*0,
     7                                                -22,    -5,    10,
     +                                   -2,     0,     0,     0,   9*0,
     8                                                       -18,    -2,
     +                                    4,     0,     0,     0,  10*0,
     9                                                                2/
      DATA (HY1D(I),I=821,900) /
     +                                    0,     0,     0,     0,  11*0,
     O                                   -6,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          h(n,m) for 1920
      DATA (HY1D(I),I=901,1045) /16*0,
     1    5845,  -1259,  -445,   220,  -122,     0,   -49,     8,   -20,
     +                                    2,     0,     0,     0,   3*0,
     2             823,   103,  -134,    58,    96,   -14,   -15,    14,
     +                                    1,     0,     0,     0,   4*0,
     3                    293,  -153,   -38,    11,   -13,     6,     5,
     +                                    2,     0,     0,     0,   5*0,
     4                           -57,  -125,   -22,     4,   -14,    -3,
     +                                    6,     0,     0,     0,   6*0,
     5                                   43,    18,    28,     5,    -2,
     +                                   -4,     0,     0,     0,   7*0,
     6                                         -57,   -16,    17,     9,
     +                                    0,     0,     0,     0,   8*0,
     7                                                -22,    -5,    10,
     +                                   -2,     0,     0,     0,   9*0,
     8                                                       -19,    -2,
     +                                    4,     0,     0,     0,  10*0,
     9                                                                2/
      DATA (HY1D(I),I=1046,1125) /
     +                                    0,     0,     0,     0,  11*0,
     O                                   -6,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          h(n,m) for 1925
      DATA (HY1D(I),I=1126,1270) /16*0,
     1    5817,  -1334,  -462,   216,   -96,     3,   -50,     8,   -20,
     +                                    2,     0,     0,     0,   3*0,
     2             728,   119,  -163,    58,    99,   -14,   -15,    14,
     +                                    1,     0,     0,     0,   4*0,
     3                    229,  -130,   -44,    14,   -14,     6,     5,
     +                                    2,     0,     0,     0,   5*0,
     4                           -70,  -122,   -18,     5,   -14,    -3,
     +                                    6,     0,     0,     0,   6*0,
     5                                   51,    13,    29,     5,    -2,
     +                                   -4,     0,     0,     0,   7*0,
     6                                         -52,   -17,    17,     9,
     +                                    0,     0,     0,     0,   8*0,
     7                                                -21,    -5,    10,
     +                                   -2,     0,     0,     0,   9*0,
     8                                                       -19,    -2,
     +                                    4,     0,     0,     0,  10*0,
     9                                                                2/
      DATA (HY1D(I),I=1271,1350) /
     +                                    0,     0,     0,     0,  11*0,
     O                                   -6,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          h(n,m) for 1930
      DATA (HY1D(I),I=1351,1495) /16*0,
     1    5808,  -1424,  -480,   205,   -72,     4,   -51,     8,   -20,
     +                                    2,     0,     0,     0,   3*0,
     2             644,   133,  -195,    60,   102,   -15,   -15,    14,
     +                                    1,     0,     0,     0,   4*0,
     3                    166,  -109,   -53,    19,   -14,     5,     5,
     +                                    2,     0,     0,     0,   5*0,
     4                           -90,  -118,   -16,     6,   -14,    -3,
     +                                    6,     0,     0,     0,   6*0,
     5                                   58,     8,    29,     5,    -2,
     +                                   -4,     0,     0,     0,   7*0,
     6                                         -46,   -18,    18,     9,
     +                                    0,     0,     0,     0,   8*0,
     7                                                -20,    -5,    10,
     +                                   -2,     0,     0,     0,   9*0,
     8                                                       -19,    -2,
     +                                    4,     0,     0,     0,  10*0,
     9                                                                2/
      DATA (HY1D(I),I=1496,1575) /
     +                                    0,     0,     0,     0,  11*0,
     O                                   -6,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          h(n,m) for 1935
      DATA (HY1D(I),I=1576,1720) /16*0,
     1    5812,  -1520,  -494,   188,   -51,     4,   -52,     8,   -20,
     +                                    2,     0,     0,     0,   3*0,
     2             586,   146,  -226,    64,   104,   -17,   -15,    15,
     +                                    1,     0,     0,     0,   4*0,
     3                    101,   -90,   -64,    25,   -14,     5,     5,
     +                                    2,     0,     0,     0,   5*0,
     4                          -114,  -115,   -15,     7,   -15,    -3,
     +                                    6,     0,     0,     0,   6*0,
     5                                   64,     4,    29,     5,    -3,
     +                                   -4,     0,     0,     0,   7*0,
     6                                         -40,   -19,    18,     9,
     +                                    0,     0,     0,     0,   8*0,
     7                                                -19,    -5,    11,
     +                                   -1,     0,     0,     0,   9*0,
     8                                                       -19,    -2,
     +                                    4,     0,     0,     0,  10*0,
     9                                                                2/
      DATA (HY1D(I),I=1721,1800) /
     +                                    0,     0,     0,     0,  11*0,
     O                                   -6,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          h(n,m) for 1940
      DATA (HY1D(I),I=1801,1945) /16*0,
     1    5821,  -1614,  -499,   169,   -33,     4,   -52,     8,   -21,
     +                                    2,     0,     0,     0,   3*0,
     2             528,   163,  -252,    71,   105,   -18,   -14,    15,
     +                                    1,     0,     0,     0,   4*0,
     3                     43,   -72,   -75,    33,   -14,     5,     5,
     +                                    2,     0,     0,     0,   5*0,
     4                          -141,  -113,   -15,     7,   -15,    -3,
     +                                    6,     0,     0,     0,   6*0,
     5                                   69,     0,    29,     5,    -3,
     +                                   -4,     0,     0,     0,   7*0,
     6                                         -33,   -20,    19,     9,
     +                                    0,     0,     0,     0,   8*0,
     7                                                -19,    -5,    11,
     +                                   -1,     0,     0,     0,   9*0,
     8                                                       -19,    -2,
     +                                    4,     0,     0,     0,  10*0,
     9                                                                2/
      DATA (HY1D(I),I=1946,2025) /
     +                                    0,     0,     0,     0,  11*0,
     O                                   -6,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          h(n,m) for 1945
      DATA (HY1D(I),I=2026,2170) /16*0,
     1    5810,  -1702,  -499,   144,   -12,     6,   -45,    12,   -27,
     +                                    5,     0,     0,     0,   3*0,
     2             477,   186,  -276,    95,   100,   -18,   -21,    17,
     +                                    1,     0,     0,     0,   4*0,
     3                    -11,   -55,   -67,    16,     2,   -12,    29,
     +                                  -20,     0,     0,     0,   5*0,
     4                          -178,  -119,    -9,     6,    -7,    -9,
     +                                   -1,     0,     0,     0,   6*0,
     5                                   82,   -16,    28,     2,     4,
     +                                   -6,     0,     0,     0,   7*0,
     6                                         -39,   -17,    18,     9,
     +                                    6,     0,     0,     0,   8*0,
     7                                                -22,     3,     6,
     +                                   -4,     0,     0,     0,   9*0,
     8                                                       -11,     1,
     +                                   -2,     0,     0,     0,  10*0,
     9                                                                8/
      DATA (HY1D(I),I=2171,2250) /
     +                                    0,     0,     0,     0,  11*0,
     O                                   -2,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          h(n,m) for 1950
      DATA (HY1D(I),I=2251,2395) /16*0,
     1    5815,  -1810,  -476,   136,     3,    -1,   -35,     5,   -24,
     +                                   13,     0,     0,     0,   3*0,
     2             381,   206,  -278,   103,    99,   -17,   -22,    19,
     +                                   -2,     0,     0,     0,   4*0,
     3                    -46,   -37,   -87,    33,     0,     0,    12,
     +                                  -10,     0,     0,     0,   5*0,
     4                          -210,  -122,   -12,    10,   -21,     2,
     +                                    2,     0,     0,     0,   6*0,
     5                                   80,   -12,    36,    -8,     2,
     +                                   -3,     0,     0,     0,   7*0,
     6                                         -30,   -18,    17,     8,
     +                                    6,     0,     0,     0,   8*0,
     7                                                -16,    -4,     8,
     +                                   -3,     0,     0,     0,   9*0,
     8                                                       -17,   -11,
     +                                    6,     0,     0,     0,  10*0,
     9                                                               -7/
      DATA (HY1D(I),I=2396,2475) /
     +                                   11,     0,     0,     0,  11*0,
     O                                    8,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          h(n,m) for 1955
      DATA (HY1D(I),I=2476,2620) /16*0,
     1    5820,  -1898,  -462,   133,    15,    -9,   -50,    10,   -11,
     +                                   -4,     0,     0,     0,   3*0,
     2             291,   216,  -274,   110,    96,   -24,   -15,    12,
     +                                    0,     0,     0,     0,   4*0,
     3                    -83,   -23,   -98,    48,    -4,     5,     7,
     +                                   -8,     0,     0,     0,   5*0,
     4                          -230,  -121,   -16,     8,   -23,     6,
     +                                   -2,     0,     0,     0,   6*0,
     5                                   78,   -12,    28,     3,    -2,
     +                                   -4,     0,     0,     0,   7*0,
     6                                         -24,   -20,    23,    10,
     +                                    1,     0,     0,     0,   8*0,
     7                                                -18,    -4,     7,
     +                                   -3,     0,     0,     0,   9*0,
     8                                                       -13,    -6,
     +                                    7,     0,     0,     0,  10*0,
     9                                                                5/
      DATA (HY1D(I),I=2621,2700) /
     +                                   -1,     0,     0,     0,  11*0,
     O                                   -3,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          h(n,m) for 1960
      DATA (HY1D(I),I=2701,2845) /16*0,
     1    5791,  -1967,  -414,   135,    16,   -10,   -55,    11,   -18,
     +                                    4,     0,     0,     0,   3*0,
     2             206,   224,  -278,   125,    99,   -28,   -14,    12,
     +                                    1,     0,     0,     0,   4*0,
     3                   -130,     3,  -117,    60,    -6,     7,     2,
     +                                    0,     0,     0,     0,   5*0,
     4                          -255,  -114,   -20,     7,   -18,     0,
     +                                    2,     0,     0,     0,   6*0,
     5                                   81,   -11,    23,     4,    -3,
     +                                   -5,     0,     0,     0,   7*0,
     6                                         -17,   -18,    23,     9,
     +                                    1,     0,     0,     0,   8*0,
     7                                                -17,     1,     8,
     +                                   -1,     0,     0,     0,   9*0,
     8                                                       -20,     0,
     +                                    6,     0,     0,     0,  10*0,
     9                                                                5/
      DATA (HY1D(I),I=2846,2925) /
     +                                    0,     0,     0,     0,  11*0,
     O                                   -7,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          h(n,m) for 1965
      DATA (HY1D(I),I=2926,3070) /16*0,
     1    5776,  -2016,  -404,   148,    19,   -11,   -61,     7,   -22,
     +                                    2,     0,     0,     0,   3*0,
     2             114,   240,  -269,   128,   100,   -27,   -12,    15,
     +                                    1,     0,     0,     0,   4*0,
     3                   -165,    13,  -126,    68,    -2,     9,     7,
     +                                    2,     0,     0,     0,   5*0,
     4                          -269,   -97,   -32,     6,   -16,    -4,
     +                                    6,     0,     0,     0,   6*0,
     5                                   81,    -8,    26,     4,    -5,
     +                                   -4,     0,     0,     0,   7*0,
     6                                          -7,   -23,    24,    10,
     +                                    0,     0,     0,     0,   8*0,
     7                                                -12,    -3,    10,
     +                                   -2,     0,     0,     0,   9*0,
     8                                                       -17,    -4,
     +                                    3,     0,     0,     0,  10*0,
     9                                                                1/
      DATA (HY1D(I),I=3071,3150) /
     +                                    0,     0,     0,     0,  11*0,
     O                                   -6,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          h(n,m) for 1970
      DATA (HY1D(I),I=3151,3295) /16*0,
     1    5737,  -2047,  -366,   167,    26,   -12,   -70,     7,   -21,
     +                                    1,     0,     0,     0,   3*0,
     2              25,   251,  -266,   139,   100,   -27,   -15,    16,
     +                                    1,     0,     0,     0,   4*0,
     3                   -196,    26,  -139,    72,    -4,     6,     6,
     +                                    3,     0,     0,     0,   5*0,
     4                          -279,   -91,   -37,     8,   -17,    -4,
     +                                    4,     0,     0,     0,   6*0,
     5                                   83,    -6,    23,     6,    -5,
     +                                   -4,     0,     0,     0,   7*0,
     6                                           1,   -23,    21,    10,
     +                                    0,     0,     0,     0,   8*0,
     7                                                -11,    -6,    11,
     +                                   -1,     0,     0,     0,   9*0,
     8                                                       -16,    -2,
     +                                    3,     0,     0,     0,  10*0,
     9                                                                1/
      DATA (HY1D(I),I=3296,3375) /
     +                                    1,     0,     0,     0,  11*0,
     O                                   -4,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          h(n,m) for 1975
      DATA (HY1D(I),I=3376,3520) /16*0,
     1    5675,  -2067,  -333,   191,    31,   -13,   -77,     6,   -21,
     +                                    1,     0,     0,     0,   3*0,
     2             -68,   262,  -265,   148,    99,   -26,   -16,    16,
     +                                    1,     0,     0,     0,   4*0,
     3                   -223,    39,  -152,    75,    -5,     4,     7,
     +                                    3,     0,     0,     0,   5*0,
     4                          -288,   -83,   -41,    10,   -19,    -4,
     +                                    4,     0,     0,     0,   6*0,
     5                                   88,    -4,    22,     6,    -5,
     +                                   -4,     0,     0,     0,   7*0,
     6                                          11,   -23,    18,    10,
     +                                   -1,     0,     0,     0,   8*0,
     7                                                -12,   -10,    11,
     +                                   -1,     0,     0,     0,   9*0,
     8                                                       -17,    -3,
     +                                    3,     0,     0,     0,  10*0,
     9                                                                1/
      DATA (HY1D(I),I=3521,3600) /
     +                                    1,     0,     0,     0,  11*0,
     O                                   -5,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          h(n,m) for 1980
      DATA (HY1D(I),I=3601,3745) /16*0,
     1    5604,  -2129,  -336,   212,    46,   -15,   -82,     7,   -21,
     +                                    1,     0,     0,     0,   3*0,
     2            -200,   271,  -257,   150,    93,   -27,   -18,    16,
     +                                    0,     0,     0,     0,   4*0,
     3                   -252,    53,  -151,    71,    -5,     4,     9,
     +                                    3,     0,     0,     0,   5*0,
     4                          -297,   -78,   -43,    16,   -22,    -5,
     +                                    6,     0,     0,     0,   6*0,
     5                                   92,    -2,    18,     9,    -6,
     +                                   -4,     0,     0,     0,   7*0,
     6                                          17,   -23,    16,     9,
     +                                    0,     0,     0,     0,   8*0,
     7                                                -10,   -13,    10,
     +                                   -1,     0,     0,     0,   9*0,
     8                                                       -15,    -6,
     +                                    4,     0,     0,     0,  10*0,
     9                                                                2/
      DATA (HY1D(I),I=3746,3825) /
     +                                    0,     0,     0,     0,  11*0,
     O                                   -6,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          h(n,m) for 1985
      DATA (HY1D(I),I=3826,3970) /16*0,
     1    5500,  -2197,  -310,   232,    47,   -16,   -83,     8,   -21,
     +                                    1,     0,     0,     0,   3*0,
     2            -306,   284,  -249,   150,    88,   -27,   -19,    15,
     +                                    0,     0,     0,     0,   4*0,
     3                   -297,    69,  -154,    69,    -2,     5,     9,
     +                                    3,     0,     0,     0,   5*0,
     4                          -297,   -75,   -48,    20,   -23,    -6,
     +                                    6,     0,     0,     0,   6*0,
     5                                   95,    -1,    17,    11,    -6,
     +                                   -4,     0,     0,     0,   7*0,
     6                                          21,   -23,    14,     9,
     +                                    0,     0,     0,     0,   8*0,
     7                                                 -7,   -15,     9,
     +                                   -1,     0,     0,     0,   9*0,
     8                                                       -11,    -7,
     +                                    4,     0,     0,     0,  10*0,
     9                                                                2/
      DATA (HY1D(I),I=3971,4050) /
     +                                    0,     0,     0,     0,  11*0,
     O                                   -6,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          h(n,m) for 1990
      DATA (HY1D(I),I=4051,4195) /16*0,
     1    5406,  -2279,  -284,   247,    46,   -16,   -80,    10,   -20,
     +                                    2,     0,     0,     0,   3*0,
     2            -373,   293,  -240,   154,    82,   -26,   -19,    15,
     +                                    1,     0,     0,     0,   4*0,
     3                   -352,    84,  -153,    69,     0,     6,    11,
     +                                    3,     0,     0,     0,   5*0,
     4                          -299,   -69,   -52,    21,   -22,    -7,
     +                                    6,     0,     0,     0,   6*0,
     5                                   97,     1,    17,    12,    -7,
     +                                   -4,     0,     0,     0,   7*0,
     6                                          24,   -23,    12,     9,
     +                                    0,     0,     0,     0,   8*0,
     7                                                 -4,   -16,     8,
     +                                   -2,     0,     0,     0,   9*0,
     8                                                       -10,    -7,
     +                                    3,     0,     0,     0,  10*0,
     9                                                                2/
      DATA (HY1D(I),I=4196,4275) /
     +                                   -1,     0,     0,     0,  11*0,
     O                                   -6,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          h(n,m) for 1995
      DATA (HY1D(I),I=4276,4420) /16*0,
     1    5306,  -2366,  -262,   262,    46,   -17,   -69,    11,   -20,
     +                                    1,     0,     0,     0,   3*0,
     2            -413,   302,  -236,   165,    72,   -25,   -21,    15,
     +                                    0,     0,     0,     0,   4*0,
     3                   -427,    97,  -143,    67,     4,     8,    12,
     +                                    4,     0,     0,     0,   5*0,
     4                          -306,   -55,   -58,    24,   -23,    -6,
     +                                    5,     0,     0,     0,   6*0,
     5                                  107,     1,    17,    15,    -8,
     +                                   -5,     0,     0,     0,   7*0,
     6                                          36,   -24,    11,     8,
     +                                   -1,     0,     0,     0,   8*0,
     7                                                 -6,   -16,     5,
     +                                   -2,     0,     0,     0,   9*0,
     8                                                        -4,    -8,
     +                                    1,     0,     0,     0,  10*0,
     9                                                                3/
      DATA (HY1D(I),I=4421,4500) /
     +                                   -2,     0,     0,     0,  11*0,
     O                                   -7,     0,     0,     0,  12*0,
     1                                           0,     0,     0,  13*0,
     2                                                  0,     0,  14*0,
     3                                                         0,  16*0/
C          h(n,m) for 2000
      DATA (HY1D(I),I=4501,4645) /16*0,
     1  5186.1,-2481.6,-227.6, 272.6,  43.8, -17.4, -64.6,  11.9, -19.7,
     +                                  1.7,   0.1,  -0.4,  -0.9,   3*0,
     2          -458.0, 293.4,-231.9, 171.9,  63.7, -24.2, -21.5,  13.4,
     +                                  0.0,   1.3,   0.3,   0.2,   4*0,
     3                 -491.1, 119.8,-133.1,  65.1,   6.2,   8.5,  12.5,
     +                                  4.0,  -0.9,   2.5,   1.8,   5*0,
     4                        -303.8, -39.3, -61.2,  24.0, -21.5,  -6.2,
     +                                  4.9,  -2.6,  -2.6,  -0.4,   6*0,
     5                                106.3,   0.7,  14.8,  15.5,  -8.4,
     +                                 -5.9,   0.9,   0.7,  -1.0,   7*0,
     6                                        43.8, -25.4,   8.9,   8.4,
     +                                 -1.2,  -0.7,   0.3,  -0.1,   8*0,
     7                                               -5.8, -14.9,   3.8,
     +                                 -2.9,  -2.8,   0.0,   0.7,   9*0,
     8                                                      -2.1,  -8.2,
     +                                  0.2,  -0.9,   0.0,   0.3,  10*0,
     9                                                              4.8/
      DATA (HY1D(I),I=4646,4725) /
     +                                 -2.2,  -1.2,   0.3,   0.6,  11*0,
     O                                 -7.4,  -1.9,  -0.9,   0.3,  12*0,
     1                                        -0.9,  -0.4,  -0.2,  13*0,
     2                                                0.8,  -0.5,  14*0,
     3                                                      -0.9,  16*0/
C          h(n,m) for 2005
      DATA (HY1D(I),I=4726,4870) /16*0,
     1 5077.99,-2594.50,-198.86,282.07, 42.72,-20.33,-61.14, 11.20,
     +                        -20.11, 2.19,  0.26, -0.55, -0.76,   3*0,
     2         -515.43,269.72,-225.23,180.25, 54.75,-22.57,-20.88,12.69,
     +                                 0.10,  1.44,  0.23,  0.33,   4*0,
     3                 -524.72,145.15,-123.45, 63.63,  6.82, 9.83,12.67,
     +                                 4.46, -0.77,  2.38,  1.72,   5*0,
     4                        -305.36,-19.57,-63.53, 25.35,-19.71,-6.72,
     +                                 4.76, -2.27, -2.63, -0.54,   6*0,
     5                               103.85,  0.24, 10.93, 16.22, -8.16,
     +                                -6.58,  0.90,  0.61, -1.07,   7*0,
     6                                       50.94,-26.32,  7.61,  8.10,
     +                                -1.01, -0.58,  0.40, -0.04,   8*0,
     7                                              -4.64,-12.76,  2.92,
     +                                -3.47, -2.69,  0.01,  0.63,   9*0,
     8                                                     -0.06, -7.73,
     +                                -0.86, -1.08,  0.02,  0.21,  10*0,
     9                                                             6.01/
      DATA (HY1D(I),I=4871,4950) /
     +                                -2.31, -1.58,  0.28,  0.53,  11*0,
     O                                -7.93, -1.90, -0.87,  0.38,  12*0,
     1                                       -1.39, -0.34, -0.22,  13*0,
     2                                               0.88, -0.57,  14*0,
     3                                                     -0.82,  16*0/
C          h(n,m) for 2010
      DATA (HY1D(I),I=4951,5095) /16*0,
     1 4944.26, -2708.54,-160.4,286.48,44.58,-20.9, -57.8,10.84, -20.54,
     +                                  2.73, 0.13, -0.87, -0.87,   3*0,
     2    -575.73, 251.75, -211.03, 189.01, 44.18, -21.2, -20.03, 11.51,
     +                                     -0.1, 1.67, 0.27, 0.3,   4*0,
     3              -537.03, 164.46, -118.06, 61.54, 6.54, 11.83, 12.75,
     +                                   4.71, -0.66, 2.13, 1.66,   5*0,
     4                     -309.72, -0.01, -66.26, 24.96, -17.41, -7.14,
     +                                 4.44, -1.76, -2.49, -0.59,   6*0,
     5                                 101.04, 3.02, 7.03, 16.71, -7.42,
     +                                  -7.22, 0.85, 0.49, -1.14,   7*0,
     6                                         55.4, -27.61, 6.96, 7.97,
     +                                 -0.96, -0.39, 0.59, -0.07,   8*0,
     7                                              -3.28, -10.74, 2.14,
     +                                   -3.95, -2.51, 0.0, 0.54,   9*0,
     8                                                      1.64, -6.08,
     +                                   -1.99, -1.27, 0.13, 0.1,  10*0,
     9                                                             7.01/
      DATA (HY1D(I),I=5096,5175) /
     +                                  -1.97, -2.11, 0.27, 0.49,  11*0,
     O                                 -8.31, -1.94, -0.86, 0.44,  12*0,
     1                                       -1.86, -0.23, -0.25,  13*0,
     2                                               0.87, -0.53,  14*0,
     3                                                     -0.79,  16*0/
C          h(n,m) for 2015
      DATA (HY1D(I),I=5176,5320) /16*0,
     1  4797.1, -2845.6, -115.3, 283.3, 47.3, -20.8, -54.1, 10.1, -21.6,
     +                                     3.2, -0.1, -1.1, -0.9,   3*0,
     2             -641.9, 244.9, -188.7, 197, 33.2, -19.5, -18.3, 10.8,
     +                                       -0.4, 2.0, 0.4, 0.4,   4*0,
     3                     -538.4, 180.9, -119.3, 58.9, 5.7, 13.3, 11.8,
     +                                       4.6, -0.7, 1.9, 1.6,   5*0,
     4                           -329.5, 16.0, -66.7, 24.4, -14.6, -6.8,
     +                                     4.4, -1.1, -2.2, -0.5,   6*0,
     5                                      100.2, 7.3, 3.4, 16.2, -6.9,
     +                                      -7.9, 0.8, 0.3, -1.2,   7*0,
     6                                            62.6, -27.4, 5.7, 7.8,
     +                                     -0.6, -0.2, 0.7, -0.1,   8*0,
     7                                                  -2.2, -9.1, 1.0,
     +                                     -4.2, -2.2, -0.1, 0.4,   9*0,
     8                                                        2.1, -4.0,
     +                                     -2.8, -1.4, 0.3, -0.1,  10*0,
     9                                                              8.4/
      DATA (HY1D(I),I=5321,5400) /
     +                                      -1.2, -2.5, 0.2, 0.4,  11*0,
     O                                     -8.7, -2.0, -0.9, 0.5,  12*0,
     1                                          -2.4, -0.1, -0.3,  13*0,
     2                                                 0.7, -0.4,  14*0,
     3                                                      -0.8,  16*0/
C          Secular variation rates are nominally okay through 2020
      DATA (GT1D(I),I=1,145) /0,
     O    10.3,   -8.7,   3.4,  -0.7,  -0.2,  -0.3,   0.3,   0.2,   0.0,
     +                                  0.0,   0.0,   0.0,   0.0,   2*0,
     1    18.1,   -3.3,  -5.5,   0.2,   0.5,  -0.1,  -0.2,   0.0,   0.0,
     +                                  0.0,   0.0,   0.0,   0.0,   3*0,
     2             2.1,  -0.7,  -9.1,  -1.3,  -0.7,  -0.5,  -0.6,   0.0,
     +                                  0.0,   0.0,   0.0,   0.0,   4*0,
     3                  -10.1,   4.1,  -0.1,   2.1,   1.3,   0.5,   0.0,
     +                                  0.0,   0.0,   0.0,   0.0,   5*0,
     4                          -4.3,   1.4,  -1.2,   0.1,  -0.2,   0.0,
     +                                  0.0,   0.0,   0.0,   0.0,   6*0,
     5                                  3.9,   0.3,  -0.6,   0.4,   0.0,
     +                                  0.0,   0.0,   0.0,   0.0,   7*0,
     6                                         1.6,  -0.8,   0.1,   0.0,
     +                                  0.0,   0.0,   0.0,   0.0,   8*0,
     7                                                0.2,  -0.4,   0.0,
     +                                  0.0,   0.0,   0.0,   0.0,   9*0,
     8                                                       0.3,   0.0,
     +                                  0.0,   0.0,   0.0,   0.0,  10*0,
     9                                                              0.0/
      DATA (GT1D(I),I=146,225) /
     +                                  0.0,   0.0,   0.0,   0.0,  11*0,
     O                                  0.0,   0.0,   0.0,   0.0,  12*0,
     1                                         0.0,   0.0,   0.0,  13*0,
     2                                                0.0,   0.0,  14*0,
     3                                                       0.0,  16*0/
      DATA (HT1D(I),I=1,145) /16*0,
     1 -26.6,   -27.4,   8.2,   -1.3,   0.6,   0.0,   0.8,  -0.3,   0.0,
     +                                  0.0,   0.0,   0.0,   0.0,   3*0,
     2           -14.1,  -0.4,   5.3,   1.7,  -2.1,   0.4,   0.3,   0.0,
     +                                  0.0,   0.0,   0.0,   0.0,   4*0,
     3                    1.8,   2.9,  -1.2,  -0.7,  -0.2,   0.1,   0.0,
     +                                  0.0,   0.0,   0.0,   0.0,   5*0,
     4                          -5.2,   3.4,   0.2,  -0.3,   0.5,   0.0,
     +                                  0.0,   0.0,   0.0,   0.0,   6*0,
     5                                  0.0,   0.9,  -0.6,  -0.2,   0.0,
     +                                  0.0,   0.0,   0.0,   0.0,   7*0,
     6                                         1.0,   0.1,  -0.3,   0.0,
     +                                  0.0,   0.0,   0.0,   0.0,   8*0,
     7                                               -0.2,   0.3,   0.0,
     +                                  0.0,   0.0,   0.0,   0.0,   9*0,
     8                                                       0.0,   0.0,
     +                                  0.0,   0.0,   0.0,   0.0,  10*0,
     9                                                              0.0/
      DATA (HT1D(I),I=146,225) /
     +                                  0.0,   0.0,   0.0,   0.0,  11*0,
     O                                  0.0,   0.0,   0.0,   0.0,  12*0,
     1                                         0.0,   0.0,   0.0,  13*0,
     2                                                0.0,   0.0,  14*0,
     3                                                       0.0,  16*0/

C          Do not need to load new coefficients if date has not changed
      ICHG = 0
      IF (DATE .EQ. DATEL) GO TO 300
      DATEL = DATE
      ICHG = 1
 
C          Trap out of range date:
      IF (DATE .LT. EPOCH(1)) GO TO 9100
      IF (DATE .GT. EPOCH(NEPO)+5.) WRITE(0,9200) DATE, EPOCH(NEPO) + 5.
 
      DO 100 I=1,NEPO
      IF (DATE .LT. EPOCH(I)) GO TO 110
      IY = I
  100 CONTINUE
  110 CONTINUE
 
      NMAX  = NMXE(IY)
      TIME  = DATE
      T     = TIME-EPOCH(IY)
      TO5   = T/5.
      IY1   = IY + 1
      GB(1) = 0.0
      GV(1) = 0.0
      I  = 2
      F0 = -1.0D-5
      DO 200 N=1,NMAX
      F0 = F0 * REAL(N)/2.
      F  = F0 / SQRT(2.0)
      NN = N+1
      MM = 1
      IF (IY .LT. NEPO) GB(I) = (GYR(NN,MM,IY) +                             ! interpolate (m=0 terms)
     +                          (GYR(NN,MM,IY1)-GYR(NN,MM,IY))*TO5) * F0
      IF (IY .EQ. NEPO) GB(I) = (GYR(NN,MM,IY) + GT(NN,MM)    *T  ) * F0     ! extrapolate (m=0 terms)
      GV(I) = GB(I) / REAL(NN)
      I = I+1
      DO 200 M=1,N
      F  = F / SQRT( REAL(N-M+1) / REAL(N+M) )
      NN = N+1
      MM = M+1
      I1 = I+1
      IF (IY .LT. NEPO) THEN                                                ! interpolate (m>0 terms)
	GB(I)  = (GYR(NN,MM,IY) +
     +           (GYR(NN,MM,IY1)-GYR(NN,MM,IY))*TO5) * F
	GB(I1) = (HYR(NN,MM,IY) +
     +           (HYR(NN,MM,IY1)-HYR(NN,MM,IY))*TO5) * F
      ELSE                                                                  ! extrapolate (m>0 terms)
	GB(I)  = (GYR(NN,MM,IY) +GT (NN,MM)    *T  ) * F
	GB(I1) = (HYR(NN,MM,IY) +HT (NN,MM)    *T  ) * F
      ENDIF
      RNN = REAL(NN)
      GV(I)  = GB(I)  / RNN
      GV(I1) = GB(I1) / RNN
  200 I = I+2
 
  300 CONTINUE

      RETURN
 
C          Error trap diagnostics:
 9100 WRITE (0,'(''COFRM:  DATE'',F9.3,'' preceeds earliest available ('
     +',F6.1,'')'')') DATE, EPOCH(1)
      CALL EXIT (1)
 9200 FORMAT('COFRM:  DATE',F9.3,' is after the last recommended for ext
     +rapolation (',F6.1,')')
      END
 
      SUBROUTINE DYPOL (COLAT,ELON,VP)
C          Computes parameters for dipole component of geomagnetic field.
C          COFRM must be called before calling DYPOL!
C          940504 A. D. Richmond
C
C          INPUT from COFRM through COMMON /MAGCOF/ NMAX,GB(255),GV(225),ICHG
C            NMAX = Maximum order of spherical harmonic coefficients used
C            GB   = Coefficients for magnetic field calculation
C            GV   = Coefficients for magnetic potential calculation
C            ICHG = Flag indicating when GB,GV have been changed
C
C          RETURNS:
C            COLAT = Geocentric colatitude of geomagnetic dipole north pole
C                    (deg)
C            ELON  = East longitude of geomagnetic dipole north pole (deg)
C            VP    = Magnitude, in T.m, of dipole component of magnetic
C                    potential at geomagnetic pole and geocentric radius
C                    of 6371.2 km
 
      PARAMETER (RTOD = 57.2957795130823, RE = 6371.2)
      COMMON /MAGCOF/ NMAX,GB(255),GV(225),ICHG
 
C          Compute geographic colatitude and longitude of the north pole of
C          earth centered dipole
      GPL   = SQRT (GB(2)**2 + GB(3)**2 + GB(4)**2)
      CTP   = GB(2) / GPL
      STP   = SQRT (1. - CTP*CTP)
      COLAT = ACOS (CTP) * RTOD
      ELON  = ATAN2 (GB(4),GB(3)) * RTOD
 
C          Compute magnitude of magnetic potential at pole, radius Re.
      VP = .2*GPL*RE
C          .2 = 2*(10**-4 T/gauss)*(1000 m/km) (2 comes through F0 in COFRM).
 
      RETURN
      END
 
      SUBROUTINE FELDG (IENTY,GLAT,GLON,ALT, BNRTH,BEAST,BDOWN,BABS)
C          Compute the DGRF/IGRF field components at the point GLAT,GLON,ALT.
C          COFRM must be called to establish coefficients for correct date
C          prior to calling FELDG.
C
C          IENTY is an input flag controlling the meaning and direction of the
C                remaining formal arguments:
C          IENTY = 1
C            INPUTS:
C              GLAT = Latitude of point (deg)
C              GLON = Longitude (east=+) of point (deg)
C              ALT  = Ht of point (km)
C            RETURNS:
C              BNRTH  north component of field vector (Gauss)
C              BEAST  east component of field vector  (Gauss)
C              BDOWN  downward component of field vector (Gauss)
C              BABS   magnitude of field vector (Gauss)
C
C          IENTY = 2
C            INPUTS:
C              GLAT = X coordinate (in units of earth radii 6371.2 km )
C              GLON = Y coordinate (in units of earth radii 6371.2 km )
C              ALT  = Z coordinate (in units of earth radii 6371.2 km )
C            RETURNS:
C              BNRTH = X component of field vector (Gauss)
C              BEAST = Y component of field vector (Gauss)
C              BDOWN = Z component of field vector (Gauss)
C              BABS  = Magnitude of field vector (Gauss)
C          IENTY = 3
C            INPUTS:
C              GLAT = X coordinate (in units of earth radii 6371.2 km )
C              GLON = Y coordinate (in units of earth radii 6371.2 km )
C              ALT  = Z coordinate (in units of earth radii 6371.2 km )
C            RETURNS:
C              BNRTH = Dummy variable
C              BEAST = Dummy variable
C              BDOWN = Dummy variable
C              BABS  = Magnetic potential (T.m)
C
C          INPUT from COFRM through COMMON /MAGCOF/ NMAX,GB(255),GV(225),ICHG
C            NMAX = Maximum order of spherical harmonic coefficients used
C            GB   = Coefficients for magnetic field calculation
C            GV   = Coefficients for magnetic potential calculation
C            ICHG = Flag indicating when GB,GV have been changed
C
C          HISTORY:
C          Apr 1983: written by Vincent B. Wickwar (Utah State Univ.).
C
C          May 1994 (A.D. Richmond): Added magnetic potential calculation
C
C          Oct 1995 (Barnes): Added ICHG
 
      PARAMETER (DTOR = 0.01745329251994330, RE = 6371.2)
      COMMON /MAGCOF/ NMAX,GB(255),GV(225),ICHG
      DIMENSION G(255), H(255), XI(3)
      SAVE IENTYP, G
      DATA IENTYP/-10000/
 
      IF (IENTY .EQ. 1) THEN
	IS   = 1
	RLAT = GLAT * DTOR
	CT   = SIN (RLAT)
	ST   = COS (RLAT)
	RLON = GLON * DTOR
	CP   = COS (RLON)
	SP   = SIN (RLON)
        CALL GD2CART (GLAT,GLON,ALT,XXX,YYY,ZZZ)
        XXX = XXX/RE
        YYY = YYY/RE
        ZZZ = ZZZ/RE
      ELSE
        IS   = 2
        XXX  = GLAT
        YYY  = GLON
        ZZZ  = ALT
      ENDIF
      RQ    = 1./(XXX**2+YYY**2+ZZZ**2)
      XI(1) = XXX*RQ
      XI(2) = YYY*RQ
      XI(3) = ZZZ*RQ
      IHMAX = NMAX*NMAX+1
      LAST  = IHMAX+NMAX+NMAX
      IMAX  = NMAX+NMAX-1
 
      IF (IENTY .NE. IENTYP .OR. ICHG .EQ. 1) THEN
        IENTYP = IENTY
	ICHG = 0
        IF (IENTY .NE. 3) THEN
	  DO 10 I=1,LAST
   10     G(I) = GB(I)
        ELSE
	  DO 20 I=1,LAST
   20     G(I) = GV(I)
        ENDIF
      ENDIF
 
      DO 30 I=IHMAX,LAST
   30 H(I) = G(I)

      MK = 3
      IF (IMAX .EQ. 1) MK=1

      DO 100 K=1,MK,2
      I  = IMAX
      IH = IHMAX

   60 IL = IH-I
      F = 2./FLOAT(I-K+2)
      X = XI(1)*F
      Y = XI(2)*F
      Z = XI(3)*(F+F)

      I = I-2
      IF (I .LT. 1) GO TO 90
      IF (I .EQ. 1) GO TO 80

      DO 70 M=3,I,2
      IHM = IH+M
      ILM = IL+M
      H(ILM+1) = G(ILM+1)+ Z*H(IHM+1) + X*(H(IHM+3)-H(IHM-1))
     +                                        -Y*(H(IHM+2)+H(IHM-2))
   70 H(ILM)   = G(ILM)  + Z*H(IHM)   + X*(H(IHM+2)-H(IHM-2))
     +                                        +Y*(H(IHM+3)+H(IHM-1))

   80 H(IL+2) = G(IL+2) + Z*H(IH+2) + X*H(IH+4) - Y*(H(IH+3)+H(IH))
      H(IL+1) = G(IL+1) + Z*H(IH+1) + Y*H(IH+4) + X*(H(IH+3)-H(IH))

   90 H(IL)   = G(IL)   + Z*H(IH)   + 2.*(X*H(IH+1)+Y*H(IH+2))
      IH = IL
      IF (I .GE. K) GO TO 60
  100 CONTINUE
 
      S = .5*H(1)+2.*(H(2)*XI(3)+H(3)*XI(1)+H(4)*XI(2))
      T = (RQ+RQ)*SQRT(RQ)
      BXXX = T*(H(3)-S*XXX)
      BYYY = T*(H(4)-S*YYY)
      BZZZ = T*(H(2)-S*ZZZ)
      BABS = SQRT(BXXX**2+BYYY**2+BZZZ**2)
      IF (IS .EQ. 1) THEN            ! (convert back to geodetic)
        BEAST = BYYY*CP-BXXX*SP
        BRHO  = BYYY*SP+BXXX*CP
        BNRTH =  BZZZ*ST-BRHO*CT
        BDOWN = -BZZZ*CT-BRHO*ST
      ELSEIF (IS .EQ. 2) THEN        ! (leave in earth centered cartesian)
        BNRTH = BXXX
        BEAST = BYYY
        BDOWN = BZZZ
      ENDIF
 
C          Magnetic potential computation makes use of the fact that the
C          calculation of V is identical to that for r*Br, if coefficients
C          in the latter calculation have been divided by (n+1) (coefficients
C          GV).  Factor .1 converts km to m and gauss to tesla.
      IF (IENTY.EQ.3) BABS = (BXXX*XXX + BYYY*YYY + BZZZ*ZZZ)*RE*.1
 
      RETURN
      END
 
      SUBROUTINE GD2CART (GDLAT,GLON,ALT,X,Y,Z)
C          Convert geodetic to cartesian coordinates by calling CONVRT
C          940503 A. D. Richmond
      PARAMETER (DTOR = 0.01745329251994330)
      CALL CONVRT (1,GDLAT,ALT,RHO,Z)
      ANG = GLON*DTOR
      X = RHO*COS(ANG)
      Y = RHO*SIN(ANG)
      RETURN
      END
 
      SUBROUTINE CONVRT (I,GDLAT,ALT,X1,X2)
C          Convert space point from geodetic to geocentric or vice versa.
C
C          I is an input flag controlling the meaning and direction of the
C            remaining formal arguments:
C
C          I = 1  (convert from geodetic to cylindrical geocentric)
C            INPUTS:
C              GDLAT = Geodetic latitude (deg)
C              ALT   = Altitude above reference ellipsoid (km)
C            RETURNS:
C              X1    = Distance from Earth's rotation axis (km)
C              X2    = Distance above (north of) Earth's equatorial plane (km)
C
C          I = 2  (convert from geodetic to spherical geocentric)
C            INPUTS:
C              GDLAT = Geodetic latitude (deg)
C              ALT   = Altitude above reference ellipsoid (km)
C            RETURNS:
C              X1    = Geocentric latitude (deg)
C              X2    = Geocentric distance (km)
C
C          I = 3  (convert from cylindrical geocentric to geodetic)
C            INPUTS:
C              X1    = Distance from Earth's rotation axis (km)
C              X2    = Distance from Earth's equatorial plane (km)
C            RETURNS:
C              GDLAT = Geodetic latitude (deg)
C              ALT   = Altitude above reference ellipsoid (km)
C
C          I = 4  (convert from spherical geocentric to geodetic)
C            INPUTS:
C              X1    = Geocentric latitude (deg)
C              X2    = Geocentric distance (km)
C            RETURNS:
C              GDLAT = Geodetic latitude (deg)
C              ALT   = Altitude above reference ellipsoid (km)
C
C
C          HISTORY:
C          940503 (A. D. Richmond):  Based on a routine originally written
C          by V. B. Wickwar.
C
C          Mar 2004: (Barnes) Revise spheroid definition to WGS-1984 to conform
C          with IGRF-9 release (EOS Volume 84 Number 46 November 18 2003).
C
C          REFERENCE: ASTRON. J. VOL. 66, p. 15-16, 1961
 
C          E2  = square of eccentricity of ellipse
C          REP = earth's polar      radius (km)
C          REQ = earth's equatorial radius (km)
      PARAMETER (RTOD = 57.2957795130823, DTOR = 0.01745329251994330,
     +           REP  = 6356.752, REQ = 6378.137, E2 = 1.-(REP/REQ)**2,
     +     E4 = E2*E2, E6 = E4*E2, E8 = E4*E4, OME2REQ = (1.-E2)*REQ,
     +     A21 =     (512.*E2 + 128.*E4 + 60.*E6 + 35.*E8)/1024. ,
     +     A22 =     (                        E6 +     E8)/  32. ,
     +     A23 = -3.*(                     4.*E6 +  3.*E8)/ 256. ,
     +     A41 =    -(           64.*E4 + 48.*E6 + 35.*E8)/1024. ,
     +     A42 =     (            4.*E4 +  2.*E6 +     E8)/  16. ,
     +     A43 =                                   15.*E8 / 256. ,
     +     A44 =                                      -E8 /  16. ,
     +     A61 =  3.*(                     4.*E6 +  5.*E8)/1024. ,
     +     A62 = -3.*(                        E6 +     E8)/  32. ,
     +     A63 = 35.*(                     4.*E6 +  3.*E8)/ 768. ,
     +     A81 =                                   -5.*E8 /2048. ,
     +     A82 =                                   64.*E8 /2048. ,
     +     A83 =                                 -252.*E8 /2048. ,
     +     A84 =                                  320.*E8 /2048. )
 
      IF (I .GE. 3) GO TO 300
 
C          Geodetic to geocentric
 
C          Compute RHO,Z
      SINLAT = SIN(GDLAT*DTOR)
      COSLAT = SQRT(1.-SINLAT*SINLAT)
      D      = SQRT(1.-E2*SINLAT*SINLAT)
      Z      = (ALT+OME2REQ/D)*SINLAT
      RHO    = (ALT+REQ/D)*COSLAT
      X1 = RHO
      X2 = Z
      IF (I .EQ. 1) RETURN
 
C          Compute GCLAT,RKM
      RKM   = SQRT(Z*Z + RHO*RHO)
      GCLAT = RTOD*ATAN2(Z,RHO)
      X1 = GCLAT
      X2 = RKM
      RETURN
 
C          Geocentric to geodetic
  300 IF (I .EQ. 3) THEN
         RHO = X1
         Z = X2
         RKM = SQRT(Z*Z+RHO*RHO)
         SCL = Z/RKM
         GCLAT = ASIN(SCL)*RTOD
      ELSEIF (I .EQ. 4) THEN
         GCLAT = X1
         RKM = X2
         SCL = SIN(GCLAT*DTOR)
      ELSE
         RETURN
      ENDIF
 
      RI = REQ/RKM
      A2 = RI*(A21+RI*(A22+RI* A23))
      A4 = RI*(A41+RI*(A42+RI*(A43+RI*A44)))
      A6 = RI*(A61+RI*(A62+RI* A63))
      A8 = RI*(A81+RI*(A82+RI*(A83+RI*A84)))
      CCL = SQRT(1.-SCL*SCL)
      S2CL = 2.*SCL*CCL
      C2CL = 2.*CCL*CCL-1.
      S4CL = 2.*S2CL*C2CL
      C4CL = 2.*C2CL*C2CL-1.
      S8CL = 2.*S4CL*C4CL
      S6CL = S2CL*C4CL+C2CL*S4CL
      DLTCL = S2CL*A2+S4CL*A4+S6CL*A6+S8CL*A8
      GDLAT = DLTCL*RTOD+GCLAT
      SGL = SIN(GDLAT*DTOR)
      ALT = RKM*COS(DLTCL)-REQ*SQRT(1.-E2*SGL*SGL)
      RETURN
      END
