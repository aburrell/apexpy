C  FILE NAME: magfld.f

      SUBROUTINE COFRM (DATE, FILENAME)
C          Define the International Geomagnetic Reference Field (IGRF) as a
C          scalar potential field using a truncated series expansion with
C          Schmidt semi-normalized associated Legendre functions of degree n and
C          order m.  The polynomial coefficients are a function of time and are
C          interpolated between five year epochs or extrapolated at a constant
C          rate after the last epoch.
C
C          INPUTS:
C            DATE = yyyy.fraction (UT)
C            FILENAME = filename for IGRF coefficient file
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
C
C          May 2020 (Achim Morschhauser): Update with routine to read
C          IGRF coefficients file directly.
C
      use igrf
C
c      implicit none
C
      REAL(4) DATE,DATEL
      REAL(4) F,F0
      REAL(4), ALLOCATABLE :: EPOCH(:),NMXE(:)
      COMMON /MAGCOF/ NMAX,GB(255),GV(225),ICHG
      DATA ICHG /-99999/
 
C          NEPO = Number of epochs
C          NGH  = Single dimensioned array size of 2D version (GYR or HYR)
C          NGHT = Single dimensioned array size of 2D version (GT  or HT)
!!      PARAMETER (NEPO = 24, NGH = 225*NEPO, NGHT = 225)
      INTEGER NEPO,NGHT,NGH
      REAL*8, ALLOCATABLE :: GYR(:,:,:),HYR(:,:,:) 
      REAL*8, ALLOCATABLE :: GT(:,:),HT(:,:) 
     
      CHARACTER(LEN=1000) :: FILENAME

      SAVE DATEL, GYR, HYR, GT, HT, NEPO, EPOCH, NGHT, NMXE
      DATA DATEL /-999./

C          Do not need to load new coefficients if date has not changed
      ICHG = 0
      IF (DATE .EQ. DATEL) GO TO 300
      DATEL = DATE
      ICHG = 1

c          Load coefficients
      if (.not. allocated(GYR)) then
        call read_igrf(FILENAME,GYR,HYR,GT,HT,NEPO,NGHT,EPOCH,NMXE)
      endif
      NGH=NGHT*NEPO
 
C          Trap out of range date:
      IF (DATE .LT. EPOCH(1)) GO TO 9100
      IF (DATE .GT. EPOCH(NEPO)+5.) WRITE(0,9200) DATE, EPOCH(NEPO) + 5.
 
      DO 100 I=1,NEPO
      IF (DATE .LT. EPOCH(I)) GO TO 110
      IY = I
  100 CONTINUE
  110 CONTINUE
 
      NGH=NGHT*NEPO
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
