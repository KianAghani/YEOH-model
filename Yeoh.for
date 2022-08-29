C *******************************************************************
C  VUMAT for 3D Generalized Yeoh Strain-Energy Function 
C  Compatible with Abaqus 2020 & R2018
C  ---------------------------------------------------
C  This work is based on:
C  Hohenberger et. al (2019), "A CONSTITUTIVE MODEL FOR BOTH LOW AND 
C  HIGH STRAIN NONLINEARITIES IN HIGHLY FILLED ELASTOMERS AND  
C  IMPLEMENTATION WITH USER-DEFINED MATERIAL SUBROUTINES IN ABAQUS",
C  Rubber Chem. & Tech., 92 (4), 653 
C *******************************************************************
C
C  Strain-energy function:
C
C  W = A1*(I1-3)^m + A2*(I1-3)^p + A3*(I1-3)^q + (1/D1)*(J-1)^2
C
C **********************************************************************
C
       SUBROUTINE VUMAT(
     1     nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2     stepTime, totalTime, dt, cmname, coordMp, charLength,
     3     props, density, strainInc, relSpinInc,
     4     tempOld, stretchOld, defgradOld, fieldOld,
     5     stressOld, stateOld, enerInternOld, enerInelasOld,
     6     tempNew, stretchNew, defgradNew, fieldNew,
     7     stressNew, stateNew, enerInternNew, enerInelasNew )
C
      INCLUDE 'vaba_param.inc'
C
      DIMENSION props(nprops), density(nblock), coordMp(nblock,*),
     1          charLength(nblock), strainInc(nblock,ndir+nshr),
     2          relSpinInc(nblock,nshr), tempOld(nblock),
     3          stretchOld(nblock,ndir+nshr),
     4          defgradOld(nblock,ndir+nshr),
     5          fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6          stateOld(nblock,nstatev), enerInternOld(nblock),
     7          enerInelasOld(nblock), tempNew(nblock),
     8          stretchNew(nblock,ndir+nshr),
     9          defgradNew(nblock,ndir+nshr),
     1          fieldNew(nblock,nfieldv),
     2          stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     3          enerInternNew(nblock), enerInelasNew(nblock)
C
      CHARACTER*80 cmname
C
C     PARAMETERS
C     ----------
      REAL*8    oneThrd, half, twoThrd, one, two, three, toler
      PARAMETER(oneThrd=1.d0/3.d0, half=0.5d0, twoThrd=2.d0/3.d0,
     1          one=1.d0, two=2.d0, three=3.d0, toler=10.d0**-12.d0)
C
C     LOCAL VARIABLES
C     ---------------
      REAL*8 k0      , twoG   , Lambda , trace  , EP     , QU      ,
     1       A1      , A2     , A3     , E      , J      , Scale   , 
     2       Bxx     , Byy    , Bzz    , Bxy    , Bxz    , Byz     ,
     3       BbarXX  , BbarYY , BbarZZ , BbarXY , BbarXZ , BbarYZ  ,
     4       dBbarXX , dBbarYY, dBbarZZ, dBbarXY, dBbarXZ, dBbarYZ ,
     5       dWdI1   , dWdJ   , I1     , p0        
C
      A1 = props(1)
      A2 = props(2)
      A3 = props(3)
      E  = props(4)
      EP = props(5)
      QU = props(6)
      D1 = props(7)
C
      EG = two * A1
      Ebulk = two / D1
C
      twoG = two * EG
      Lambda = Ebulk - twoG * oneThrd
C
      IF (totalTime.EQ.0.0) THEN
C
         DO k = 1,nblock
            trace = strainInc(k,1) + strainInc(k,2) + strainInc(k,3)
            stressNew(k,1) = stressOld(k,1) + twoG*strainInc(k,1) +
     1                       Lambda*trace
            stressNew(k,2) = stressOld(k,2) + twoG*strainInc(k,2) +
     1                       Lambda*trace
            stressNew(k,3) = stressOld(k,3) + twoG*strainInc(k,3) +
     1                       Lambda*trace
            stressNew(k,4) = stressOld(k,4) + twoG*strainInc(k,4)
            stressNew(k,5) = stressOld(k,5) + twoG*strainInc(k,5)
            stressNew(k,6) = stressOld(k,6) + twoG*strainInc(k,6)
         END DO
C
         RETURN
C
      END IF
C

C
      DO k = 1,nblock
C
C        !CALCULATE LEFT CAUCHY-GREEN STRAIN TENSOR, B = U*U

         Bxx = stretchNew(k,1) * stretchNew(k,1) +
     1         stretchNew(k,4) * stretchNew(k,4) +
     2         stretchNew(k,6) * stretchNew(k,6)
         Byy = stretchNew(k,2) * stretchNew(k,2) +
     1         stretchNew(k,4) * stretchNew(k,4) +
     2         stretchNew(k,5) * stretchNew(k,5)
         Bzz = stretchNew(k,3) * stretchNew(k,3) +
     1         stretchNew(k,5) * stretchNew(k,5) +
     2         stretchNew(k,6) * stretchNew(k,6)
         Bxy = stretchNew(k,1) * stretchNew(k,4) +
     1         stretchNew(k,4) * stretchNew(k,2) +
     2         stretchNew(k,6) * stretchNew(k,5)
         Bxz = stretchNew(k,1) * stretchNew(k,6) +
     1         stretchNew(k,4) * stretchNew(k,5) +
     2         stretchNew(k,6) * stretchNew(k,3)
         Byz = stretchNew(k,4) * stretchNew(k,6) +
     1         stretchNew(k,2) * stretchNew(k,5) +
     2         stretchNew(k,5) * stretchNew(k,3)
C
C        !CALCULATE J = |F| = |U|

C
         J =    stretchNew(k,1) *
     1        ( stretchNew(k,2) * stretchNew(k,3)   -
     2          stretchNew(k,5) * stretchNew(k,5) ) +
     3          stretchNew(k,4) *
     4        ( stretchNew(k,5) * stretchNew(k,6)   -
     5          stretchNew(k,3) * stretchNew(k,4) ) +
     6          stretchNew(k,6) *
     7        ( stretchNew(k,4) * stretchNew(k,5)   -
     8          stretchNew(k,2) * stretchNew(k,6) )
C
C        !CALCULATE MODIFIED STRAIN TENSOR, Bbar = J^(-2/3)*B 

         Scale = J**(-twoThrd)
C
         BbarXX = Scale * Bxx
         BbarYY = Scale * Byy
         BbarZZ = Scale * Bzz
         BbarXY = Scale * Bxy
         BbarXZ = Scale * Bxz
         BbarYZ = Scale * Byz
C
C        !FIRST INVARIANT OF Bbar = tr(Bbar)

         I1 = BbarXX + BbarYY + BbarZZ
C
C        !DEVIATORIC PART OF Bbar

         p0 = oneThrd * I1
C
         dBbarXX = BbarXX - p0
         dBbarYY = BbarYY - p0
         dBbarZZ = BbarZZ - p0
         dBbarXY = BbarXY
         dBbarXZ = BbarXZ
         dBbarYZ = BbarYZ
C
C        !DERIVATIVES OF STRAIN-ENERGY FUNCTION

         IF ((I1 - three).LT.toler) THEN
            dWdI1 = zero
         ELSE
            dWdI1 = E  * A1 * (I1 - three)**(E-one) +
     1              EP * A2 * (I1 - three)**(EP-one)+
     2              QU * A3 * (I1 - three)**(QU-one)
         END IF
C
         dWdJ = two/D1 * (J - one)
C
C        !COROTATIONAL CAUCHY (TRUE) STRESSES

         g1 = two/J * dWdI1
C
         stressNew(k,1) = g1 * dBbarXX + dWdJ
         stressNew(k,2) = g1 * dBbarYY + dWdJ
         stressNew(k,3) = g1 * dBbarZZ + dWdJ
         stressNew(k,4) = g1 * dBbarXY
         stressNew(k,5) = g1 * dBbarYZ
         stressNew(k,6) = g1 * dBbarXZ

C
      END DO
C
      RETURN
	  end
