      SUBROUTINE PMABS
C///////////////////////////////////////////////////////////////////
C
C            ONE DIMENSIONAL LAMINAR PREMIXED FLAME CODE
C
C     WRITTEN BY:
C         ROBERT J. KEE
C         COMPUTATIONAL MECHANICS DIVISION
C         SANDIA NATIONAL LABORATORIES
C         LIVERMORE, CA  94550
C         (415) 294-3272
C
C/////////////////////////////////////////////////////////////////////
C
C     VERSION 2.5b
C
C   CHANGES FROM VERSION 1.0
C     1.  CHANGED REAL*8 TO DOUBLE PRECISION
C   CHANGES FROM VERSION 1.1
C     1.  READ LENICK, LENCK, LENCCK FROM LINKING FILE
C         (SUBROUTINE POINT) INSTEAD OF CALCULATING
C     2.  ALLOW UPPER/LOWER CASE INPUT
C   CHANGES FROM VERSION 1.2
C     1.  ALLOW MIXED CASE SPECIES NAMES TO REMAIN MIXED
C   CHANGES FROM VERSION 1.3
C     1.  LINKING FILE HAS ADDITIONAL VARIABLES KERR, MAXTB
C   CHANGES FROM VERSION 1.4
C     1.  CONVERSION TO MULTICOMPONENT GAS TRANSPORT
C     2.  INCLUDE SPECIES V CORRECTION AS OPTION
C
C   CHANGES TO VERSION 1.6
C     1. Modify POINT to use new subroutines CKLEN and TPLEN instead
C        of reading linking file.
C   CHANGES TO VERSION 1.7
C     1. Multicomponent version
C     2. Solution has additional work array records
C   CHANGES TO VERSION 1.8
C     1. Refine some DO loops
C   CHANGES TO VERSION 1.9
C     1. Restore file has additioal work array records, as in solution
C        file
C   CHANGES TO VERSION 2.0
C     1. Do not use linking file information stored in a restart
C        solution, as mechanism may have changed.
C     2. Call list for TWOPNT requires additional input; optional
C        use of new keywords reset default values:
C        'ISTP' n - sets NINIT initial time steps before newton
C                   (default is 0)
C        'IRET' n - set retirement age IRETIR of old time step
C                   (default 50)
C        'NJAC' n - set retirement age NJAC of Jacobian during
C                   steady state newton (default 20)
C        'TJAC' n - set retirement age ITJAC of Jacobian during
C                   time stepping (default 20)
C        'DTMX' x - set maximum time step DTMAX (default 1.0E-4)
C     3. Keyword SPOS sets a minimum value of species solutions.
C     4. Keyword NOFT allows skipping of the fixed temperature
C        problem.
C     5. Keywork GFAC allows scaling of production rates.
C     CHANGES FOR VERSION 2.1 (1/15/91 F. Rupley per Bob Kee)
C     1. Set default value of LVCOR to .TRUE. in Subroutine RDKEY
C     CHANGES FOR VERSION 2.3 (4/29/91 F. Rupley per Bob Kee)
C     1. Correction in MTNRPR for TDIF option.
C     CHANGES FOR VERSION 2.4 (5/9/91 F. Rupley per R. Kee)
C     1. Subroutine REHSEN and keyword "HSEN" for sensitivity to
C        heats of formation.
C     CHANGES FOR VERSION 2.5 (7/16/91 F. Rupley per R. Kee)
C     1. Correction for velocity calculation in Subroutine PRINT.
C     CHANGES FOR VERSION 2.5b (12/4/92 F. Rupley)
C     1. Upgrad to CKLIB.40
      END
C
C/////////////////////////////////////////////////////////////////////
C
      SUBROUTINE PREMIX (NMAX, LIN, LOUT, LINKCK, LINKMC, LREST, LSAVE,
     1                   LRCRVR, LENLWK, L, LENIWK, I, LENRWK, R,
     2                   LENCWK, C)
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION I(*), R(*)
      LOGICAL   L(*)
      CHARACTER C(*)*(*)
C
      COMMON /FLFLFL/ NCKW, NMCW, NEPS, NWT, NRE, NSCH, NX, NCON, NTGV,
     1                NXGV, ND, NDKJ, NTDR, NYV, NSCL, NABV, NBLW, NBUF,
     2                NTWP, NS, NSN, NF, NFN, NDS, NSSV, NKA6, NA, ICKW,
     3                IMCW, IKS, ICC, IKR, IKI, IKP, IIP, LAC, LMK
C
      COMMON /LOCS/ NT, NM, NYS, NY, NTR
C
C          WRITE VERSION NUMBER
C
      WRITE (LOUT,15)
   15 FORMAT (
     1/' PREMIX:  One-dimensional steady premixed laminar flame code',
     2/'          CHEMKIN-II Version 2.5b, December 1992',
C*****precision > double
     3/'          DOUBLE PRECISION')
C*****END precision > double
C*****precision > single
C     3/'          SINGLE PRECISION')
C*****END precision > single
C
C          SET UP INTERNAL WORK POINTERS
C
      CALL POINT (LINKCK, LINKMC, LENIWK, LENRWK, LENCWK, NMAX, LOUT,
     1            LSAVE, KK, II, NATJ, LENTWP, LTOT, ITOT, NTOT,
     2            ICTOT, I, R, C)
C
C           CHECK FOR ENOUGH SPACE
C
      WRITE (LOUT, 7000) LENLWK, LTOT, LENIWK, ITOT, LENRWK, NTOT,
     1                   LENCWK, ICTOT
7000  FORMAT (/,'                WORKING SPACE REQUIREMENTS',
     1        /,'                 PROVIDED        REQUIRED ',
     2        /,' LOGICAL  ' , 2I15,
     3        /,' INTEGER  ' , 2I15,
     4        /,' REAL     ' , 2I15,
     5        /,' CHARACTER' , 2I15,/)
C
      IF (LTOT.GT.LENLWK .OR. ITOT.GT.LENIWK .OR. NTOT.GT.LENRWK
     1                   .OR. ICTOT.GT.LENCWK) THEN
         WRITE (LOUT, *) '  FATAL ERROR, NOT ENOUGH WORK SPACE PROVIDED'
         STOP
      ENDIF
C
      CALL FLDRIV (LENTWP, LIN, LOUT, LREST, LSAVE, LRCRVR, 
     1             KK, II,
     1             NATJ, NMAX, R(NCKW), R(NMCW), R(NEPS), R(NWT),
     2             R(NRE), R(NSCH), R(NX), R(NCON), R(NTGV), R(NXGV),
     3             R(ND), R(NDKJ), R(NTDR), R(NYV), R(NSCL), R(NABV),
     4             R(NBLW), R(NBUF), R(NTWP), R(NS), R(NSN), R(NF),
     5             R(NFN), R(NDS), R(NSSV), R(NKA6), R(NA), I(ICKW),
     6             I(IMCW), C(IKS), C(ICC), I(IKR), I(IKI), I(IKP),
     7             I(IIP), L(LAC),  L(LMK))
C
      RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE POINT (LINKCK, LINKMC, LENIWK, LENRWK, LENCWK, NMAX,
     1                  LOUT, LSAVE, KK, II, NATJ, LENTWP, LTOT, ITOT,
     2                  NTOT, ICTOT, I, R, C)
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION I(*), R(*)
      CHARACTER C(*)*(*)
C
      COMMON /FLFLFL/ NCKW, NMCW, NEPS, NWT, NRE, NSCH, NX, NCON, NTGV,
     1                NXGV, ND, NDKJ, NTDR, NYV, NSCL, NABV, NBLW, NBUF,
     2                NTWP, NS, NSN, NF, NFN, NDS, NSSV, NKA6, NA, ICKW,
     3                IMCW, IKS, ICC, IKR, IKI, IKP, IIP, LAC, LMK
C
      COMMON /LOCS/ NT, NM, NYS, NY, NTR
C
C          SET THE POINTERS INTO THE SOLUTION VECTOR
C
      NT  = 1
      NYS = 1
      NY  = 2
      NM  = KK+2
C
      CALL CKLEN (LINKCK, LOUT, LENICK, LENRCK, LENCCK)
      CALL MCLEN (LINKMC, LOUT, LENIMC, LENRMC)
C
C     real chemkin work space
      NCKW = 1
C     real transport work space
      NMCW = NCKW + LENRCK
      NTOT = NMCW + LENRMC
C
C     integer chemkin work space
      ICKW = 1
C     integer transport work space
      IMCW = ICKW + LENICK
      ITOT = IMCW + LENIMC
C
C     character chemkin work space
      ICC = 1
      ICTOT = ICC + LENCCK
C
C
      IF (ITOT.LT.LENIWK .AND. NTOT.LT.LENRWK .AND. ICTOT.LT.LENCWK)
     1   THEN
         CALL CKINIT (LENICK, LENRCK, LENCCK, LINKCK, LOUT, I, R, C)
         CALL CKINDX (I, R, MM, KK, II, NFIT)
C
         CALL CKMXTP (I, MAXTP)
         NTR = MAXTP - 1
C
         CALL MCINIT (LINKMC, LOUT, LENIMC, LENRMC, I(IMCW), R(NMCW))
         REWIND LSAVE
         CALL PRSAVE (I, R, C, I(IMCW), R(NMCW),
     1                LOUT, LSAVE)
      ENDIF
C
      NATJ = KK+2
      LENTWP = (7*NATJ+2)*NMAX
C
C          APPORTION THE BALANCE OF THE FLOATING POINT SPACE
C
      NEPS = NTOT
      NWT  = NEPS + KK
      NRE  = NWT  + KK
      NSCH = NRE  + KK
      NX   = NSCH + KK*6
      NCON = NX   + NMAX
      NTGV = NCON + NMAX
      NXGV = NTGV + NMAX
      ND   = NXGV + NMAX
      NDKJ = ND   + KK*NMAX
      NTDR = NDKJ + KK*KK*NMAX
      NYV  = NTDR + KK*NMAX
      NSCL = NYV  + KK*NMAX
      NABV = NSCL + NATJ
      NBLW = NABV + NATJ*NMAX
      NBUF = NBLW + NATJ*NMAX
      NTWP = NBUF + NATJ*NMAX
      NS   = NTWP + LENTWP
      NSN  = NS   + NATJ*NMAX
      NF   = NSN  + NATJ*NMAX
      NFN  = NF   + NATJ*NMAX
      NDS  = NFN  + NATJ*NMAX
      NSSV = NDS  + NMAX
C     thermodynamic coefficients a6
      NKA6 = NSSV + NMAX
      NA   = NKA6 + NTR
      NTOT = NA   + (6*NATJ-2) * (NATJ*NMAX) - 1
C
C           APPORTION THE BALANCE OF THE INTEGER SPACE
C
      IKR  = ITOT
      IKI  = IKR  + KK
      IKP  = IKI  + KK
      IIP  = IKP  + KK
      ITOT  = IIP  + NATJ*NMAX - 1
C
C           APPORTION THE LOGICAL SPACE
C
      LAC  = 1
      LMK  = LAC  + NATJ
      LTOT = LMK  + NMAX - 1
C
C
C           APPORTION THE BALANCE OF THE CHARACTER SPACE
C
      IKS = ICTOT
      ICTOT = IKS + KK - 1
      RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE FLDRIV (LENTWP, LIN, LOUT, LREST, LSAVE, LRCRVR, 
     1                   KK, II, NATJ, NMAX, RCKWRK, 
     2                   RMCWRK, EPS, WT, REAC, SCRCHK, X, COND, 
     3                   TGIVEN, XGIVEN, D, DKJ, TDR, YV, SCAL, ABOVE, 
     4                   BELOW, BUFFER, TWPWK, S, SN, F, FN, DS, SSAVE, 
     5                   A6, A, ICKWRK, IMCWRK, KSYM, CCKWRK, KR, KI, 
     6                   KP, IP, ACTIVE, MARK)
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
C
      COMMON /LOCS/ NT, NM, NYS, NY, NTR
C
      DIMENSION ICKWRK(*), RCKWRK(*), IMCWRK(*), RMCWRK(*)
C
      DIMENSION EPS(*), WT(*), REAC(*), X(*), COND(*), D(KK,*),
     1          TDR(KK,*), YV(KK,*), SCAL(*), S(NATJ,*), SN(NATJ,*),
     2          F(NATJ,*), XGIVEN(*), TGIVEN(*), SCRCHK(KK,*), KR(*),
     3          KI(*), KP(*), DKJ(KK,KK,*), A6(*)
C
      CHARACTER*(*) CCKWRK(*), KSYM(*)
      CHARACTER*16 ICHR, ISOLUT, ISENSI, ICKLNK, IMCLNK, IHSENS
C
C         DIMENSION NEWTON SPACE
C
      DIMENSION ABOVE(NATJ,*), BELOW(NATJ,*), BUFFER(NATJ,*), TWPWK(*),
     1          A(6 * NATJ - 2, *), IP(NATJ,*), FN(NATJ,*), DS(*),
     2          SSAVE(*), LEVEL(2)
C
      LOGICAL LBURNR, LTIME, LTIME2, LMOLE, LENRGY, LTDIF, LUMESH,
     1        LRSTRT, LCNTUE, LASEN, LRSEN, LESEN, LMSEN, LPSEN, LMULTI,
     2        LVARMC, RSTCNT, ERROR, FUNCTN, JACOBN, REENTR, SOLVE,
     3        STORE, SUCCES, ADAPT, SHOW, SAVE, UPDATE, ACTIVE(*),
     4        MARK(*), ENERGY, BURNER, LUSTGV, IERR, LVCOR, LHSEN
C
      INTEGER CALL, CALLS
C
      DATA ADAFLR /1.0E-08/, LCNTUE/.FALSE./
      DATA ISOLUT/'SOLUTION        '/, ISENSI/'SENSITIVITY     '/,
     1     ICKLNK/'CKLINK          '/, IMCLNK/'MCLINK          '/,
     2     IHSENS/'HSENSITIVITY    '/

C
C          COMPUTE THE UNIT ROUNDOFF, AND THE RELATIVE AND ABSOLUTE
C            PERTURBATIONS FOR THE JACOBIAN EVALUATION.
C
      U = 1.0
10    CONTINUE
      U = U*0.5
      COMP = 1.0 + U
      IF (COMP .NE. 1.0) GO TO 10
      ABSOL = SQRT(2.0*U)
      RELAT = SQRT(2.0*U)
C
      RSTCNT = .FALSE.
C
      NT = 1
      NYS = 1
      NY = 2
      NM = KK+2
C
      CALL CKSYMS (CCKWRK, LOUT, KSYM, IERR)
      CALL CKWT   (ICKWRK, RCKWRK, WT)
      CALL CKRP   (ICKWRK, RCKWRK, RU, RUC, PATM)
C
200   CONTINUE
C
C         CALL THE KEYWORD INPUT
C
      CALL RDKEY (KK, NMAX, NATJ, LIN, LOUT, KSYM, PATM, LBURNR, LTIME,
     1            LTIME2, LMOLE, LUSTGV, LENRGY, LMULTI, LVCOR, LTDIF,
     2            LUMESH, LRSTRT, LCNTUE, MFILE, IPRNT, UFAC, DFAC,
     3            DTMAX, DTMIN, LASEN, LHSEN, LRSEN, LESEN, LMSEN,
     4            LPSEN, LCSEN, LDSEN, P, JJ, NTOT, NADP, X, SCAL,
     5            NREAC, NINTM, NPROD, REAC, SCRCHK(1,2), SCRCHK(1,3),
     6            KR, KI, KP, XSTR, XCEN, XEND, WMIX, FLRT, GRAD,
     7            CURV, SFLR, NTEMP, XGIVEN, TGIVEN, TFIXT, ATOL,
     8            RTOL, ATIM, RTIM, NJAC, ITJAC, NINIT, NUMDT, IRETIR,
     9            DT1, NUMDT2, DT2, WNDFAC, GFAC, SPOS, N1CALL)
C
C         SET THE SOLUTION BOUNDS
C
      DO 100 J = 1, NMAX
        BELOW(NT,J) = 50.
        ABOVE(NT,J) = 6000.0
        BELOW(NM,J) = 0.0
        ABOVE(NM,J) = 10000.0
  100 CONTINUE
      DO 175 J = 1, NMAX
         DO 150 K = 1, KK
            BELOW (NYS+K, J) = SFLR
            ABOVE (NYS+K, J) = 10.
150      CONTINUE
175   CONTINUE
C
C         SET THE VARIABLES ON WHICH TO ADAPT
C
      ACTIVE(NM) = .FALSE.
      ACTIVE(NT) = LENRGY
      DO 300 K = 1, KK
         ACTIVE(NYS+K) = .TRUE.
300   CONTINUE
C
      IF (IPRNT .LT. 10) THEN
         LEVEL(1) = MIN (2, IPRNT)
         LEVEL(2) = LEVEL(1)
      ELSE
         LEVEL(1) = IPRNT/10
         LEVEL(2) = IPRNT - 10*LEVEL(1)
         LEVEL(1) = MAX( LEVEL(1), LEVEL(2) )
      ENDIF
C
      IF (LRSTRT) THEN
         IF (.NOT. RSTCNT) THEN
C
C          THIS IS A RESTART
C
            NREST = 0
  320       CONTINUE
            READ (LREST, END=350, ERR=350) ICHR
            IF (ICHR .EQ. ICKLNK) THEN
               DO 321 L = 1, 4
                  READ (LREST, END=350, ERR=350)
  321          CONTINUE
            ELSEIF (ICHR .EQ. IMCLNK) THEN
               DO 322 L = 1, 3
                  READ (LREST, END=350, ERR=350)
  322          CONTINUE
            ELSEIF (ICHR .EQ. ISOLUT) THEN
               READ (LREST, END=350, ERR=350) NNNN, JJ, PDUM, DMFLRT
               READ (LREST, END=350, ERR=350) (X(J), J=1,JJ)
               READ (LREST, END=350, ERR=350) ((S(N,J), N=1,NNNN),
     1               J=1,JJ)
               IF (NNNN .NE. NATJ) THEN
                  WRITE (LOUT, *)
     1                   ' FATAL ERROR, INCOMPATIBLE RESTART FILE'
                  STOP
               ENDIF
               NREST = NREST + 1
               IF (NREST .EQ. MFILE) GO TO 350
            ELSEIF (ICHR .EQ. ISENSI) THEN
               DO 310 I = 1, II
                  READ (LREST, END=350, ERR=350)
     1            IS, ((SN(N,J),N=1,NATJ), J=1,JJ)
  310          CONTINUE
            ELSEIF (ICHR .EQ. IHSENS) THEN
               DO 315 K = 1, II
                  READ (LREST, END=350, ERR=350)
     1            KS, ((SN(N,J),N=1,NATJ), J=1,JJ)
  315          CONTINUE
            ELSE
               WRITE (LOUT, *)
     1         'FATAL ERROR, NOT A SOLUTION ON RESTART FILE'
               STOP
            ENDIF
            GO TO 320
  350       CONTINUE
            IF (NREST .NE. MFILE) THEN
               WRITE (LOUT, *) ' Error reading solution file...'
               STOP
            ENDIF
         ENDIF
C
         CALL RESTRT (KK, NMAX, NATJ, JJ, LOUT, LMOLE, LUMESH, LBURNR,
     1                LUSTGV, NREAC, NINTM, NPROD, REAC, SCRCHK(1,2),
     2                SCRCHK(1,3), KR, KI, KP, FLRT, NTEMP, XGIVEN,
     3                TGIVEN, XSTR, XCEN, XEND, WMIX, ICKWRK, RCKWRK,
     4                SCRCHK(1,4), COND, EPS, JFIXT, TFIXT, X, S)
C
      ELSE
C
C         SET THE STARTING PROFILES
C          (in the following, COND is used as length njj scratch space)
C
         CALL START (KK, NMAX, NATJ, JJ, LOUT, LMOLE, LUMESH, LBURNR,
     1               NREAC, NINTM, NPROD, REAC, SCRCHK(1,2),
     2               SCRCHK(1,3), KR, KI, KP, FLRT, NTEMP, XGIVEN,
     3               TGIVEN, XSTR, XCEN, XEND, WMIX, ICKWRK, RCKWRK,
     4               SCRCHK(1,4), COND, EPS, JFIXT, TFIXT, X, S)
C
      ENDIF
C
C///  DECIDE HOW MANY TIMES TO CALL TWOPNT.
C
      ENERGY = LENRGY
      BURNER = LBURNR
      IF (ENERGY) THEN
         CALLS = 2
      ELSE
         CALLS = 1
      ENDIF
C
C///  TOP OF THE LOOP OVER CALLS TO TWOPNT.
C
C          SET XFIXT AT TFIXT
C
      IF (.NOT. LBURNR) XFIXT = X(JFIXT)
C
      DO 1100 CALL = N1CALL, CALLS
C
         IF (CALL .EQ. 2) WRITE (LOUT,'(/A/)')
     1   ' FLDRIV: FINISHED FIXED TEMPERATURE, ADDING ENERGY EQUATION'
C
         IF (ENERGY) THEN
            IF (CALL .EQ. 1) THEN
               LENRGY = .FALSE.
               LBURNR = .TRUE.
               ADAPT = .FALSE.
            ELSE
               LENRGY = .TRUE.
               LBURNR = BURNER
               ADAPT = .TRUE.
            ENDIF
         ELSE
            ADAPT = .TRUE.
         ENDIF
C
         IF (LTIME2 .AND. CALL.GT.1) THEN
            NUMDT = NUMDT2
            DT1 = DT2
         ENDIF
C
         ACTIVE(NT) = LENRGY
         REENTR = .FALSE.
C
         IPASSS = 1
500      CONTINUE
C
         CALL TWOPNT (ERROR, LOUT, LEVEL, NATJ*NMAX, IP, LENTWP, TWPWK,
     1                ABOVE, ACTIVE, ADAPT, BELOW, 0, BUFFER, NATJ,
     2                CONDIT, FUNCTN, JACOBN, MARK, X, IPASS, IPASSS,
     3                NADP, NMAX, JJ, REENTR, IREPRT, SAVE, 0,
     4                SHOW, SOLVE, ATOL, NJAC, RTOL, NINIT, NUMDT,
     5                IRETIR, STORE, SUCCES, ATIM, ITJAC, DFAC, RTIM,
     6                LTIME, UFAC, DTMAX, DTMIN, ADAFLR, GRAD, CURV,
     7                DT1, DT, UPDATE, S)
C
         IF (ERROR) THEN
            STOP
         ELSEIF (.NOT.SUCCES .AND. .NOT.REENTR) THEN
            IF (IREPRT .EQ. 2) WRITE (LOUT,*)
     1   '     TWOPNT requires more mesh points, but NMAX too small'
C
         ELSEIF (REENTR) THEN
C
            IF (FUNCTN) THEN
C
               LVARMC = .TRUE.
               CALL FUN (KK, JJ, NATJ, LBURNR, LENRGY, LMULTI, LVCOR,
     1                   LTDIF, LVARMC, LTIME, JFIXT, TFIXT, FLRT, P,
     2                   WT, EPS, DT, NTEMP, XGIVEN, TGIVEN, X, SN,
     3                   BUFFER, WNDFAC, SCRCHK(1,1), YV, SCRCHK(1,2),
     4                   SCRCHK(1,3), SCRCHK(1,4), COND, D, DKJ, TDR,
     5                   ICKWRK, RCKWRK, IMCWRK, RMCWRK, F, SCRCHK(1,5),
     6                   GFAC)
C
C*****precision > double
              CALL DCOPY (NATJ*JJ, F, 1, BUFFER, 1)
C*****END precision > double
C
C*****precision > single
C              CALL SCOPY (NATJ*JJ, F, 1, BUFFER, 1)
C*****END precision > single
C
           ELSEIF (JACOBN) THEN
C
              CALL JACOB (KK, JJ, NATJ, LBURNR, LENRGY, LMULTI, LVCOR,
     1                    LTDIF, LVARMC, LTIME, JFIXT, TFIXT, FLRT, P,
     2                    WT, EPS, DT, NTEMP, XGIVEN, TGIVEN, X, SN,
     3                    BUFFER, WNDFAC, ABSOL, RELAT, SCRCHK, YV,
     4                    COND, D, DKJ, TDR, ICKWRK, RCKWRK, IMCWRK,
     5                    RMCWRK, F, FN, A, DS, SSAVE, GFAC)
C
C*****precision > double
             CALL DGBCO
C*****END precision > double
C*****precision > single
C             CALL SGBCO
C*****END precision > single
     +          (A, 6 * NATJ - 2, NATJ * JJ, 2 * NATJ - 1, 2 * NATJ - 1,
     +          IP, RCOND, FN)
C
             IF (RCOND .LE. 0.0) THEN
                WRITE (LOUT,*) ' FATAL ERROR, SINGULAR JACOBIAN '
                STOP
             ENDIF
             CONDIT = 1.0 / RCOND
C
           ELSEIF (SOLVE) THEN
C
C*****precision > double
              CALL DGBSL
C*****END precision > double
C*****precision > single
C             CALL SGBSL
C*****END precision > single
     +        (A, 6 * NATJ - 2, NATJ * JJ, 2 * NATJ - 1, 2 * NATJ - 1,
     +        IP, BUFFER, 0)
C
           ELSEIF (STORE) THEN
C
C*****precision > double
              CALL DCOPY (NATJ*JJ, BUFFER, 1, SN, 1)
C*****END precision > double
C
C*****precision > single
C              CALL SCOPY (NATJ*JJ, BUFFER, 1, SN, 1)
C*****END precision > single
C
              IF (SPOS .GE. 0.0) THEN
                 DO 1040 J = 1, JJ
                    DO 1035 K = 1, KK
                       SN(NYS+K,J) = MAX (SN(NYS+K,J), SPOS)
 1035               CONTINUE
 1040            CONTINUE
              ENDIF
C
           ELSEIF (SHOW) THEN
              CALL PRINT (LOUT, KK, JJ, NATJ, LMOLE, P, X, BUFFER,
     1                    SCRCHK(1,1), SCRCHK(1,2), KSYM, ICKWRK,
     2                    RCKWRK)
C
           ELSEIF (SAVE) THEN
              REWIND LRCRVR
              CALL PRSAVE (ICKWRK, RCKWRK, CCKWRK, 
     1                     IMCWRK, RMCWRK, LOUT, LRCRVR)
              WRITE (LRCRVR) ISOLUT
              WRITE (LRCRVR) NATJ, JJ, P, BUFFER(NM,1)
              WRITE (LRCRVR) (X(J), J=1,JJ)
              WRITE (LRCRVR) ((BUFFER(N,J), N=1,NATJ), J=1,JJ)
C
           ELSEIF (UPDATE) THEN
              IF (.NOT. LENRGY) THEN
                 DO 900 J = 1, JJ
                    CALL TEMP (NTEMP, X(J), XGIVEN, TGIVEN, TI)
                    BUFFER(NT,J) = TI
900              CONTINUE
              ENDIF
              IF (.NOT. LBURNR) THEN
                 DO 950 J = 1, JJ
                    IF (X(J) .EQ. XFIXT) JFIXT=J
950              CONTINUE
              ENDIF
              DO 1000 J = 2, JJ
                 BUFFER(NM,J) = BUFFER(NM,1)
1000          CONTINUE
           ENDIF
C
           GO TO 500
C
         ENDIF
C
1100  CONTINUE
C
C              WRITE TO LSAVE WHEN SOLUTION IS COMPLETE
C
      WRITE (LSAVE) ISOLUT
      WRITE (LSAVE) NATJ, JJ, P, S(NM,1)
      WRITE (LSAVE) (X(J), J=1,JJ)
      WRITE (LSAVE) ((S(N,J), N=1,NATJ), J=1,JJ)
C
C
      IF (LASEN) THEN
C
         LDA = 6*NATJ - 2
         CALL   REASEN (II, KK, JJ, NATJ, LDA, LBURNR, LENRGY, LMULTI,
     1                  LTDIF, LMOLE, LSAVE, LRCRVR, LOUT, LVARMC,
     2                  LTIME, JFIXT, TFIXT, FLRT, P, WT, EPS, DT,
     3                  NTEMP, XGIVEN, TGIVEN, X, SN, S, WNDFAC, ABSOL,
     4                  RELAT, SCRCHK, YV, COND, D, DKJ, TDR, ICKWRK,
     5                  RCKWRK, IMCWRK, RMCWRK, F, FN, A, DS, SSAVE,
     6                  IP, GFAC)
C
         WRITE (LOUT,'(/A/)')' SENSITIVITY CALCULATION COMPLETE'
C
      ENDIF
C
      IF (LHSEN) THEN
        LDA = 6*NATJ - 2
        CALL REHSEN (II, KK, JJ, NATJ, LDA, LBURNR, LENRGY, LMULTI,
     1               LTDIF, LMOLE, LSAVE, LRCRVR, LOUT, LVARMC,
     2               LTIME, JFIXT, TFIXT, FLRT, P, WT, EPS, DT,
     3               NTEMP, XGIVEN, TGIVEN, X, SN, S, WNDFAC, ABSOL,
     4               RELAT, SCRCHK, YV, COND, D, DKJ, TDR, ICKWRK,
     5               RCKWRK, IMCWRK, RMCWRK, F, FN, A, DS, SSAVE,
     6               IP, GFAC, A6)
C
         WRITE (LOUT,'(/A/)')' H SENSITIVITY CALCULATION COMPLETE'
C
      ENDIF
C
C            CHECK FOR CONTINUATION
C
      IF (LCNTUE) THEN
C
         WRITE (LOUT,'(/////)')
         DO 1210 L = 1, 5
            WRITE (LOUT,*)
     1    ' ////////////////// CONTINUING TO NEW PROBLEM /////////////'
1210     CONTINUE
         WRITE (LOUT,'(/////)')
C
         RSTCNT = .TRUE.
         LRSTRT = .TRUE.
         GO TO 200
      ENDIF
C
      STOP
      END
C
C---------------------------------------------------------------------
C
      SUBROUTINE JACOB (KK, POINTS, COMPS, LBURNR, LENRGY, LMULTI,
     1                  LVCOR, LTDIF, LVARMC, LTIME, JFIXT, TFIXT, FLRT,
     2                  P, WT, EPS, DT, NTEMP, XGIVEN, TGIVEN, MESH,
     3                  SN, X0, WNDFAC, ABSOL, RELAT, SCRTCH, YV, COND,
     4                  D, DKJ, TDR, ICKWRK, RCKWRK, IMCWRK, RMCWRK, Y1,
     5                  Y0, A, PERTRB, SAVE, GFAC)
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INTEGER CLASS, COMP, COMPS, INDEX, LOWER, POINT, POINTS, ROWS,
     1        UPPER
      LOGICAL LBURNR, LENRGY, LTDIF, LVARMC, LTIME, LMULTI, LVCOR
      REAL MESH
C
      PARAMETER (LOWER = 1, UPPER = 1, MODULS = LOWER + 1 + UPPER)
C
      DIMENSION A((2 * LOWER + UPPER + 3) * COMPS - 2, COMPS * POINTS),
     1   MESH(POINTS), PERTRB(POINTS), SAVE(POINTS), SCRTCH(KK, 5),
     2   X0(COMPS, POINTS), Y0(COMPS, POINTS), Y1(COMPS, POINTS)
C
      KOFSET = (LOWER + UPPER + 2) * COMPS - 1
      ROWS = (2 * LOWER + UPPER + 3) * COMPS - 2
C
C///  ZERO THE MATRIX STORAGE SPACE.
C
C
      ZERO = 0.0
C*****precision > double
      CALL DCOPY (ROWS * COMPS * POINTS, ZERO, 0, A, 1)
C*****END precision > double
C
C*****precision > single
C      CALL SCOPY (ROWS * COMPS * POINTS, ZERO, 0, A, 1)
C*****END precision > single
C
C
C///  CALL THE FUNCTION AT X0 AND STORE IN Y0.
C
      CALL FUN (KK, POINTS, COMPS, LBURNR, LENRGY, LMULTI, LVCOR,
     1          LTDIF, .TRUE., LTIME, JFIXT, TFIXT, FLRT, P, WT, EPS,
     2          DT, NTEMP, XGIVEN, TGIVEN, MESH, SN, X0, WNDFAC,
     3          SCRTCH(1,1), YV, SCRTCH(1,2), SCRTCH(1,3), SCRTCH(1,4),
     4          COND, D, DKJ, TDR, ICKWRK, RCKWRK, IMCWRK, RMCWRK, Y0,
     5          SCRTCH(1,5), GFAC)
C
C///  TOP OF THE LOOPS OVER THE RESIDUE CLASSES AND SOLUTION COMPONENTS.
C
      DO 0400 CLASS = 1, MODULS
         DO 0400 COMP = 1, COMPS
C
C///  FOR A GIVEN RESIDUE CLASS AND A GIVEN SOLUTION COMPONENT,
C///  PERTRB THE X0 VECTOR AT POINTS IN THE SAME RESIDUE CLASS.
C
            INDEX = 0
            DO 0100 POINT = CLASS, POINTS, MODULS
               INDEX = INDEX + 1
               SAVE(INDEX) = X0(COMP, POINT)
               PERTRB(INDEX) = ABS(X0(COMP, POINT)) * RELAT + ABSOL
               X0(COMP, POINT) = X0(COMP, POINT) + PERTRB(INDEX)
0100        CONTINUE
C
C///  CALL THE FUNCTION AT THE PERTRBED X0 AND STORE THE RESULT IN Y1.
C
         CALL FUN (KK, POINTS, COMPS, LBURNR, LENRGY, LMULTI, LVCOR,
     1             LTDIF, .FALSE., LTIME, JFIXT, TFIXT, FLRT, P, WT,
     2             EPS, DT, NTEMP, XGIVEN, TGIVEN, MESH, SN, X0,
     3             WNDFAC, SCRTCH(1,1), YV, SCRTCH(1,2), SCRTCH(1,3),
     4             SCRTCH(1,4), COND, D, DKJ, TDR, ICKWRK, RCKWRK,
     5             IMCWRK, RMCWRK, Y1, SCRTCH(1,5), GFAC)
C
C///  RESTORE X0 TO ITS ORIGINAL VALUE.
C
            INDEX = 0
            DO 0200 POINT = CLASS, POINTS, MODULS
               INDEX = INDEX + 1
               X0(COMP, POINT) = SAVE(INDEX)
0200        CONTINUE
C
C///  DIFFERENCE TO GET THE COLUMNS OF THE JACOBIAN.
C
            INDEX = 0
            DO 0300 POINT = CLASS, POINTS, MODULS
               INDEX = INDEX + 1
               TEMP = 1.0 / PERTRB(INDEX)
               K = COMP + (POINT - 1) * COMPS
               DO 0300 J2 = MAX (POINT - LOWER, 1),
     +                      MIN (POINT + UPPER, POINTS)
                  JOFSET = (J2 - 1) * COMPS - K + KOFSET
                  DO 0250 J1 = 1, COMPS
                     A(J1 + JOFSET, K) =
     1                  (Y1(J1, J2) - Y0(J1, J2)) * TEMP
0250              CONTINUE
0300        CONTINUE
C
C///  BOTTOM OF THE LOOPS OVER THE RESIDUE CLASSES AND SOLUTION
C///  COMPONENTS.
C
0400  CONTINUE
C
      RETURN
      END
C
C---------------------------------------------------------------------
C
      SUBROUTINE FUN (KK, JJ, NATJ, LBURNR, LENRGY, LMULTI, LVCOR,
     1                LTDIF, LVARMC, LTIME, JFIXT, TFIXT, FLRT, P, WT,
     2                EPS, DT, NTEMP, XGIVEN, TGIVEN, X, SN, S, WNDFAC,
     3                YAV, YV, WDOT, CP, H, COND, D, DKJ, TDR, ICKWRK,
     4                RCKWRK, IMCWRK, RMCWRK, F, XAV, GFAC)
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /LOCS/ NT, NM, NYS, NY, NTR
C
      DIMENSION WT(*), EPS(*), XGIVEN(NTEMP), TGIVEN(NTEMP), X(*),
     1          S(NATJ,*), SN(NATJ,*), YAV(*), YV(KK,*), WDOT(*),
     2          CP(*), H(*), COND(*), D(KK,*), TDR(KK,*), ICKWRK(*),
     3          RCKWRK(*), IMCWRK(*), RMCWRK(*), F(NATJ,*),
     4          DKJ(KK,KK,*), XAV(*)
C
      LOGICAL LBURNR, LENRGY, LTDIF, LVARMC, LTIME, LMULTI, LVCOR
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   KK     - NUMBER OF CHEMICAL SPECIES.
C   JJ     - NUMBER OF MESH POINTS.
C   NATJ   - NUMBER OF DEPENDENT VARIABLES AT EACH MESH POINT, AND
C            THE EXACT FIRST DIMENSION OF S(NATJ,*), SN(NATJ,*), AND
C            F(NATJ,*).  NATJ=KK+2.
C   LBURNR - IF LBURNR=.TRUE. THEN THIS IS A BURNER STABILIZED FLAME
C            IF LBURNR=.FALSE. THIS IS A FREELY PROPAGATING ADIABATIC.
C   LENRGY - IF LENRGY=.TRUE. THEN THE ENERGY EQUATION IS TO BE SOLVED,
C            OTHERWISE A FIXED TEMPERATURE PROFILE IS USED.
C   LTDIF  - IF LTDIF=.TRUE. THEN EVALUATE THERMAL DIFFUSION RATIOS AS
C            WELL AS DIFFUSION COEFFICIENTS.
C   LMULTI - IF LMULTI=.TRUE. THEN MULTICOMPONENT FORMULAS USED
C            IF LMULTI=.FALSE. THEN MIXTURE-AVERAGED FORMULAS USED
C   LVCOR  - IF LVCOR=.TRUE. THEN USE CORRECTION VELOCITY FORMULISM
C            IF LVCOR=.FALSE. THEN USE 'TRACE' APPROXIMATION, LUMPING
C            ALL TRANSPORT ERRORS INTO THE 'LAST' SPECIES.
C   LVARMC - IF LVARMC=.TRUE. THEN NEW TRANSPORT PROPERTIES WILL BE
C            COMPUTED, OTHERWISE THE PREVIOUSLY STORED VALUES WILL
C            BE USED.
C   LTIME  - IF LTIME=.TRUE. THEN A TIME STEP OF DT WILL BE ADDED INTO
C            THE RESIDUAL.
C   JFIXT  - THE MESH POINT AT WHICH THE TEMPERATURE IS FIXED IN AN
C            ADIABATIC FLAME COMPUTATION.
C   TFIXT  - THE TEMPERATURE THAT IS FIXED AT THE JFIXT MESH POINT IN
C            AN ADIABATIC FLAME COMPUTATION.
C   FLRT   - THE MASS FLOW RATE IN A BURNER STABILIZED PROBLEM.
C              CGS UNITS - GM/(CM**2-SEC)
C   P      - PRESSURE.
C              CGS UNITS - DYNES/CM**2
C   WT     - THE ARRAY OF SPECIES MOLECULAR WEIGHTS.
C              CGS UNITS - GM/MOLE
C              DIMENSION WT(*) AT LEAST KK.
C   EPS    - THE MASS FLUX FRACTIONS AT THE BURNER.
C              DIMENSION EPS(*) AT LEAST KK.
C   DT     - THE TIME STEP THAT IS USED IF LTIME=.TRUE.
C              CGS UNITS - SEC
C   NTEMP  - THE NUMBER OF TEMPERATURE-POSITION PAIRS GIVEN FOR BURNER
C            STABILIZED PROBLEMS WHERE THE ENERGY EQUATION IS NOT
C            SOLVED, AND THE TEMPERATURES ARE INTERPOLATED FROM FIXED
C            PROFILES.
C   XGIVEN - THE DISTANCES FROM THE BURNER AT WHICH TEMPERATURES ARE
C            SPECIFIED.
C              DIMENSION XGIVEN(*) AT LEAST NTEMP.
C   TGIVEN - THE TEMPERATURES AT THE XGIVEN LOCATIONS.
C              DIMENSION TGIVEN(*) AT LEAST NTEMP.
C   X      - THE ARRAY OF MESH POINT LOCATIONS.
C              DIMENSION X(*) AT LEAST JJ
C              CGS UNITS - CM
C   SN     - DEPENDENT VARIABLE MATRIX AT PREVIOUS TIMESTEP.  STORED
C            IN THE SAME FASHION AS S (BELOW)
C   S      - DEPENDENT VARIABLE MATRIX.  THE TEMPERATURES ARE STORED IN
C            T(J)=S(NT,J), THE MASS FRACTIONS ARE IN Y(K,J)=S(NYS+K,J),
C            AND THE FLOW RATES ARE IN FLRT(J)=S(NM,J).
C              DIMENSION S(NATJ,*) EXACTLY NATJ FOR THE FIRST DIMENSION,
C              AND AT LEAST JJ FOR THE SECOND.
C   WNDFAC - IF WINDFAC=1, THEN WINDWARD DIFFERENCING ON CONVECTION.
C            IF WINDFAC=0, THEN CENTRAL DIFFERENCING ON CONVECTION.
C
C WORK AND SCRATCH SPACE
C   YAV    - ARRAY OF MASS FRACTIONS AT MESH MIDPOINTS.  YAV(J) IS THE
C            MASS FRACTION BETWEEN J AND J+1.
C              DIMENSION Y(*) AT LEAST KK.
C   XAV    - ARRAY OF MOLE FRACTIONS AT MESH MIDPOINTS.  XAV(J) IS THE
C            MOLE FRACTION BETWEEN J AND J+1.
C              DIMENSION X(*) AT LEAST KK.
C   YV     - MATRIX OF MASS FRACTIONS TIMES DIFFUSION VELOCITIES AT THE
C            MESH MIDPOINTS.  YV(K,J) IS THE FLUX OF KTH SPECIES BETWEEN
C            J AND J+1.
C   WDOT   - ARRAY OF CHEMICAL PRODUCTION RATES.
C              CGS UNITS - MOLES/(CM**3-SEC)
C              DIMENSION WDOT(*) AT LEAST KK.
C   CP     - ARRAY OF SPECIES SPECIFIC HEATS.
C              CGS UNITS - ERGS/GM-K
C              DIMENSION CP(*) AT LEAST KK
C   H      - ARRAY OF SPECIES ENTHALPIES.
C              CGS UNITS - ERGS/MOLE
C              DIMENSION H(*) AT LEAST KK.
C   COND   - ARRAY OF THERMAL CONDUCTIVITIES AT THE MESH MIDPOINTS.
C            IF LVARMC=.TRUE. THESE ARE COMPUTED EACH TIME THE
C            FUNCTION IS CALLED, OTHERWISE THE STORED VALUES ARE USED.
C              CGS UNITS - ERG/CM-K-SEC
C              DIMENSION COND(*) AT LEAST JJ.
C   D      - MATRIX OF SPECIES DIFFUSION COEFFICIENTS AT THE MESH
C            MIDPOINTS.  IF LVARMC=.TRUE. THESE ARE COMPUTED EACH
C            TIME THE FUNCTION IS CALLED, OTHERWISE THE STORED VALUES
C            ARE USED.
C              CGS UNITS - CM**2/SEC
C              DIMENSION D(KK,*) EXACTLY KK FOR THE FIRST DIMENSION
C              AND AT LEAST JJ FOR THE SECOND.
C   TDR    - MATRIX OF SPECIES THERMAL DIFFUSION RATIOS AT THE MESH
C            MIDPOINTS.  IF LVARMC=.TRUE. THESE ARE COMPUTED EACH
C            TIME THE FUNCTION IS CALLED, OTHERWISE THE STORED VALUES
C            ARE USED.
C              CGS UNITS - NONE
C              DIMENSION TDR(KK,*) EXACTLY KK FOR THE FIRST DIMENSION
C              AND AT LEAST JJ FOR THE SECOND.
C   DKJ    - ARRAY OF MULTICOMPONENT DIFFUSION COEFFICIENTS AT THE
C             MESH MIDPOINTS. DIMENSION D(KK,KK,*) EXACTLY KK FOR THE
C             FIRST AND SECOND DIMENSION, AND AT LEAST JJ FOR THE THIRD.
C             CGS UNITS - CM**2/SEC
C   ICKWRK - INTEGER CHEMKIN WORK SPACE.
C              DIMENSIONING - SEE CHEMKIN DOCUMENTATION.
C   RCKWRK  - FLOATING POINT CHEMKIN WORK SPACE.
C              DIMENSIONING - SEE CHEMKIN DOCUMENTATION.
C   IMCWRK - INTEGER TRANSPORT PROPERTY WORK SPACE.
C              DIMENSIONING - SEE TRANSPORT DOCUMENTATION.
C   RMCWRK  - FLOATING POINT TRANSPORT PROPERTY WORK SPACE.
C              DIMENSIONING - SEE TRANSPORT DOCUMENTATION.
C OUTPUT-
C   F      - RESIDUALS OF THE GOVERNING EQUATIONS AS EVALUATED AT
C            S(N,J).  THE RESIDUAL OF THE KK SPECIES EQUATIONS IS IN
C            F(NYS+K,J), THE MASS EQUATION IN F(NM,J), AND THE
C            ENERGY EQUATION IN F(NT,J).  IF INERGY=.FALSE. THEN
C            THE ENERGY EQUATION IS REPLACED BY INTERPOLATING
C            TEMPERATURES FROM THE XGIVEN-TGIVEN TABLE (SEE ABOVE).
C               DIMENSION F(NATJ,*) EXACTLY NATJ FOR THE FIRST
C               DIMENSION AND AT LEAST JJ FOR THE SECOND.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C            timelimit
C
C*****unicos timelimit
C      CALL TREMAIN (SECONDS)
C         IF (SECONDS .LT. 60.0) STOP
C*****END unicos timelimit
C
C             EVALUATE AND STORE THE TRANSPORT COEFFICIENTS
C
      IF (LVARMC) CALL MTRNPR (KK, JJ, NATJ, LENRGY, LMULTI, LTDIF, P,
     1                         X, S, YAV, ICKWRK, RCKWRK, IMCWRK,
     2                         RMCWRK, WT, H, CP, XAV, COND, D, TDR,
     3                         DKJ)
C
C           EVALUATE AND STORE THE DIFFUSION VELOCITIES
C             (in the following call h(*) and cp(*) are used
C              temporarily for scratch space)
C
      CALL MDIFV (KK, JJ, NATJ, LMULTI, LVCOR, LTDIF, X, S, WT, YAV, H,
     1            CP, P, D, TDR, ICKWRK, RCKWRK, YV, DKJ)
C
C                       LEFT BOUNDARY
C
      XMID = 0.5 * (X(1)+X(2))
      TAV = 0.5 * (S(NT,1) + S(NT,2))
      DO 10 K = 1, KK
         YAV(K) = 0.5 * (S(NYS+K,1) + S(NYS+K,2))
10    CONTINUE
      CALL AREA (XMID, AREAP)
      CALL CKRHOY (P, TAV, YAV, ICKWRK, RCKWRK, RHOP)
C
      SUMYK = 0.
      RHOPA = RHOP * AREAP
      DO 30 K = 1, KK-1
         SUMYK = SUMYK + S(NYS+K,1)
         F(NYS+K,1) = EPS(K) - S(NYS+K,1) -
     1                RHOPA * YV(K,1) / S(NM,1)
30    CONTINUE
C
      IF (LVCOR) THEN
         F(NYS+KK,1) = EPS(KK) - S(NYS+KK,1) -
     1                 RHOPA * YV(KK,1) / S(NM,1)
      ELSE
         F(NYS+KK,1) = 1.0 - SUMYK - S(NYS+KK,1)
      ENDIF
C
      F(NT,1) = S(NT,1) - TGIVEN(1)
C
      IF (LBURNR) THEN
         F(NM,1) = S(NM,1) - FLRT
      ELSE
         F(NM,1) = S(NM,2) - S(NM,1)
      ENDIF
C
C                   INTERIOR MESH POINTS
C
      DO 1000 J = 2, JJ-1
C
        TAV = 0.5 * (S(NT,J) + S(NT,J+1))
        DO 100 K = 1, KK
           YAV(K) = 0.5 * (S(NYS+K,J) + S(NYS+K,J+1))
100     CONTINUE
C
        XMDOT = S(NM,J)
        RHOM = RHOP
        CALL CKRHOY (P, TAV, YAV, ICKWRK, RCKWRK, RHOP)
C
        XMID = 0.5 * (X(J)+X(J+1))
        AREAM = AREAP
        CALL AREA (XMID, AREAP)
        AREAJ = 0.5 * (AREAP+AREAM)
C
C             FORM THE CHEMICAL RATE TERMS
C
        CALL CKWYP (P, S(NT,J), S(NY,J), ICKWRK, RCKWRK, WDOT)
        DO 35 K = 1, KK
           WDOT(K) = WDOT(K)*GFAC
   35   CONTINUE
        CALL CKHML (S(NT,J), ICKWRK, RCKWRK, H)
        CALL CKCPBS (S(NT,J), S(NY,J), ICKWRK, RCKWRK, CPB)
        CALL CKCPMS (S(NT,J), ICKWRK, RCKWRK, CP)
C
C               FORM THE MESH DIFFERENCES
C
        DXP =        (X(J+1) - X(J)  )
        DXM =        (X(J)   - X(J-1))
        DXAV = 0.5 * (X(J+1) - X(J-1))
        DXPM =       (X(J+1) - X(J-1))
C
C               FORM THE COEFFICIENTS FOR CENTRAL DIFFERENCES
C
        CENDFM = - DXP / (DXM*DXPM)
        CENDFC =   (DXP-DXM) / (DXP*DXM)
        CENDFP =   DXM / (DXP*DXPM)
C
C               SPECIES CONSERVATION EQUATION
C
        RHOPA = RHOP * AREAP
        RHOMA = RHOM * AREAM
        IF (WNDFAC .EQ. 1.0) THEN
          SUMYK = 0.
          XMDXM = XMDOT / DXM
          DO 200 K = 1, KK-1
             SUMYK = SUMYK + S(NYS+K,J)
             F(NYS+K,J) = XMDXM * ( S(NYS+K,J) - S(NYS+K,J-1) ) +
     1             (RHOPA*YV(K,J) - RHOMA*YV(K,J-1)) / DXAV -
     2               AREAJ * WDOT(K) * WT(K)
200       CONTINUE
          IF (LVCOR) THEN
             F(NYS+KK,J) = XMDXM * (S(NYS+KK,J) - S(NYS+KK,J-1)) +
     1            (RHOPA*YV(KK,J) - RHOMA*YV(KK,J-1))/ DXAV -
     2               AREAJ * WDOT(KK) * WT(KK)
          ELSE
             F(NYS+KK,J) = 1.0 - SUMYK - S(NYS+KK,J)
          ENDIF
        ELSE
          SUMYK = 0.
          DO 210 K = 1, KK-1
             SUMYK = SUMYK + S(NYS+K,J)
             F(NYS+K,J) = XMDOT *
     1                    (CENDFP*S(NYS+K,J+1) + CENDFC*S(NYS+K,J) +
     2                     CENDFM*S(NYS+K,J-1) )  +
     3             (RHOPA*YV(K,J) - RHOMA*YV(K,J-1)) / DXAV -
     4               AREAJ * WDOT(K) * WT(K)
210       CONTINUE
          IF (LVCOR) THEN
             F(NYS+KK,J) = XMDOT *
     1                    (CENDFP*S(NYS+KK,J+1) + CENDFC*S(NYS+KK,J) +
     2                     CENDFM*S(NYS+KK,J-1) )  +
     3           (RHOPA*YV(KK,J) - RHOMA*YV(KK,J-1)) / DXAV -
     4               AREAJ * WDOT(KK) * WT(KK)
          ELSE
             F(NYS+KK,J) = 1.0 - SUMYK - S(NYS+KK,J)
          ENDIF
        ENDIF
C
C               MASS FLOW RATE EQUATION
C
        IF (LBURNR) THEN
           F(NM,J) = S(NM,J) - S(NM,J-1)
        ELSE
           IF (J .GT. JFIXT) THEN
              F(NM,J) = S(NM,J) - S(NM,J-1)
           ELSE
              F(NM,J) = S(NM,J+1) - S(NM,J)
           ENDIF
           IF (J .EQ. JFIXT) F(NM,J) = S(NT,J) - TFIXT
        ENDIF
C
C               ENERGY EQUATION
C
        IF (LENRGY) THEN
C
           SUM = 0.0
           TDOT = 0.0
           DO 400 K = 1, KK
              TDOT = TDOT + WDOT(K)*H(K)
              SUM = SUM + 0.5 * (RHOP*YV(K,J) + RHOM*YV(K,J-1)) *
     1                    CP(K) * (CENDFP*S(NT,J+1) + CENDFC*S(NT,J) +
     2                             CENDFM*S(NT,J-1) )
400        CONTINUE
C
           IF (WNDFAC .EQ. 1.0) THEN
             F(NT,J) = XMDOT * ( S(NT,J) - S(NT,J-1) ) / DXM -
     1            ( COND(J)  *AREAP*(S(NT,J+1)-S(NT,J))/DXP -
     2              COND(J-1)*AREAM*(S(NT,J)-S(NT,J-1))/DXM ) /
     3                                                   (CPB*DXAV) +
     4               AREAJ * (SUM + TDOT) / CPB
           ELSE
             F(NT,J) = XMDOT * (CENDFP*S(NT,J+1) + CENDFC*S(NT,J) +
     1                          CENDFM*S(NT,J-1) ) -
     2            ( COND(J)  *AREAP*(S(NT,J+1)-S(NT,J))/DXP -
     3              COND(J-1)*AREAM*(S(NT,J)-S(NT,J-1))/DXM ) /
     4                                                   (CPB*DXAV) +
     5               AREAJ * (SUM + TDOT) / CPB
           ENDIF
C
        ELSE
           CALL TEMP (NTEMP, X(J), XGIVEN, TGIVEN, TI)
           F(NT, J) = S(NT,J) - TI
        ENDIF
C
1000  CONTINUE
C
C            RIGHT BOUNDARY
C
      SUMYK = 0.
      DO 1100 K = 1, KK-1
         SUMYK = SUMYK + S(NYS+K,JJ)
         F(NYS+K,JJ) = S(NYS+K,JJ) - S(NYS+K,JJ-1)
1100  CONTINUE
C
      IF (LVCOR) THEN
         F(NYS+KK,JJ) = S(NYS+KK,JJ) - S(NYS+KK,JJ-1)
      ELSE
         F(NYS+KK,JJ) = 1.0 - SUMYK - S(NYS+KK,JJ)
      ENDIF
C
      IF (LENRGY) THEN
         F(NT, JJ) = S(NT,JJ) - S(NT,JJ-1)
      ELSE
         CALL TEMP (NTEMP, X(JJ), XGIVEN, TGIVEN, TI)
         F(NT, JJ) = S(NT,JJ) - TI
      ENDIF
C
      F(NM, JJ) = S(NM,JJ) - S(NM,JJ-1)
C
C           ADD THE TIME STEP, IF NEEDED
C
      IF (.NOT. LTIME) RETURN

      DO 2500 J = 2, JJ-1
C
         CALL CKRHOY (P, S(NT,J), S(NY,J), ICKWRK, RCKWRK, RHO)
         CALL AREA (X(J), AREAJ)
C
         IF (LVCOR) THEN
            KK1 = KK
         ELSE
            KK1 = KK - 1
         ENDIF
C
         RHOA = RHO * AREAJ
         DO 2200 K = 1, KK1
            DYDT = (S(NYS+K,J) - SN(NYS+K,J)) / DT
            F(NYS+K,J) = F(NYS+K,J) + RHOA*DYDT
2200     CONTINUE
C
         IF (LENRGY) THEN
            DTDT = (S(NT,J) - SN(NT,J)) / DT
            F(NT,J) = F(NT,J) + RHOA*DTDT
         ENDIF
C
2500  CONTINUE
C
      RETURN
      END
C--------------------------------------------------------------------
      SUBROUTINE MTRNPR (KK, JJ, NATJ, LENRGY, LMULTI, LTDIF,
     1                   P, X, S, YAV, ICKWRK, RCKWRK, IMCWRK, RMCWRK,
     2                   WT, XMF, XMFP, XAV, COND, D, TDR, DKJ)
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /LOCS/ NT, NM, NYS, NY, NTR
C
      DIMENSION X(*), S(NATJ,*), YAV(*), ICKWRK(*), RCKWRK(*),
     1          IMCWRK(*), RMCWRK(*), COND(*), D(KK,*), TDR(KK,*),
     2          DKJ(KK,KK,*), XAV(*), WT(*), XMF(*), XMFP(*)
C
      LOGICAL LENRGY, LTDIF, LMULTI
C
      DATA EPS/1.0E-30/
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   KK     - NUMBER OF CHEMICAL SPECIES.
C   JJ     - NUMBER OF MESH POINTS.
C   NATJ   - NUMBER OF DEPENDENT VARIABLES AT EACH MESH POINT.  ALSO,
C            NATJ MUST BE THE EXACT FIRST DIMENSION OF S(NATJ,*).
C   LENRGY - IF LENRGY=.TRUE. THEN EVALUATE CONDUCTIVITIES AS WELL AS
C            DIFFUSION COEFFICIENTS.
C   LTDIF  - IF LTDIF=.TRUE. THEN EVALUATE THERMAL DIFFUSION RATIOS AS
C            WELL AS DIFFUSION COEFFICIENTS.
C   LMULTI - IF LMULTI=.TRUE. THEN MULTICOMPONENT FORMULAS USED
C            IF LMULTI=.FALSE. THEN MIXTURE-AVERAGED FORMULAS USED
C   P      - PRESSURE.
C              CGS UNITS - DYNES/CM**2
C   X      - THE ARRAY OF MESH POINT LOCATIONS.
C              DIMENSION X(*) AT LEAST JJ
C              CGS UNITS - CM
C   S      - DEPENDENT VARIABLE MATRIX.  THE TEMPERATURES ARE STORED IN
C            T(J)=S(NT,J), THE MASS FRACTIONS ARE IN Y(K,J)=S(NYS+K,J),
C            AND THE FLOW RATES ARE IN FLRT(J)=S(NM,J).
C              DIMENSION S(NATJ,*) EXACTLY NATJ FOR THE FIRST DIMENSION,
C              AND AT LEAST JJ FOR THE SECOND.
C
C WORK AND SCRATCH SPACE
C   YAV    - ARRAY OF MASS FRACTIONS AT MESH MIDPOINTS.
C              DIMENSION Y(*) AT LEAST KK.
C   XAV    - ARRAY OF MOLE FRACTIONS AT MESH MIDPOINTS.
C              DIMENSION X(*) AT LEAST KK
C   ICKWRK - INTEGER CHEMKIN WORK SPACE.
C              DIMENSIONING - SEE CHEMKIN DOCUMENTATION.
C   RCKWRK  - FLOATING POINT CHEMKIN WORK SPACE.
C              DIMENSIONING - SEE CHEMKIN DOCUMENTATION.
C   IMCWRK - INTEGER TRANSPORT PROPERTY WORK SPACE.
C              DIMENSIONING - SEE TRANSPORT DOCUMENTATION.
C   RMCWRK  - FLOATING POINT TRANSPORT PROPERTY WORK SPACE.
C              DIMENSIONING - SEE TRANSPORT DOCUMENTATION.
C OUTPUT-
C   COND   - ARRAY OF CONDUCTIVITIES AT THE MESH MID-POINTS.
C   D      - MATRIX OF DIFFUSION COEFFICIENTS AT THE MESH MID-POINTS.
C              DIMENSION D(KK,*) EXACTLY KK FOR THE FIRST DIMENSION,
C              AND AT LEAST JJ FOR THE SECOND.
C   TDR    - MATRIX OF THERMAL DIFFUSION COEFFICIENTS AT THE MESH
C            MID-POINTS.
C              DIMENSION TDD(KK,*) EXACTLY KK FOR THE FIRST DIMENSION,
C              AND AT LEAST JJ FOR THE SECOND.
C   DKJ    - ARRAY OF MULTICOMPONENT DIFFUSION COEFFICIENTS AT THE
C             MESH MIDPOINTS. DIMENSION D(KK,KK,*) EXACTLY KK FOR THE
C             FIRST AND SECOND DIMENSION, AND AT LEAST JJ FOR THE THIRD.
C             CGS UNITS - CM**2/SEC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IF ( LMULTI) THEN
C-----------------MULTICOMPONENT FORMULAS---------------------
         CALL CKYTX (S(NY,1), ICKWRK, RCKWRK, XMFP)
         DO 200 J = 1, JJ-1
            TAV = 0.5 * (S(NT,J) + S(NT,J+1))
            XDIF = X(J+1) - X(J)
C
C              DIMENSIONAL TEMPERATURE AT THE GRID POINTS
C
            DO 100 K = 1, KK
               YAV(K) = 0.5 * (S(NYS+K,J) + S(NYS+K,J+1))
               XMF(K) = XMFP(K)
100         CONTINUE
            CALL CKYTX (YAV, ICKWRK, RCKWRK, XAV)
            CALL CKYTX (S(NY,J+1), ICKWRK, RCKWRK, XMFP)
            CALL CKMMWX ( XAV, ICKWRK, RCKWRK, WTMAV)
            CALL MCMDIF(P, TAV, XAV, KK, IMCWRK, RMCWRK, DKJ(1,1,J) )
            DO 75 K = 1, KK
               SUMN = 0.0
               DO 50 L = 1, KK
                  SUMN = SUMN + DKJ(K,L,J) *
     1                 (S(NYS+L,J+1) - S(NYS+L,J)) / XDIF
   50          CONTINUE
               DENOM = - (S(NYS+K,J+1) - S(NYS+K,J)) / XDIF
               D(K,J) = (SUMN + EPS) / ( WTMAV * (DENOM + EPS))
   75       CONTINUE
C
C              DETERMINE THE MIXTURE CONDUCTIVITY AND
C              THERMAL DIFFUSION COEFFICIENT AT J
C
            IF (LENRGY .OR. LTDIF) CALL MCMCDT
     1         ( P, TAV, XAV, IMCWRK, RMCWRK, ICKWRK, RCKWRK,
     2           TDR(1,J), COND(J) )
C
200      CONTINUE
C
      ELSE
C-----------------MIXTURE-AVERAGED FORMULAS---------------------
C
         DO 400 J = 1, JJ-1
            TAV = 0.5 * (S(NT,J) + S(NT,J+1))
C
C              DIMENSIONAL TEMPERATURE AT THE GRID POINTS
C
            DO 300 K = 1, KK
               YAV(K) = 0.5 * (S(NYS+K,J) + S(NYS+K,J+1))
300         CONTINUE
            CALL CKYTX (YAV, ICKWRK, RCKWRK, XAV)
            CALL MCADIF(P, TAV, XAV, RMCWRK, D(1,J) )
C
C              DETERMINE THE MIXTURE CONDUCTIVITY AT J
C
            IF (LENRGY) CALL MCACON( TAV, XAV, RMCWRK, COND(J) )
C
            IF (LTDIF) THEN
               CALL MCATDR( TAV, XAV, IMCWRK, RMCWRK, TDR(1,J) )
               CALL CKRHOY (P, TAV, YAV, ICKWRK, RCKWRK, RHOAV)
               CALL CKMMWY (YAV, ICKWRK, RCKWRK, WTM)
               DO 350 K = 1, KK
                  TDR(K,J) = D(K,J) * TDR(K,J) * RHOAV * WT(K)/WTM
  350          CONTINUE
            ENDIF
C
400      CONTINUE
C
      ENDIF
C
      RETURN
      END
C--------------------------------------------------------------------
      SUBROUTINE MDIFV (KK, JJ, NATJ, LMULTI, LVCOR, LTDIF, X, S, WT,
     1                  YAV, XMF, XMFP, P, D, TDR, ICKWRK,
     2                  RCKWRK, YV, DKJ)
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /LOCS/ NT, NM, NYS, NY, NTR
C
      DIMENSION X(*), S(NATJ,*), WT(*), YAV(*), XMF(*), XMFP(*),
     1          D(KK,*), TDR(KK,*), YV(KK,*), DKJ(KK,KK,*)
C
      LOGICAL LTDIF, LMULTI, LVCOR
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   KK     - NUMBER OF CHEMICAL SPECIES.
C   JJ     - NUMBER OF MESH POINTS.
C   NATJ   - NUMBER OF DEPENDENT VARIABLES AT EACH MESH POINT, AND
C            THE EXACT FIRST DIMENSION OF S(NATJ,*), SN(NATJ,*), AND
C            F(NATJ,*).  NATJ=KK+2.
C   LTDIF  - IF LTDIF=.TRUE. THEN EVALUATE THERMAL DIFFUSION RATIOS AS
C            WELL AS DIFFUSION COEFFICIENTS.
C              CGS UNITS - GM/(CM**2-SEC)
C   LMULTI - IF LMULTI=.TRUE. THEN MULTICOMPONENT FORMULAS USED
C            IF LMULTI=.FALSE. THEN MIXTURE-AVERAGED FORMULAS USED
C   LVCOR  - IF LVCOR=.TRUE. THEN USE CORRECTION VELOCITY FORMULISM
C            IF LVCOR=.FALSE. THEN USE 'TRACE' APPROXIMATION, LUMPING
C            ALL TRANSPORT ERRORS INTO THE 'LAST' SPECIES.
C   X      - THE ARRAY OF MESH POINT LOCATIONS.
C              DIMENSION X(*) AT LEAST JJ
C              CGS UNITS - CM
C   S      - DEPENDENT VARIABLE MATRIX.  THE TEMPERATURES ARE STORED IN
C            T(J)=S(NT,J), THE MASS FRACTIONS ARE IN Y(K,J)=S(NYS+K,J),
C            AND THE FLOW RATES ARE IN FLRT(J)=S(NM,J).
C              DIMENSION S(NATJ,*) EXACTLY NATJ FOR THE FIRST DIMENSION,
C              AND AT LEAST JJ FOR THE SECOND.
C   WT     - THE ARRAY OF SPECIES MOLECULAR WEIGHTS.
C              CGS UNITS - GM/MOLE
C              DIMENSION WT(*) AT LEAST KK.
C
C WORK AND SCRATCH SPACE
C   YAV    - ARRAY OF MASS FRACTIONS AT MESH MIDPOINTS.  YAV(J) IS THE
C            MASS FRACTION BETWEEN J AND J+1.
C              DIMENSION YAV(*) AT LEAST KK.
C   XMF    - ARRAY OF MOLE FRACTIONS AT MESH POINT J.
C              DIMENSION XMF(*) AT LEAST KK.
C   XMFP   - ARRAY OF MOLE FRACTIONS AT MESH POINT J+1.
C              DIMENSION XMFP(*) AT LEAST KK.
C   D      - MATRIX OF SPECIES DIFFUSION COEFFICIENTS AT THE MESH
C            MIDPOINTS.  IF LVARMC=.TRUE. THESE ARE COMPUTED EACH
C            TIME THE FUNCTION IS CALLED, OTHERWISE THE STORED VALUES
C            ARE USED.
C              CGS UNITS - CM**2/SEC
C              DIMENSION D(KK,*) EXACTLY KK FOR THE FIRST DIMENSION
C              AND AT LEAST JJ FOR THE SECOND.
C   TDR    - MATRIX OF SPECIES THERMAL DIFFUSION RATIOS AT THE MESH
C            MIDPOINTS.  IF LVARMC=.TRUE. THESE ARE COMPUTED EACH
C            TIME THE FUNCTION IS CALLED, OTHERWISE THE STORED VALUES
C            ARE USED.
C              CGS UNITS - NONE
C              DIMENSION TDR(KK,*) EXACTLY KK FOR THE FIRST DIMENSION
C              AND AT LEAST JJ FOR THE SECOND.
C   DKJ    - ARRAY OF MULTICOMPONENT DIFFUSION COEFFICIENTS AT THE
C             MESH MIDPOINTS. DIMENSION D(KK,KK,*) EXACTLY KK FOR THE
C             FIRST AND SECOND DIMENSION, AND AT LEAST JJ FOR THE THIRD.
C             CGS UNITS - CM**2/SEC
C   ICKWRK - INTEGER CHEMKIN WORK SPACE.
C              DIMENSIONING - SEE CHEMKIN DOCUMENTATION.
C   RCKWRK  - FLOATING POINT CHEMKIN WORK SPACE.
C              DIMENSIONING - SEE CHEMKIN DOCUMENTATION.
C OUTPUT-
C   YV     - MATRIX OF MASS FRACTIONS TIMES DIFFUSION VELOCITIES AT THE
C            MESH MIDPOINTS.  YV(K,J) IS THE FLUX OF KTH SPECIES BETWEEN
C            J AND J+1.
C              DIMENSION YV(KK,*) EXACTLY KK FOR THE FIRST DIMENSION
C              AND AT LEAST JJ FOR THE SECOND.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      CALL CKYTX (S(NY,1), ICKWRK, RCKWRK, XMFP)
C
C           LOOP OVER ALL MESH POINTS, COMPUTING THE DIFFUSION
C         VELOCITY AT THE MID POINTS.  THE INDEXING IS SUCH THAT
C         YV(K,J) IS THE DIFFUSION VELOCITY OF THE KTH SPECIES
C         MIDWAY BETWEEN NODES J AND J+1.
C
      DO 1000 J = 1, JJ-1
C
         TAV = 0.5 * (S(NT,J) + S(NT,J+1))
C
         DO 300 K = 1, KK
            YAV(K) = 0.5 * (S(NYS+K,J) + S(NYS+K,J+1))
300      CONTINUE
C
         DO 400 K = 1, KK
            XMF(K) = XMFP(K)
400      CONTINUE
         CALL CKMMWY (YAV, ICKWRK, RCKWRK, WTM)
         CALL CKRHOY ( P, TAV, YAV, ICKWRK, RCKWRK, RHOAV)
         CALL CKYTX (S(NY,J+1), ICKWRK, RCKWRK, XMFP)
C
         XDIF = X(J+1) - X(J)
         IF ( LMULTI ) THEN
C           EVALUATE THE MULTICOMPONENT DIFFUSION VELOCITY DIRECTLY,
C           RATHER THAN USE THE MIXTURE-AVERAGED FORM FOR D(K,J)
            DO 475 K = 1, KK
               SUM = 0.0
               DO 450 L = 1, KK
                  SUM = SUM + WT(L) * DKJ(K,L,J) *
     1                     (XMFP(L)-XMF(L)) / XDIF
  450          CONTINUE
               YV(K,J) = (WT(K)/WTM**2) * SUM
  475       CONTINUE
         ELSE
C           USE MIXTURE-AVERAGED FORM FOR FICKIAN DIFFUSION,
C           WHETHER WE ARE USING THE MULTICOMPONENT FORMALISM
C           OR MIXTURE-AVERAGED
            DO 500 K = 1, KK
               YV(K,J) = - D(K,J) * (WT(K)/WTM) *
     1                     (XMFP(K)-XMF(K)) / XDIF
500         CONTINUE
         ENDIF
C
C            ADD THE THERMAL DIFFUSION, IF REQUESTED
C
         IF (LTDIF) THEN
C
            TDIF = S(NT, J+1) - S(NT,J)
            DO 600 K = 1, KK
               YV(K,J) = YV(K,J) -
     1                   (TDR(K,J) / (TAV*RHOAV)) * TDIF/XDIF
600         CONTINUE
C
         ENDIF
C
C           COMPUTE AND ADD THE CORRECTION VELOCITY
C
         IF (LVCOR) THEN
            SUM = 0.0
            DO 700 K = 1, KK
               SUM = SUM + YV(K,J)
700         CONTINUE
C
            VC = - SUM
C
            DO 800 K = 1, KK
               YV(K,J) = YV(K,J) + YAV(K)*VC
800         CONTINUE
         ENDIF
C
1000  CONTINUE
C
      RETURN
      END
C------------------------------------------------------------------
      SUBROUTINE PRINT (LOUT, KK, JJ, NATJ, LMOLE, P, X, S, Y, XMF,
     1                  KSYM, ICKWRK, RCKWRK)
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /LOCS/ NT, NM, NYS, NY, NTR
C
      DIMENSION X(*), S(NATJ, *), Y(*), XMF(*), ICKWRK(*), RCKWRK(*)
      CHARACTER KSYM(*)*(*)
C
      LOGICAL LMOLE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C INPUT-
C   LOUT   - UNIT NUMBER FOR OUTPUT.
C   KK     - NUMBER OF CHEMICAL SPECIES.
C   JJ     - NUMBER OF MESH POINTS.
C   NATJ   - NUMBER OF DEPENDENT VARIABLES AT EACH MESH POINT.  ALSO,
C            NATJ MUST BE THE EXACT FIRST DIMENSION OF S(NATJ,*).
C   LMOLE  - FLAG INDICATING WHETHER OUTPUT IS IN MOLE OR MASS FRACTION.
C              LMOLE=.TRUE. MOLE FRACTION.
C              LMOLE=.FALSE. MASS FRACTION.
C   P      - PRESSURE.
C              CGS UNITS - DYNES/CM**2
C   X      - THE ARRAY OF MESH POINT LOCATIONS.
C              DIMENSION X(*) AT LEAST JJ
C              CGS UNITS - CM
C   S      - DEPENDENT VARIABLE MATRIX.  THE TEMPERATURES ARE STORED IN
C            T(J)=S(NT,J)*SF(NT), THE MASS FLOW RATES ARE STORED IN
C            FLRT(J)=S(NM,J), AND THE MASS FRACTIONS  ARE IN
C            Y(K,J)=S(NYS+K,J)*SF(NYS+K).
C              DIMENSION S(NATJ,*) EXACTLY NATJ FOR THE FIRST DIMENSION,
C              AND AT LEAST JJ FOR THE SECOND.
C WORK AND SCRATCH SPACE
C   Y      - ARRAY OF MASS FRACTIONS.
C              DIMENSION Y(*) AT LEAST KK.
C   XMF    - ARRAY OF MOLE FRACTIONS.
C              DIMENSION XMF(*) AT LEAST KK.
C   KSYM   - ARRAY OF CHARACTER*(*) CHEMKIN SPECIES NAMES.
C              DIMENSION KSYM AT LEAST KK.
C   ICKWRK - INTEGER CHEMKIN WORK SPACE.
C              DIMENSIONING - SEE CHEMKIN DOCUMENTATION.
C   RCKWRK  - FLOATING POINT CHEMKIN WORK SPACE.
C              DIMENSIONING - SEE CHEMKIN DOCUMENTATION.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DATA KPERLN /10/
C
      K2 = KPERLN-3
      LS = K2+1
      K2 = MIN (K2, KK)
C
C        PRINT THE FIRST LINE
C
      WRITE (LOUT, 7070) (KSYM(K), K = 1, K2)
      DO 100 J = 1, JJ
         DO 50 K = 1, KK
            Y(K) = S(NYS+K,J)
50       CONTINUE
         CALL CKRHOY (P, S(NT,J), Y, ICKWRK, RCKWRK, RHO)
         CALL AREA (X(J), A)
         FLMSPD = S(NM, J) / (RHO * A)
         IF (LMOLE) THEN
            CALL CKYTX (S(NY,J), ICKWRK, RCKWRK, XMF)
            DO 60 K = 1, KK
               Y(K) = XMF(K)
60          CONTINUE
         ENDIF
         WRITE (LOUT, 7020) J, X(J), S(NT,J), FLMSPD, RHO,
     1                      (Y(K), K=1,K2)
100   CONTINUE
C
      IF (K2 .EQ. KK) RETURN
C
      DO 200 L = LS, KK, KPERLN
         K2   = MIN (L+KPERLN-1, KK)
         WRITE (LOUT, 7060) (KSYM(K), K=L,K2)
         DO 200 J = 1, JJ
            DO 150 K = 1, KK
               Y(K) = S(NYS+K, J)
150         CONTINUE
            IF (LMOLE) THEN
               CALL CKYTX (Y, ICKWRK, RCKWRK, XMF)
               DO 160 K = 1, KK
                  Y(K) = XMF(K)
160            CONTINUE
            ENDIF
            WRITE (LOUT, 7020) J, X(J), (Y(K), K=L,K2)
200     CONTINUE
C
      RETURN
C
7020  FORMAT (I4, 2X, F7.4, 10(1PE11.3))
7060  FORMAT (/8X,'X',8X,10(A10, 1X))
7070  FORMAT (/8X,'X',8X,'T',12X,'V',10X,'RHO',8X,7(A10,1X))
C
      END
C
C---------------------------------------------------------------
C
      SUBROUTINE PRSAVE (ICKWRK, RCKWRK, CCKWRK, 
     1                   IMCWRK, RMCWRK, LOUT, LSAVE)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION ICKWRK(*), RCKWRK(*), IMCWRK(*), RMCWRK(*)
      CHARACTER*(*) CCKWRK(*)
      CHARACTER*16 ILINK
C
      ILINK = 'CKLINK          '
      WRITE (LSAVE) ILINK
      CALL CKSAVE (LOUT, LSAVE, ICKWRK, RCKWRK, CCKWRK)
C
      ILINK = 'MCLINK          '
      WRITE (LSAVE) ILINK
      CALL MCSAVE (LOUT, LSAVE, IMCWRK, RMCWRK)
C
      RETURN
      END
C
C---------------------------------------------------------------------
C
      SUBROUTINE START (KK, NMAX, NATJ, JJ, LOUT, LMOLE, LUMESH, LBURNR,
     1                  NREAC, NINTM, NPROD, REAC, XINTM, PROD, KR, KI,
     2                  KP, FLRT, NTEMP, XGIVEN, TGIVEN, XSTR, XCEN,
     3                  XEND, WMIX, ICKWRK, RCKWRK, Y, SI, EPS, JFIXT,
     4                  TFIXT, X, S)
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /LOCS/ NT, NM, NYS, NY, NTR
C
      DIMENSION REAC(*), XINTM(*), PROD(*), KR(*), KI(*), KP(*),
     1          ICKWRK(*), RCKWRK(*), Y(*), SI(*), EPS(*),
     2          XGIVEN(*), TGIVEN(*), X(*), S(NATJ, *)
C
      LOGICAL LMOLE, LUMESH, LBURNR
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C INPUT-
C   KK     - NUMBER OF CHEMICAL SPECIES.
C   NMAX   - MAXIMUM NUMBER OF MESH POINTS.
C   NATJ   - NUMBER OF DEPENDENT VARIABLES AT EACH MESH POINT, AND
C            THE EXACT FIRST DIMENSION OF S(NATJ,*), SN(NATJ,*), AND
C            F(NATJ,*).  NATJ=KK+2.
C   JJ     - NUMBER OF MESH POINTS IN THE STARTING MESH.
C   LMOLE  - IF LMOLE=.TRUE. THEN THE INPUT AND OUTPUT IS IN TERMS OF
C            MOLE FRACTIONS.
C   LUMESH - IF LUMESH=.TRUE. THEN A UNIFORM STARTING MESH IS USED.
C   LBURNR - IF LBURNR=.TRUE. THIS IS A BURNER STABILIZED FLAME PROBLEM.
C            IF LBURNR=.FALSE. THIS IS A FREELY PROPAGATING ADIABATIC
C            FLAME.
C   NREAC  - NUMBER OF REACTANT SPECIES SPECIFIED.
C   NINTM  - NUMBER OF INTERMEDIATE SPECIES SPECIFIED.
C   NPROD  - NUMBER OF PRODUCT SPECIES SPECIFIED.
C   REAC   - ARRAY OF REACTANT SPECIES INPUT MOLE (MASS) FRACTIONS.
C            DIMENSION REAC(*) AT LEAST KK.
C   XINTM  - ARRAY OF INTERMEDIATE SPECIES INPUT MOLE (MASS) FRACTIONS.
C            DIMENSION XINTM(*) AT LEAST KK.
C   PROD   - ARRAY OF PRODUCT SPECIES INPUT MOLE (MASS) FRACTIONS.
C            DIMENSION PROD(*) AT LEAST KK.
C   KR     - ARRAY OF CHEMKIN SPECIES INDICIES FOR THE REACTANTS AS
C            SPECIFIED IN THE REAC ARRAY.
C            DIMENSION KP(*) AT LEAST KK.
C   KI     - ARRAY OF CHEMKIN SPECIES INDICIES FOR THE INTERMEDIATES AS
C            SPECIFIED IN THE XINTM ARRAY.
C            DIMENSION KI(*) AT LEAST KK.
C   KP     - ARRAY OF CHEMKIN SPECIES INDICIES FOR THE PRODUCTS AS
C            SPECIFIED IN THE PROD ARRAY.
C            DIMENSION KP(*) AT LEAST KK.
C   FLRT   - THE MASS FLOW RATE IN A BURNER STABILIZED PROBLEM.
C              CGS UNITS - GM/(CM**2-SEC)
C   NTEMP  - THE NUMBER OF TEMPERATURE-POSITION PAIRS GIVEN FOR BURNER
C            STABILIZED PROBLEMS WHERE THE ENERGY EQUATION IS NOT
C            SOLVED, AND THE TEMPERATURES ARE INTERPOLATED FROM FIXED
C            PROFILES.
C   XGIVEN - THE DISTANCES FROM THE BURNER AT WHICH TEMPERATURES ARE
C            SPECIFIED.
C              DIMENSION XGIVEN(*) AT LEAST NTEMP.
C   TGIVEN - THE TEMPERATURES AT THE XGIVEN LOCATIONS.
C              DIMENSION TGIVEN(*) AT LEAST NTEMP.
C   XSTR   - BEGINNING X POSITION FOR THE MESH.
C              CGS UNITS - CM
C   XCEN   - X POSITION WHERE THE INITIAL STARTING ESTIMATES ARE
C            CENTERED.
C              CGS UNITS - CM
C   XEND   - ENDING POSITION FOR THE MESH.
C              CGS UNITS - CM
C   WMIX   - WIDTH OF THE MIXING REGION OVER WHICH THE STARTING
C            ESTIMATES ARE FIT.
C              CGS UNITS - CM
C   TFIXT  - TEMPERATURE THAT IS SPECIFIED TO BE HELD FIXED FOR FREE
C            FLAMES.
C              CGS UNITS - K
C
C WORK AND SCRATCH SPACE-
C   ICKWRK - INTEGER CHEMKIN WORK SPACE.
C   RCKWRK  - FLOATING POINT CHEMKIN WORK SPACE.
C   Y      - ARRAY OF MASS FRACTIONS.
C              DIMENSION Y(*) AT LEAST KK.
C   SI     - ARRAY OF THE SUM OF THE INTERMEDIATES.  USED FOR
C            DECREMENTING THE REACTANTS AND PRODUCTS TO KEEP
C            THE MOLE FRACTIONS SUMMED TO 1.
C              DIMENSION SI(*) AT LEAST JJ.
C
C OUTPUT-
C   EPS    - INLET MASS FLUX FRACTONS.
C              DIMENSION EPS(*) AT LEAST KK.
C   JFIXT  - MESH POINT THAT HAS A FIXED TEMPERATURE IN AN ADIABATIC
C            FLAME PROBLEM.
C   X      - STARTING MESH LOCATIONS.
C              DIMENSION X(*) AT LEAST JJ.
C   S      - STARTING SOLUTION ESTIMATES.  S(N,J) IS THE STARTING
C            ESTIMATE FOR COMPONENT N AT NODE J.
C              DIMENSION S(NATJ,*) EXACTLY NATJ FOR THE FIRST DIMENSION
C              AND AT LEAST JJ FOR THE SECOND.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C        INITIALIZE MASS FLUX FRACTIONS AND MASS FRACTIONS TO 0.
C
      DO 90 K = 1, KK
         EPS(K) = 0.0
 90   CONTINUE
      DO 100 K = 1, KK
         DO 95 J = 1, NMAX
            S(NYS+K, J) = 0.0
 95      CONTINUE
100   CONTINUE
C
      IF (LUMESH) THEN
C
C           SET UNIFORM X MESH COORDINATES
C
         DX = (XEND-XSTR) / (JJ-1)
         DO 200 J = 1,JJ
            X(J) = XSTR + DX*(J-1)
200      CONTINUE
      ENDIF
C
C          FOR FREE FLAMES, ADD THE FIXED TEMPERATURE POINT TO THE MESH
C
      IF (.NOT. LBURNR) THEN
C
         DO 300 N = 1, NTEMP
            NFIXT = N
            IF (TGIVEN(N) .GE. TFIXT) GO TO 350
300      CONTINUE
         WRITE (LOUT,*) ' ERROR...NO USABLE MESH LOCATION FOR TFIXT'
         STOP
350      CONTINUE
C
         IF (TGIVEN(NFIXT) .EQ. TFIXT) THEN
            XFIXT = XGIVEN(NFIXT)
         ELSE
            XFIXT = XGIVEN(NFIXT-1) + (TFIXT-TGIVEN(NFIXT-1)) *
     1              (XGIVEN(NFIXT) - XGIVEN(NFIXT-1)) /
     2              (TGIVEN(NFIXT) - TGIVEN(NFIXT-1))
         ENDIF
C
         DO 400 J = 1, JJ
            IF (XFIXT .EQ. X(J)) THEN
               JFIXT = J
               GO TO 700
            ENDIF
400      CONTINUE
C
         DO 500 J = 1, JJ-1
            IF (XFIXT.GT.X(J) .AND. XFIXT.LT.X(J+1)) JFIXT = J+1
500      CONTINUE
C
         JJ = JJ + 1
         DO 600 J = JJ, JFIXT+1, -1
            X(J) = X(J-1)
600      CONTINUE
C
         X(JFIXT) = XFIXT
C
700      CONTINUE
C
      ENDIF
C
C         SET INTERMEDIATE GAUSSIANS
C
      DO 800 N = 1,NINTM
         GBAS = 0.0
         GM   = XINTM(KI(N))
         GMIX = 0.15*GM + GBAS
         W    = -LOG((GMIX-GBAS)/GM)/(WMIX/2.)**2
         DO 750 J = 1, JJ
            S(NYS+KI(N), J) = GM*EXP(-W*(X(J)-XCEN)**2) + GBAS
750      CONTINUE
800   CONTINUE
C
C            SUM THE INTERMEDIATES AT EACH J
C
      DO 1000 J = 1, JJ
         SI(J) = 0.
         DO 900 N = 1, NINTM
            SI(J) = SI(J) + S(NYS+KI(N), J)
900      CONTINUE
1000  CONTINUE
C
C            SET STARTING SPECIES PROFILES
C
      DO 1400 J = 1, JJ
         CALL LINWMX (WMIX, XCEN, X(J), XRE, XPD)
         FAC = 1.0 - SI(J)
         DO 1100 N = 1, NREAC
            S(NYS+KR(N), J) = (XPD*PROD(KR(N)) + XRE*REAC(KR(N))) * FAC
1100     CONTINUE
         DO 1300 N = 1, NPROD
            DO 1200 L = 1, NREAC
               IF (KP(N) .EQ. KR(L)) GO TO 1300
1200        CONTINUE
            S(NYS+KP(N), J) = (XPD*PROD(KP(N)) + XRE*REAC(KP(N))) * FAC
1300     CONTINUE
1400  CONTINUE
C
C             SET THE MASS FLUX FRACTION BOUNDARY CONDITIONS
C
      IF (LMOLE) THEN
         CALL CKXTY (REAC, ICKWRK, RCKWRK, EPS)
      ELSE
         DO 1500 N = 1, NREAC
            EPS(N) = REAC(N)
1500     CONTINUE
      ENDIF
C
C             CONVERT STARTING ESTIMATES TO MASS FRACTION, IF NEEDED
C
      IF (LMOLE) THEN
         DO 1600 J = 1,JJ
            CALL CKXTY (S(NY, J), ICKWRK, RCKWRK, Y)
            DO 1550 K = 1, KK
               S(NYS+K, J) = Y(K)
1550        CONTINUE
1600     CONTINUE
      ENDIF
C
C        SET THE TEMPERATURE AND FLOW RATE PROFILES
C
      DO 1700 J = 1, JJ
         CALL TEMP (NTEMP, X(J), XGIVEN, TGIVEN, S(NT,J))
         S(NM,J) = FLRT
1700  CONTINUE
      IF (.NOT. LBURNR) S(NT,JFIXT)=TFIXT
C
      RETURN
      END
C
C------------------------------------------------------------------
C
      SUBROUTINE LINWMX (WMIX, XCEN, XNODE, XRE, XPD)
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      IF (XNODE .LE. (XCEN-WMIX/2.) ) THEN
         XPD = 0.0
      ELSE
         IF (XNODE .LT. (XCEN+WMIX/2.) ) THEN
            XPD = (1.0/WMIX)*(XNODE-XCEN) + 0.5
         ELSE
            XPD = 1.0
         ENDIF
      ENDIF
C
      XRE = 1.0 - XPD
C
      RETURN
      END
C
C-------------------------------------------------------------------
C
      SUBROUTINE RESTRT (KK, NMAX, NATJ, JJ, LOUT, LMOLE, LUMESH,
     1                   LBURNR, LUSTGV, NREAC, NINTM, NPROD, REAC,
     2                   XINTM, PROD, KR, KI, KP, FLRT, NTEMP, XGIVEN,
     3                   TGIVEN, XSTR, XCEN, XEND, WMIX, ICKWRK,
     4                   RCKWRK, Y, SI, EPS, JFIXT, TFIXT, X, S)
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /LOCS/ NT, NM, NYS, NY, NTR
C
      DIMENSION REAC(*), XINTM(*), PROD(*), KR(*), KI(*), KP(*),
     1          ICKWRK(*), RCKWRK(*), Y(*), SI(*), EPS(*),
     2          XGIVEN(*), TGIVEN(*), X(*), S(NATJ, *)
C
      LOGICAL LMOLE, LUMESH, LBURNR, LUSTGV
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C INPUT-
C   KK     - NUMBER OF CHEMICAL SPECIES.
C   NMAX   - MAXIMUM NUMBER OF MESH POINTS.
C   NATJ   - NUMBER OF DEPENDENT VARIABLES AT EACH MESH POINT, AND
C            THE EXACT FIRST DIMENSION OF S(NATJ,*), SN(NATJ,*), AND
C            F(NATJ,*).  NATJ=KK+2.
C   JJ     - NUMBER OF MESH POINTS IN THE STARTING MESH.
C   LOUT   - UNIT NUMBER FOR PRINTED OUTPUT.
C   LMOLE  - IF LMOLE=.TRUE. THEN THE INPUT AND OUTPUT IS IN TERMS OF
C            MOLE FRACTIONS.
C   LUMESH - IF LUMESH=.TRUE. THEN A UNIFORM STARTING MESH IS USED.
C   LBURNR - IF LBURNR=.TRUE. THIS IS A BURNER STABILIZED FLAME PROBLEM.
C            IF LBURNR=.FALSE. THIS IS A FREELY PROPAGATING ADIABATIC
C            FLAME.
C   NREAC  - NUMBER OF REACTANT SPECIES SPECIFIED.
C   NINTM  - NUMBER OF INTERMEDIATE SPECIES SPECIFIED.
C   NPROD  - NUMBER OF PRODUCT SPECIES SPECIFIED.
C   REAC   - ARRAY OF REACTANT SPECIES INPUT MOLE (MASS) FRACTIONS.
C            DIMENSION REAC(*) AT LEAST KK.
C   XINTM  - ARRAY OF INTERMEDIATE SPECIES INPUT MOLE (MASS) FRACTIONS.
C            DIMENSION XINTM(*) AT LEAST KK.
C   PROD   - ARRAY OF PRODUCT SPECIES INPUT MOLE (MASS) FRACTIONS.
C            DIMENSION PROD(*) AT LEAST KK.
C   KR     - ARRAY OF CHEMKIN SPECIES INDICIES FOR THE REACTANTS AS
C            SPECIFIED IN THE REAC ARRAY.
C            DIMENSION KP(*) AT LEAST KK.
C   KI     - ARRAY OF CHEMKIN SPECIES INDICIES FOR THE INTERMEDIATES AS
C            SPECIFIED IN THE XINTM ARRAY.
C            DIMENSION KI(*) AT LEAST KK.
C   KP     - ARRAY OF CHEMKIN SPECIES INDICIES FOR THE PRODUCTS AS
C            SPECIFIED IN THE PROD ARRAY.
C            DIMENSION KP(*) AT LEAST KK.
C   FLRT   - THE MASS FLOW RATE IN A BURNER STABILIZED PROBLEM.
C              CGS UNITS - GM/(CM**2-SEC)
C   NTEMP  - THE NUMBER OF TEMPERATURE-POSITION PAIRS GIVEN FOR BURNER
C            STABILIZED PROBLEMS WHERE THE ENERGY EQUATION IS NOT
C            SOLVED, AND THE TEMPERATURES ARE INTERPOLATED FROM FIXED
C            PROFILES.
C   XGIVEN - THE DISTANCES FROM THE BURNER AT WHICH TEMPERATURES ARE
C            SPECIFIED.
C              DIMENSION XGIVEN(*) AT LEAST NTEMP.
C   TGIVEN - THE TEMPERATURES AT THE XGIVEN LOCATIONS.
C              DIMENSION TGIVEN(*) AT LEAST NTEMP.
C   XSTR   - BEGINNING X POSITION FOR THE MESH.
C              CGS UNITS - CM
C   XCEN   - X POSITION WHERE THE INITIAL STARTING ESTIMATES ARE
C            CENTERED.
C              CGS UNITS - CM
C   XEND   - ENDING POSITION FOR THE MESH.
C              CGS UNITS - CM
C   WMIX   - WIDTH OF THE MIXING REGION OVER WHICH THE STARTING
C            ESTIMATES ARE FIT.
C              CGS UNITS - CM
C   TFIXT  - TEMPERATURE THAT IS TO BE HELD FIXED FOR FREE FLAMES.
C              CGS UNITS - K
C
C WORK AND SCRATCH SPACE-
C   ICKWRK - INTEGER CHEMKIN WORK SPACE.
C   RCKWRK  - FLOATING POINT CHEMKIN WORK SPACE.
C   Y      - ARRAY OF MASS FRACTIONS.
C              DIMENSION Y(*) AT LEAST KK.
C   SI     - ARRAY OF THE SUM OF THE INTERMEDIATES.  USED FOR
C            DECREMENTING THE REACTANTS AND PRODUCTS TO KEEP
C            THE MOLE FRACTIONS SUMMED TO 1.
C              DIMENSION SI(*) AT LEAST JJ.
C
C OUTPUT-
C   EPS    - INLET MASS FLUX FRACTONS.
C              DIMENSION EPS(*) AT LEAST KK.
C   JFIXT  - MESH POINT THAT HAS A FIXED TEMPERATURE IN AN ADIABATIC
C            FLAME PROBLEM.
C   X      - STARTING MESH LOCATIONS.
C              DIMENSION X(*) AT LEAST JJ.
C   S      - STARTING SOLUTION ESTIMATES.  S(N,J) IS THE STARTING
C            ESTIMATE FOR COMPONENT N AT NODE J.
C              DIMENSION S(NATJ,*) EXACTLY NATJ FOR THE FIRST DIMENSION
C              AND AT LEAST JJ FOR THE SECOND.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C        INITIALIZE MASS FLUX FRACTIONS TO 0.
C
      DO 100 K = 1, KK
         EPS(K) = 0.0
100   CONTINUE
C
C         SHIFT THE SOLUTION TO THE RIGHT IF A NEW XSTR IS .LT. X(1)
C
      IF (XSTR .LT. X(1)) THEN
         JJ = JJ + 1
         IF (JJ .GT. NMAX) THEN
            WRITE (LOUT, *) ' ERROR...NEW XSTR NEEDS TOO MANY POINTS'
            STOP
         ENDIF
         DO 200 I = 2, JJ
            J = JJ + 2 - I
            X(J) = X(J-1)
            DO 150 N = 1, NATJ
               S(N,J) = S(N,J-1)
150         CONTINUE
200      CONTINUE
         X(1) = XSTR
         DO 300 N = 1, NATJ
            S(N,1) = S(N,2)
300      CONTINUE
      ENDIF
C
C          IF A NEW XEND .GT. X(JJ), THEN ADD A POINT AT JJ+1
C          IF A NEW XEND .LE. X(JJ), THEN REDUCE JJ, AND SET X(JJ)=XEND
C
      IF (XEND .GT. (X(JJ)+1.E-4) ) THEN
         JJ = JJ + 1
         IF (JJ .GT. NMAX) THEN
            WRITE (LOUT, *) ' ERROR...NEW XEND NEEDS TOO MANY POINTS'
            STOP
         ENDIF
         X(JJ) = XEND
         DO 400 N = 1, NATJ
            S(N,JJ) = S(N,JJ-1)
400      CONTINUE
      ELSE
C         DO 500 J = 1, JJ
C            IF (XEND .GE. X(J)) THEN
C               X(J) = XEND
C               GO TO 550
C            ENDIF
C500      CONTINUE
C550      CONTINUE
C         JJ = J
      ENDIF
C
C             SET THE MASS FLUX FRACTION BOUNDARY CONDITIONS
C
      IF (LMOLE) THEN
         CALL CKXTY (REAC, ICKWRK, RCKWRK, EPS)
      ELSE
         DO 1000 N = 1, NREAC
            EPS(N) = REAC(N)
1000     CONTINUE
      ENDIF
C
C        IF A NOT A BURNER, OVERWRITE FLRT FROM THE RESTART FILE.
C           OTHERWISE, SET FLRT FROM THE NEW KEYWORD.
C
      IF (.NOT. LBURNR) FLRT = S(NM,1)
C
C        SET THE NEW FLOW RATE PROFILE
C
      DO 1200 J = 1, JJ
         S(NM,J) = FLRT
1200  CONTINUE
C
C         SET XGIVEN AND TGIVEN TO THE OLD SOLUTION
C
      IF (.NOT. LUSTGV) THEN
         NTEMP = JJ
         DO 1250 N = 1, NTEMP
            XGIVEN(N) = X(N)
            TGIVEN(N) = S(NT,N)
1250     CONTINUE
      ENDIF
C
C        FOR A FREE FLAME SET THE MESH POINT FOR THE FIXED TEMPERATURE.
C        SINCE ON THE RESTART FILE, TFIXT MAY NOT EXACTLY EQUAL S(NT,J)
C        WE CHECK FOR A POINT WITHIN 2K OF TFIXT, AND CHANGE TFIXT IF
C        NECESSARY.
C
      IF (.NOT. LBURNR) THEN
         DO 1300 J = 1, JJ
            IF ( ABS(S(NT,J)-TFIXT) .LT. 2.0) THEN
               JFIXT = J
               TFIXT = S(NT,J)
            ENDIF
1300     CONTINUE
      ENDIF
C
      RETURN
      END
C
C--------------------------------------------------------------------
C
      SUBROUTINE RDKEY (KK, NMAX, NATJ, LIN, LOUT, KSYM, PATM, LBURNR,
     1                  LTIME, LTIME2, LMOLE, LUSTGV, LENRGY, LMULTI,
     2                  LVCOR, LTDIF, LUMESH, LRSTRT, LCNTUE, MFILE,
     3                  IPRNT, UFAC, DFAC, DTMAX, DTMIN, LASEN, LHSEN,
     4                  LRSEN, LESEN, LMSEN, LPSEN, LCSEN, LDSEN, P,
     5                  NPTS, NTOT, NADP, X, SCAL, NREAC, NINTM,
     6                  NPROD, REAC, XINTM, PROD, KR, KI, KP, XSTR,
     7                  XCEN, XEND, WMIX, FLRT, GRAD, CURV, SFLR,
     8                  NTEMP, XX, TT, TFIXT, ATOL, RTOL, ATIM, RTIM,
     9                  NJAC, ITJAC, NINIT, NUMDT, IRETIR, DT,
     *                  NUMDT2, DT2, WNDFAC, GFAC, SPOS, N1CALL)
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /LOCS/ NT, NM, NYS, NY, NTR
C
      DIMENSION SCAL(*), REAC(*), XINTM(*), PROD(*), KR(*),
     1          KI(*), KP(*), XX(*), TT(*), X(*), VALUE(5)
C
      LOGICAL LBURNR, LTIME, LTIME2, LMOLE, LENRGY, LTDIF, LUMESH,
     1        LRSTRT, LCNTUE, LUSTGV, LASEN, LRSEN, LESEN, LMSEN, LPSEN,
     2        LCSEN, LDSEN, NEC(8), NOPT(5), CNTNUD, LFIRST, IERR, KERR,
     3        LMULTI, LVCOR, LHSEN
C
      CHARACTER KEYWRD*4, LINE*80, KSYM(*)*(*)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   KK     - NUMBER OF CHEMICAL SPECIES.
C   NMAX   - TOTAL NUMBER OF GRID POINTS THAT CAN BE USED.  NMAX IS
C            USED FOR DIMENSIONING AND ALLOWING FOR DYNAMIC STORAGE
C            ALLOCATION, IT CAN ONLY BE CHANGED IN THE DRIVER CODE.
C   LIN    - UNIT FOR READING KEYWORD INPUT.
C   LOUT   - UNIT FOR PRINTED OUTPUT.
C   KSYM   - ARRAY OF CHARACTER*(*) CHEMKIN SPECIES NAMES.
C              DIMENSION KSYM LEAST KK.
C   PATM   - PRESSURE OF ONE ATMOSPHERE.
C              CGS UNITS - DYNES/CM**2
C
C OUTPUT-
C   LBURNR - IF LBURNR=.TRUE. IT IS A BURNER STABILIZED FLAME PROBLEM.
C            IF LBURNR=.FALSE. THIS IS A FREELY PROPAGATING ADIABATIC
C            FLAME.
C   LTIME  - IF LTIME=.TRUE. THEN TIME STEPPING IS USED.
C            IF LTIME=.FALSE. THEN ONLY NEWTONS METHOD IS USED.
C   LTIME2 - IF LTIME2=.TRUE. THEN TIME STEPPING USED ON EACH NEW MESH.
C            IF LTIME2=.FALSE. THEN ONLY NEWTONS METHOD IS USED.
C   LMOLE  - FLAG INDICATING WHETHER INPUT AND OUTPUT IS IN MOLE OR
C            MASS FRACTION.
C              LMOLE=.TRUE. THEN MOLE FRACTION.
C              LMOLE=.FALSE. THEN MASS FRACTION.
C   LUSTGV - ON A RESTART, IF LUSTGV=.TRUE. THEN THE GIVEN
C            TEMPERATURE PROFILE WILL BE USED RATHER THAN THE ONE
C            ON THE RESTART FILE.
C   LENRGY - IF LENRGY=.TRUE. THEN THE ENERGY EQUATION IS SOLVED.
C            IF LENRGY=.FALSE. THEN A SPECIFIED TEMPERATURE PROFILE IS
C            USED.  THE SPECIFIED PROFILE IS IN (XX(I),TT(I)) PAIRS.
C   LTDIF  - IF LTDIF=.TRUE. THEN INCLUDE THERMAL DIFFUSION.
C            IF LTDIF=.FLASE. THEN NEGLECT THERMAL DIFFUSION.
C   LMULTI - IF LMULTI=.TRUE. THEN MULTICOMPONENT FORMULAS USED
C            IF LMULTI=.FALSE. THEN MIXTURE-AVERAGED FORMULAS USED
C   LVCOR  - IF LVCOR=.TRUE. THEN USE CORRECTION VELOCITY FORMULISM
C            IF LVCOR=.FALSE. THEN USE 'TRACE' APPROXIMATION, LUMPING
C            ALL TRANSPORT ERRORS INTO THE 'LAST' SPECIES.
C   LUMESH - IF LUMESH=.TRUE. THEN START ON A UNIFORM MESH.
C            IF LUMESH=.FALSE. THEN START ON A SPECIFIED NONUNIFORM
C            MESH AS GIVEN IN X, AS SPECIFIED WITH THE GRID KEYWORD.
C   LRSTRT - IF LRSTRT=.TRUE. THEN START FROM A PREVIOUS PROFILE.
C            IF LRSTRT=.FALSE. THEN START FRESH.
C   LCNTUE - IF LCNTUE=.TRUE. THEN A CONTINUATION PROBLEM WILL FOLLOW.
C            IF LCNTUE=.FLASE. THIS IS THE ONLY PROBLEM FOR THE RUN.
C   IPRNT  - FLAG TO SPECIFY THE AMOUNT OF PRINTING.
C              IPRNT = 0, PRINT ONLY THE FINAL SOLUTION.
C              IPRNT = 1, PRINT THE NORMS AFTER EACH ITERATION, AND
C                         THE CONVERGED SOLUTIONS ON EACH MESH.
C              IPRNT = 2, PRINT THE FULL SOLUTION AFTER EACH ITERATION.
C   LASEN  - IF LASEN=.TRUE. COMPUTE SENSITIVITIES FOR ALL REACTIONS.
C            IF LASEN=.FALSE. DO NOT COMPUTE REACTION SENSITIVITIES.
C   LHSEN  - IF LHSEN=.TRUE. COMPUTE H SENSITIVITY FOR ALL SPECIES.
C            IF LASEN=.FALSE. DO NOT COMPUTE H SENSITIVITY.
C   LRSEN  - IF LRSEN=.TRUE. COMPUTE SENSITIVITIES FOR A SET OF
C            SPECIFIED REACTIONS.
C            IF LRSEN=.FALSE. DO NOT COMPUTE REACTION SENSITIVITIES.
C   LESEN  - IF LESEN=.TRUE. COMPUTE SENSITIVITIES FOR THE INLET
C            MASS FLUX FRACTIONS.
C            IF LESEN=.FALSE. DO NOT COMPUTE INLET SENSITIVITIES.
C   LMSEN  - IF LMSEN=.TRUE. COMPUTE SENSITIVITIES FOR THE MASS
C            FLOW RATE.
C            IF LMSEN=.FALSE. DO NOT COMPUTE MASS FLOW SENSITIVITIES.
C   LPSEN  - IF LPSEN=.TRUE. COMPUTE SENSITIVITIES FOR THE PRESSURE.
C            IF LPSEN=.FALSE. DO NOT COMPUTE PRESSURE SENSITIVITIES.
C   LCSEN  - IF LCSEN=.TRUE. COMPUTE SENSITIVITIES FOR THE
C            THERMAL CONDUCTIVITIES.
C            IF LCSEN=.FALSE. DO NOT COMPUTE CONDUCTIVITY
C            SENSITIVITIES.
C   P      - PRESSURE.
C              CGS UNITS - DYNES/CM**2
C   NPTS   - NUMBER OF MESH POINTS IN THE STARTING MESH, X(*).
C   NTOT   - MAXIMUM NUMBER OF GRID POINTS ALLOWED FOR THIS PROBLEM.
C   NADP   - MAXIMUM NUMBER OF MESH POINTS THAT CAN BE ADAPTIVELY
C            ADDED DURING ANY ADAPTION STEP.
C   X      - THE ARRAY OF MESH POINT LOCATIONS.
C              DIMENSION X(*) AT LEAST NMAX
C   SCAL   - ARRAY OF SCALE FACTORS USED FOR COMPUTING NORMS OF
C            THE SOLUTION ITERATES.
C              DIMENSION SCAL(*) AT LEAST NATJ.
C   NREAC  - NUMBER OF REACTANT SPECIES SPECIFIED.
C   NINTM  - NUMBER OF INTERMEDIATE SPECIES SPECIFIED.
C   NPROD  - NUMBER OF PRODUCT SPECIES SPECIFIED.
C   REAC   - ARRAY OF REACTANT SPECIES INPUT MOLE (MASS) FRACTIONS.
C            DIMENSION REAC(*) AT LEAST KK.
C   XINTM  - ARRAY OF INTERMEDIATE SPECIES INPUT MOLE (MASS) FRACTIONS.
C            DIMENSION XINTM(*) AT LEAST KK.
C   PROD   - ARRAY OF PRODUCT SPECIES INPUT MOLE (MASS) FRACTIONS.
C            DIMENSION PROD(*) AT LEAST KK.
C   KR     - ARRAY OF CHEMKIN SPECIES INDICIES FOR THE REACTANTS AS
C            SPECIFIED IN THE REAC ARRAY.
C            DIMENSION KP(*) AT LEAST KK.
C   KI     - ARRAY OF CHEMKIN SPECIES INDICIES FOR THE INTERMEDIATES AS
C            SPECIFIED IN THE XINTM ARRAY.
C            DIMENSION KI(*) AT LEAST KK.
C   KP     - ARRAY OF CHEMKIN SPECIES INDICIES FOR THE PRODUCTS AS
C            SPECIFIED IN THE PROD ARRAY.
C            DIMENSION KP(*) AT LEAST KK.
C   XSTR   - BEGINNING X POSITION FOR THE MESH.
C              CGS UNITS - CM
C   XCEN   - X POSITION WHERE THE INITIAL STARTING ESTIMATES ARE
C            CENTERED.
C              CGS UNITS - CM
C   XEND   - ENDING POSITION FOR THE MESH.
C              CGS UNITS - CM
C   WMIX   - WIDTH OF THE MIXING REGION OVER WHICH THE STARTING
C            ESTIMATES ARE FIT.
C              CGS UNITS - CM
C   FLRT   - THE MASS FLOW RATE IN A BURNER STABILIZED PROBLEM.
C              CGS UNITS - GM/(CM**2-SEC)
C   GRAD   - GRADIENT CRITERIA FOR ADAPTIVE MESHING.
C   CURV   - CURVATURE CIRTERIA FOR ADAPTIVE MESHING.
C   NTEMP  - NUMBER OF XX,TT PAIRS FOR THE SPECIFIED TEMPERATURE
C            PROFILE.  SEE XX(*) AND TT(*)
C   XX     - ARRAY OF X LOCATIONS FOR INITIAL TEMPERATURE PROFILES.
C              DIMENSION XX(*) AT LEAST NMAX.
C              CGS UNITS - CM
C   TT     - ARRAY OF INITIAL TEMPERATURES AT XX.
C              DIMENSION TT(*) AT LEAST NMAX.
C              CGS UNITS - K
C   TFIXT  - TEMPERATURE THAT IS TO BE HELD FIXED FOR FREE FLAMES.
C              CGS UNITS - K
C   ATOL   - ABSOLUTE CONVERGENCE CRITERIA FOR THE NEWTON ITERATION.
C   RTOL   - RELATIVE CONVERGENCE CRITERIA FOR THE NEWTON ITERATION.
C   ATIM   - ABSOLUTE CONVERGENCE CRITERIA FOR THE NEWTON ITERATION
C            AS USED FOR THE TIME STEPS.
C   RTIM   - RELATIVE CONVERGENCE CRITERIA FOR THE NEWTON ITERATION
C            AS USED FOR THE TIME STEPS.
C   NUMDT  - NUMBER OF TIME STEPS TO TAKE WHEN DOING A TIME START.
C   NUMDT2 - NUMBER OF TIME STEPS TO TAKE WHEN TIME STEPPING WITH THE
C            ENERGY EQUATION INCLUDED.
C   DT     - SIZE OF THE TIME STEPS.
C              CGS UNITS - SEC
C   DT2    - SIZE OF THE TIME STEPS AFTER THE ENERGY EQUATION IS
C            INCLUDED.
C              CGS UNITS - SEC
C   WNDFAC - IF WINDFAC=1, THEN WINDWARD DIFFERENCING ON CONVECTION.
C            IF WINDFAC=0, THEN CENTRAL DIFFERENCING ON CONVECTION.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
      DATA NEC/8*.FALSE./, NOPT/5*.FALSE./, CNTNUD/.FALSE./
C
C            INITIALIZE VARIABLES
C
      KERR = .FALSE.
      IF (LCNTUE) THEN
C
         LCNTUE = .FALSE.
         CNTNUD = .TRUE.
         LFIRST = .TRUE.
         NP = 0
C
      ELSE
C
         DO 10 K = 1, KK
            SCAL(NYS+K) = 1.
            REAC(K) = 0.
            XINTM(K) = 0.
            PROD(K)  = 0.
10       CONTINUE
C
         SCAL(NT) = 1000.
         SCAL(NM) = 1.0
         NREAC = 0
         NINTM = 0
         NPROD = 0
         NPTS  = 6
         NTOT  = NMAX
         NADP  = NMAX
         GRAD  = 0.1
         CURV  = 0.5
         SFLR = -1.E-3
         ATOL   = 1.0E-9
         RTOL   = 1.0E-4
         ATIM   = 1.0E-9
         RTIM   = 1.0E-4
         NINIT = 0
         NJAC = 20
         NUMDT = 100
         IRETIR=50
         DT    = 1.0E-6
         NUMDT2= 100
         DT2   = 1.0E-6
         UFAC = 2.0
         DFAC = 2.2
         DTMIN = 1.E-10
         DTMAX = 1.0E-04
         ITJAC = 20
         N1CALL=1
         GFAC = 1.0
         SPOS = -1.0
         MFILE = 1
         NTEMP  = 0
         NP = 0
         WNDFAC = 1.0
         XSTR  = 0.0
         IPRNT = 1
         LFIRST = .TRUE.
         LUMESH = .TRUE.
         LUSTGV = .FALSE.
         LCNTUE = .FALSE.
         LRSTRT = .FALSE.
         LTDIF = .FALSE.
         LMULTI = .FALSE.
         LVCOR = .TRUE.
         LTIME = .FALSE.
         LTIME2 = .FALSE.
         LENRGY = .FALSE.
         LASEN = .FALSE.
         LHSEN = .FALSE.
         LESEN = .FALSE.
         LRSEN = .FALSE.
         LMSEN = .FALSE.
         LPSEN = .FALSE.
         LCSEN = .FALSE.
         LDSEN = .FALSE.
       ENDIF
C
C--------------------------------------------------------------
C
C         READ NEXT INPUT LINE
C
      WRITE (LOUT,'(/A/)') '           KEYWORD INPUT '
C
90    CONTINUE
      KEYWRD = ' '
      LINE = ' '
      IERR = .FALSE.
      READ  (LIN,  7000) KEYWRD, LINE
      WRITE (LOUT, 8000) KEYWRD, LINE
      CALL UPCASE (KEYWRD)
C
C               IS THIS A KEYWORD COMMENT?
C
      IF (KEYWRD(1:1) .EQ. '.' .OR. KEYWRD(1:1) .EQ. '/') GO TO 90
      IND = INDEX(LINE,'(')
      IF (IND .GT. 0) LINE(IND:) = ' '
C
C--------------PROBLEM TYPE KEYWORDS--------------------
C
C           BURNER-STABILIZED FLAME PROBLEM
C
      IF (KEYWRD .EQ. 'BURN') THEN
         LBURNR    = .TRUE.
         NEC(8)   = .TRUE.
C
C       ADIABATIC, FREELY PROPAGATING FLAME PROBLEM
C
      ELSEIF (KEYWRD .EQ. 'FREE') THEN
         LBURNR = .FALSE.
         LENRGY = .TRUE.
         NEC(8)   = .TRUE.
C
C         MOLE FRACTION INPUT/OUTPUT
C
      ELSEIF (KEYWRD .EQ. 'MOLE') THEN
         NEC(1)      = .TRUE.
         LMOLE       = .TRUE.
C
C         MASS FRACTION INPUT/OUTPUT
C
      ELSEIF (KEYWRD .EQ. 'MASS') THEN
         NEC(1)      = .TRUE.
         LMOLE       = .FALSE.
C
C         ENERGY EQUATION IS NOT INCLUDED
C
      ELSEIF (KEYWRD .EQ. 'TGIV') THEN
         NEC(2) = .TRUE.
         LENRGY   = .FALSE.
C
C         ENERGY EQUATION IS INCLUDED
C
      ELSEIF (KEYWRD .EQ. 'ENRG') THEN
         NEC(2)      = .TRUE.
         LENRGY = .TRUE.
C
C--------------METHOD OPTIONS KEYWORDS--------------------
C
C
C       ABSOLUTE NEWTON ITERATION CONVERGENCE CRITERIA
C
      ELSEIF (KEYWRD .EQ. 'ATOL') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, ATOL, IERR)
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'NJAC') THEN
C
C       RETIREMENT AGE OF JACOBIAN DURING STEADY-STATE NEWTON
C
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         NJAC = INT (VALUE(1))
         KERR = KERR.OR.IERR
C
C       RELATIVE NEWTON ITERATION CONVERGENCE CRITERIA
C
      ELSEIF (KEYWRD .EQ. 'RTOL') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, RTOL, IERR)
         KERR = KERR.OR.IERR
C
C       ABSOLUTE NEWTON CONVERGENCE CRITERIA FOR TIMESTEPS
C
      ELSEIF (KEYWRD .EQ. 'ATIM') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, ATIM, IERR)
         KERR = KERR.OR.IERR
C
C       RELATIVE NEWTON CONVERGENCE CRITERIA FOR TIMESTEPS
C
      ELSEIF (KEYWRD .EQ. 'RTIM') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, RTIM, IERR)
         KERR = KERR.OR.IERR
C
C       TIME STEP STARTING PROCEDURE
C
      ELSEIF (KEYWRD .EQ. 'TIME') THEN
         LTIME    = .TRUE.
         CALL CKXNUM (LINE, 2, LOUT, NVAL, VALUE, IERR)
         KERR = KERR.OR.IERR
         NUMDT = INT(VALUE(1))
         DT = VALUE(2)
C
      ELSEIF (KEYWRD .EQ. 'ISTP') THEN
C
C        NUMBER OF INITIAL TIME STEPS BEFORE NEWTON
C
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         NINIT = INT (VALUE(1))
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'IRET') THEN
C
C        RETIREMENT AGE OF OLD TIME STEP (DEFAULT = 50)
C
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         IRETIR = INT (VALUE(1))
         KERR = KERR.OR.IERR
C
C        TIME STEPPING, AFTER ADDING THE ENERGY EQUATION
C
      ELSEIF (KEYWRD .EQ. 'TIM2') THEN
         LTIME2    = .TRUE.
         CALL CKXNUM (LINE, 2, LOUT, NVAL, VALUE, IERR)
         KERR = KERR.OR.IERR
         NUMDT2 = INT(VALUE(1))
         DT2 = VALUE(2)
C
C       TIMESTEP INCREASE WHEN TIMESTEP DOES NOT CHANGE SOLUTION
C
      ELSEIF (KEYWRD .EQ. 'UFAC') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, UFAC, IERR)
         KERR = KERR.OR.IERR
C
C       TIMESTEP DECREASE WHEN NEWTON FAILS CONVERGENCE ON TIMESTEP
C
      ELSEIF (KEYWRD .EQ. 'DFAC') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, DFAC, IERR)
         KERR = KERR.OR.IERR
C
C       MINIMUM TIMESTEP
C
      ELSEIF (KEYWRD .EQ. 'DTMN') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, DTMIN, IERR)
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'DTMX') THEN
C
C       MAXIMUM TIMESTEP
C
         CALL CKXNUM (LINE, 1, LOUT, NVAL, DTMAX, IERR)
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'TJAC') THEN
C
C       RETIREMENT AGE OF JACOBIAN DURING TIM STEPPING
C
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         ITJAC = INT (VALUE(1))
         KERR = KERR.OR.IERR
C
C--------------GRID PARAMETER KEYWORDS--------------------
C
C
C         NUMBER OF INITIAL MESH POINTS
C                          (THIS IS OVERWRITTEN 'GRID' INPUT)
C
      ELSEIF (KEYWRD .EQ. 'NPTS') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         KERR = KERR.OR.IERR
         NPTS   = INT(VALUE(1))
C
C       INITIAL MESH
C
      ELSEIF (KEYWRD .EQ. 'GRID') THEN
         LUMESH    = .FALSE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         IF (IERR .OR. NP+1.GT.NMAX) THEN
            KERR = .TRUE.
         ELSE
            NP = NP+1
            X(NP) = VALUE(1)
         ENDIF
C
C         GRADIENT MESH ADAPTION PARAMETER
C
      ELSEIF (KEYWRD .EQ. 'GRAD') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, GRAD, IERR)
         KERR = KERR.OR.IERR
C
C         CURVATURE MESH ADAPTION PARAMETER
C
      ELSEIF (KEYWRD .EQ. 'CURV') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, CURV, IERR)
         KERR = KERR.OR.IERR
C
C            POINT FOR LEFT BOUNDARY CONDITION
C
      ELSEIF (KEYWRD .EQ. 'XSTR') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, XSTR, IERR)
         KERR = KERR.OR.IERR
C
C             CENTER OF MIXING REGION
C
      ELSEIF (KEYWRD .EQ. 'XCEN') THEN
         NOPT(2) = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, XCEN, IERR)
         KERR = KERR.OR.IERR
C
C            DISTANCE AT WHICH END BOUNDARY CONDITION IS APPLIED
C
      ELSEIF (KEYWRD .EQ. 'XEND') THEN
         NEC(5)  = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, XEND, IERR)
         KERR = KERR.OR.IERR
C
C              WIDTH OF MIXING ZONE
C
      ELSEIF (KEYWRD .EQ. 'WMIX') THEN
         NOPT(1) = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, WMIX, IERR)
         KERR = KERR.OR.IERR
C
C              WINDWARD DIFFERENCING
C
      ELSEIF (KEYWRD .EQ. 'WDIF') THEN
         WNDFAC = 1.0
C
C              CENTRAL DIFFERENCING
C
      ELSEIF (KEYWRD .EQ. 'CDIF') THEN
         WNDFAC = 0.0
C
C        FLOOR VALUE FOR THE SPECIES BOUNDS
C
      ELSEIF (KEYWRD .EQ. 'SFLR') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, SFLR, IERR)
         KERR = KERR.OR.IERR
C
C--------------FLAME DEFINITION KEYWORDS--------------------
C
C         PRESSURE
C
      ELSEIF (KEYWRD .EQ. 'PRES') THEN
         NEC(3)  = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, P, IERR)
         KERR = KERR.OR.IERR
         P = P*PATM
C
      ELSEIF (KEYWRD .EQ. 'FLRT') THEN
C
C        MASS FLOW RATE (gm/sec)
C
         NEC(4)  = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, FLRT, IERR)
         KERR = KERR.OR.IERR
C
C         REACTANT
C
      ELSEIF (KEYWRD .EQ. 'REAC') THEN
C
         IF (LFIRST) THEN
            LFIRST = .FALSE.
            NREAC = 0
            DO 1100 K = 1, KK
               REAC(K) = 0.
1100        CONTINUE
         ENDIF
         CALL CKSNUM (LINE, 1, LOUT, KSYM, KK, KSPEC, NVAL,
     1                VALUE, IERR)
         IF (IERR .OR. NREAC+1.GT.KK) THEN
            WRITE (LOUT,'(A)')
     1      ' ERROR READING DATA FOR KEYWORD '//KEYWRD
            KERR = .TRUE.
         ELSE
            NREAC       = NREAC+1
            KR(NREAC)   = KSPEC
            REAC(KSPEC) = VALUE(1)
         ENDIF
C
C         INTERMEDIATE
C
      ELSEIF (KEYWRD .EQ. 'INTM') THEN
         CALL CKSNUM (LINE, 1, LOUT, KSYM, KK, KSPEC, NVAL,
     1                VALUE, IERR)
         IF (IERR .OR. NINTM+1.GT.KK) THEN
            WRITE (LOUT,'(A)')
     1      ' ERROR READING DATA FOR KEYWORD '//KEYWRD
            KERR = .TRUE.
         ELSE
            NINTM       = NINTM + 1
            KI(NINTM)   = KSPEC
            XINTM(KSPEC) = VALUE(1)
         ENDIF
C
C         PRODUCT
C
      ELSEIF (KEYWRD .EQ. 'PROD') THEN
         CALL CKSNUM (LINE, 1, LOUT, KSYM, KK, KSPEC, NVAL,
     1                VALUE, IERR)
         IF (IERR .OR. NPROD+1.GT.KK) THEN
            WRITE (LOUT,'(A)')
     1      ' ERROR READING DATA FOR KEYWORD '//KEYWRD
            KERR = .TRUE.
         ELSE
            NPROD = NPROD + 1
            KP(NPROD) = KSPEC
            PROD(KSPEC) = VALUE(1)
         ENDIF
C
C         READ SPECIFIED TEMPERATURE PROFILE (X,T) PAIRS
C
      ELSEIF (KEYWRD .EQ. 'TEMP') THEN
         CALL CKXNUM (LINE, 2, LOUT, NVAL, VALUE, IERR)
         KERR = KERR.OR.IERR
         IF (NTEMP+1 .GT. NMAX) THEN
               WRITE (LOUT, *)
     1         ' ERROR... THE PROBLEM IS ONLY DIMENSIONED FOR ', NMAX,
     2         ' (X,T) PAIRS'
               KERR = .TRUE.
         ELSE
            NTEMP = NTEMP + 1
            XX(NTEMP) = VALUE(1)
            TT(NTEMP) = VALUE(2)
         ENDIF
C
C          ON A RESTART USE GIVEN TEMPERATURE PROFILE, NOT THE ONE ON
C               THE RESTART FILE
C
      ELSEIF (KEYWRD .EQ. 'USTG') THEN
         LUSTGV = .TRUE.
C
C       TEMPERATURE WHICH IS TO BE HELD FIXED FOR A FREE FLAME
C
      ELSEIF (KEYWRD .EQ. 'TFIX') THEN
         NEC(6)   = .TRUE.
         CALL CKXNUM (LINE, 1, LOUT, NVAL, TFIXT, IERR)
         KERR = KERR.OR.IERR
C
C--------------TRANSPORT OPTIONS KEYWORDS--------------------
C
C       MULTICOMPONENT FORMULAS USED
C
      ELSEIF (KEYWRD .EQ. 'MULT') THEN
         LMULTI = .TRUE.
C
C       MIXTURE-AVERAGED FORMULAS USED
C
      ELSEIF (KEYWRD .EQ. 'MIX') THEN
         LMULTI = .FALSE.
C
C       THERMAL DIFFUSION INCLUDED
C
      ELSEIF (KEYWRD .EQ. 'TDIF') THEN
         LTDIF    = .TRUE.
C
C       USE CORRECTION VELOCITY FORMALISM
C
      ELSEIF (KEYWRD .EQ. 'VCOR') THEN
         LVCOR = .TRUE.
C
C       USE "TRACE" APPROXIMATION , LUMP ALL TRANSPORT ERRORS
C        INTO THE "LAST" SPECIES
C
      ELSEIF (KEYWRD .EQ. 'TRCE') THEN
         LVCOR = .FALSE.
C
C--------------SENSITIVITY KEYWORDS--------------------
C
C        ALL REACTION SENSITIVITY
C
      ELSEIF (KEYWRD .EQ. 'ASEN') THEN
         LASEN = .TRUE.
C
C        SENSITIVITY TO HEATS OF FORMATION
C
      ELSEIF (KEYWRD .EQ. 'HSEN') THEN
         LHSEN = .TRUE.
C
C--------------PRINTING AND RESTARTING KEYWORDS--------------------
C
C          PRINT CONTROL
C
      ELSEIF (KEYWRD .EQ. 'PRNT') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         KERR = KERR.OR.IERR
         IPRNT    = INT(VALUE(1))
C
C          RESTART SKIPS
C
      ELSEIF (KEYWRD .EQ. 'SKIP') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         KERR = KERR.OR.IERR
         MFILE = INT(VALUE(1)) + 1
C
C          RESTART CHECK
C
      ELSEIF (KEYWRD .EQ. 'RSTR') THEN
         LRSTRT = .TRUE.
C
C          CONTINUATION FLAG
C
      ELSEIF (KEYWRD .EQ. 'CNTN') THEN
         LCNTUE   = .TRUE.
C
C        LAST CARD
C
      ELSEIF (KEYWRD .EQ. 'END ') THEN
         GO TO 6000
C
C-------------THE KEYWORDS AFTER HERE ARE NOT OPERATIONAL------------
C
C         SPECIES SCALE FACTOR
C
      ELSEIF (KEYWRD .EQ. 'SCAL') THEN
         CALL CKSNUM (LINE, 1, LOUT, KSYM, KK, KSPEC, NVAL,
     1                VALUE, IERR)
         IF (IERR) THEN
            WRITE (LOUT,'(A)')
     1      ' ERROR READING DATA FOR KEYWORD '//KEYWRD
            KERR = .TRUE.
         ENDIF
         SCAL(NYS+KSPEC) = VALUE(1)
         WRITE (LOUT,*)
     1                ' WARNING...SPECIES SCALING IS NOT IMPLEMENTED.'
C
C         SCALE TEMPERATURE
C
      ELSEIF (KEYWRD .EQ. 'SCLT') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, SCAL(NT), IERR)
         KERR = KERR.OR.IERR
         WRITE (LOUT,*) ' WARNING...TEMP SCALING IS NOT IMPLEMENTED.'
C
C        MASS FLUX FRACTION SENSITIVITY
C
      ELSEIF (KEYWRD .EQ. 'ESEN') THEN
         LESEN      = .TRUE.
         WRITE (LOUT,*)
     1              ' WARNING...FLUX SENSITIVITY IS NOT IMPLEMENTED.'
C
C        MASS FLOW RATE SENSITIVITY
C
      ELSEIF (KEYWRD .EQ. 'MSEN') THEN
         LMSEN      = .TRUE.
         WRITE (LOUT,*)
     1             ' WARNING...FLOW SENSITIVITY IS NOT IMPLEMENTED.'
C
C       PRESSURE SENSITVITY
C
      ELSEIF (KEYWRD .EQ. 'PSEN') THEN
         LPSEN      = .TRUE.
         WRITE (LOUT,*)
     1               ' WARNING...PRES SENSITIVITY IS NOT IMPLEMENTED.'
C
C       CONDUCTIVITY SENSITIVITY
C
      ELSEIF (KEYWRD .EQ. 'CSEN') THEN
         LCSEN      = .TRUE.
         WRITE (LOUT,*)
     1            ' WARNING...COND SENSITIVITY IS NOT IMPLEMENTED.'
C
C        DIFFUSION COEFFICIENT SENSITIVITY
C
      ELSEIF (KEYWRD .EQ. 'DSEN') THEN
         LDSEN      = .TRUE.
         WRITE (LOUT,*)
     1            ' WARNING...DIFFUS SENSITIVITY IS NOT IMPLEMENTED.'
C
C       TOTAL NUMBER OF POINTS ALLOWED FOR THIS RUN
C
      ELSEIF (KEYWRD .EQ. 'NTOT') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         KERR = KERR.OR.IERR
         NTOT = INT(VALUE(1))
C
C       TOTAL NUMBER OF ADAPTIVE POINTS THAT CAN BE ADDED PER ADAPTION
C
      ELSEIF (KEYWRD .EQ. 'NADP') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         KERR = KERR.OR.IERR
         NADP = INT(VALUE(1))
C
      ELSEIF (KEYWRD .EQ. 'NOFT') THEN
C
C        DO NOT DO THE FIXED TEMPERATURE PROBLEM
C
         N1CALL = 2
C
      ELSEIF (KEYWRD .EQ. 'SPOS') THEN
C
C        CONVERT NEGATIVE SPECIES SOLUTIONS
C
         CALL CKXNUM (LINE, 1, LOUT, NVAL, SPOS, IERR)
         KERR = KERR.OR.IERR
C
      ELSEIF (KEYWRD .EQ. 'GFAC') THEN
C
C        FACTOR FOR GAS-PHASE RATE CONSTANTS
C
         CALL CKXNUM (LINE, 1, LOUT, NVAL, GFAC, IERR)
         KERR = KERR.OR.IERR
C
      ELSE
C
C--------------END OF KEYWORDS--------------------
C
C        TO GET HERE, AN INVALID KEYWORD WAS READ
C
         WRITE (LOUT,*) ' ERROR...ILLEGAL KEYWORD'
         KERR = .TRUE.
      ENDIF
      GO TO 90
C
C         CHECK THE REACTANT AND PRODUCT SUMS
C
6000  CONTINUE
C
      SUMR       = 0.
      SUMP       = 0.
      DO 6100 K = 1, KK
         SUMR = SUMR + REAC(K)
         SUMP = SUMP + PROD(K)
6100  CONTINUE
C
C         NORMALIZE REACTANT AND PRODUCT FRACTIONS
C
      DO 6200 K = 1, KK
         REAC(K) = REAC(K)/SUMR
         PROD(K) = PROD(K)/SUMP
6200  CONTINUE
C
      IF (ABS(SUMR-1.0) .GT. 1.E-6) WRITE (LOUT, *)
     1                ' CAUTION...REACTANT FRACTIONS SUM TO ', SUMR
      IF ((.NOT.CNTNUD) .AND.  ABS(SUMP-1.0).GT.1.E-6)
     1    WRITE (LOUT, *)
     2                ' CAUTION...PRODUCT FRACTIONS SUM TO ',  SUMP
C
C          CHECK FOR NECESSARY INPUT
C
      IF (.NOT. NEC(8) )THEN
         WRITE (LOUT, *) ' ERROR..."BURN" OR "FREE" NOT SPECIFIED'
         KERR = .TRUE.
      ENDIF
C
      IF (.NOT. NEC(1) ) THEN
         WRITE (LOUT, *)
     1              ' ERROR...MUST SPECIFY EITHER "MOLE" OR "MASS" '
         KERR = .TRUE.
      ENDIF
C
      IF (.NOT. LBURNR) LENRGY = .TRUE.
C
      IF ((.NOT. NEC(2)) .AND. (LBURNR)) THEN
           WRITE (LOUT, *)
     1 ' ERROR..."ENRG" OR "TGIV" MUST BE PROVIDED FOR A BURNER FLAME'
         KERR = .TRUE.
      ENDIF
C
      IF (.NOT. NEC(3) ) THEN
         WRITE (LOUT, *) ' ERROR...PRESSURE NOT GIVEN'
         KERR = .TRUE.
      ENDIF
C
      IF (.NOT. NEC(4)) THEN
         WRITE (LOUT, *) ' ERROR...MASS FLOW RATE NOT SPECIFIED'
         KERR = .TRUE.
      ENDIF
C
      IF (.NOT. NEC(5)) THEN
         WRITE (LOUT, *) ' ERROR..."XEND" NOT SPECIFIED'
         KERR = .TRUE.
      ENDIF
C
      IF ((.NOT. NEC(6)) .AND. (.NOT.LBURNR)) THEN
            WRITE (LOUT, *)
     1                 ' ERROR..."TFIX" NOT GIVEN FOR A FREE FLAME'
           KERR = .TRUE.
      ENDIF
C
C           MAKE SURE THE (X,T) PAIRS ARE IN ORDER
C
      DO 5500 N = 2, NTEMP
         IF (XX(N-1) .GE. XX(N)) THEN
            WRITE (LOUT, *)
     1              ' ERROR...SPECIFIED TEMPERATURES ARE OUT OF ORDER'
            KERR = .TRUE.
         ENDIF
5500  CONTINUE
C
C           MAKE SURE THE INITIAL GRID POINTS ARE IN ORDER
C
      IF ((.NOT.CNTNUD) .AND. (.NOT.LUMESH)) THEN
         NPTS = NP
         DO 5550 N = 2, NPTS
            IF (X(N-1) .GE. X(N)) THEN
              WRITE (LOUT, *)
     1                 ' ERROR...INITIAL GRID IS OUT OF ORDER'
              KERR = .TRUE.
            ENDIF
5550     CONTINUE
      ENDIF
C
C         MAKE SURE THE GIVEN TEMPERATURES SPAN THE XEND-XSTR DOMAIN
C
      IF (.NOT.LRSTRT .OR. .NOT.CNTNUD .OR. LUSTGV) THEN
         IF(XX(1).GT.XSTR .OR. XX(NTEMP).LT.XEND) THEN
            WRITE (LOUT,*)
     1      ' ERROR...GIVEN TEMPERATURE PROFILE DOES NOT SPAN XEND-XSTR'
            KERR = .TRUE.
         ENDIF
      ENDIF
C
C         SET OPTIONAL INPUT IF NEEDED
C
      IF (.NOT. NOPT(1)) WMIX = (XEND-XSTR)*0.50
      IF (.NOT. NOPT(2)) XCEN = (XEND-XSTR)*0.35
C
      IF (KERR) STOP
C
C         formats
C
      RETURN
7000  FORMAT (A4, A)
8000  FORMAT (10X, A4, A76)
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE TEMP (NPTS, X, XX, TT, T)
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION XX(NPTS), TT(NPTS)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     THIS SUBROUTINE USES BISECTION TO LINEARLY INTERPOLATE
C     AN ARRAY OF XX,TT PAIRS.  GIVEN AN XX,TT PAIR THIS ROUTINE
C     RETURNS THE INTERPOLATED VALUE OF THE T AT THE POINT X.
C
C INPUT-
C   NPTS   - NUMBER OF XX,TT PAIRS.
C   X      - LOCATION AT WHICH INTERPOLATED T IS DESIRED.
C   XX     - ARRAY OF X POINTS AT WHICH TT ARE GIVEN.
C   TT     - ARRAY OF T VALUES AT THE XX LOCATIONS.
C
C OUTPUT-
C   T     - INTERPOLATED T AT POINT X
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C             check for x outside (1,npts)
C
      IF (X .LE. XX(2)) THEN
         N = 2
         S = (TT(N) - TT(N-1)) / (XX(N) - XX(N-1))
      ELSEIF (X .GE. XX(NPTS-1)) THEN
         N = NPTS-1
         S = (TT(N+1) - TT(N)) / (XX(N+1) - XX(N))
      ELSE
         NLO = 1
         NHI = NPTS
         S   = 0.0
C
C        bisect interval
C
50       CONTINUE
         N = (NLO+NHI)/2
         IF (X .LT. XX(N)) THEN
            IF (X .LT. XX(N-1)) THEN
               NHI = N
               GO TO 50
            ELSEIF (X .EQ. XX(N-1)) THEN
               N = N-1
            ELSE
               S = (TT(N) - TT(N-1)) / (XX(N) - XX(N-1))
            ENDIF
         ELSEIF (X .GT. XX(N)) THEN
            IF (X .GT. XX(N+1)) THEN
               NLO = N
               GO TO 50
            ELSEIF (X .EQ. XX(N+1)) THEN
               N = N + 1
            ELSE
               S = (TT(N+1) - TT(N)) / (XX(N+1) - XX(N))
            ENDIF
         ENDIF
      ENDIF
C
  100 CONTINUE
      T      = TT(N) + S * (X - XX(N))
      RETURN
      END
C--------------------------------------------------------------------
C
      SUBROUTINE REASEN (II, KK, JJ, NATJ, LDA, LBURNR, LENRGY, LMULTI,
     1                   LTDIF, LMOLE, LSAVE, LRCRVR, LOUT, LVARMC,
     2                   LTIME, JFIXT, TFIXT, FLRT, P, WT, EPS, DT,
     3                   NTEMP, XGIVEN, TGIVEN, X, SN, S, WNDFAC, ABS,
     4                   REL, SCRCHK, YV, COND, D, DKJ, TDR, ICKWRK,
     5                   RCKWRK, IMCWRK, RMCWRK, F, FN, A, DS, SSAVE,
     6                   IP, GFAC)
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION WT(*), EPS(*), XGIVEN(NTEMP), TGIVEN(NTEMP), X(*),
     1          SN(NATJ,*), S(NATJ, *), F(NATJ,*), FN(NATJ,*),
     2          SCRCHK(KK,*), YV(KK,*), COND(*), D(KK,*), TDR(KK,*),
     3          DS(*), SSAVE(*), IP(*), ICKWRK(*), RCKWRK(*),
     4          IMCWRK(*), RMCWRK(*), A(LDA, *), DKJ(KK,KK,*)
C
      CHARACTER*16 ISENSI
      COMMON /LOCS/ NT, NM, NYS, NY, NTR
      LOGICAL LBURNR, LENRGY, LTDIF, LVARMC, LTIME, LMOLE, LMULTI,
     1        LVCOR
      DATA ISENSI/'SENSITIVITY     '/
C
C/////////////////////////////////////////////////////////////////
C
      LTIME = .FALSE.
      LVARMC = .TRUE.
      ML = 2*NATJ - 1
      MU = 2*NATJ - 1
C
C          FORM THE JACOBIAN, AND AT THE SAME TIME EVALUATE THE
C               UNPERTURBED FUNCTION FN
C
      CALL JACOB  (KK, JJ, NATJ, LBURNR, LENRGY, LMULTI, LVCOR, LTDIF,
     1             LVARMC, LTIME, JFIXT, TFIXT, FLRT, P, WT, EPS, DT,
     2             NTEMP, XGIVEN, TGIVEN, X, SN, S, WNDFAC, ABS, REL,
     3             SCRCHK, YV, COND, D, DKJ, TDR, ICKWRK, RCKWRK,
     4             IMCWRK, RMCWRK, F, FN, A, DS, SSAVE, GFAC)
C
C          FACTOR THE JACOBIAN
C
C*****precision > double
      CALL DGBFA (A, LDA, NATJ*JJ, ML, MU, IP, INFO)
C*****END precision > double
C
C*****precision > single
C      CALL SGBFA (A, LDA, NATJ*JJ, ML, MU, IP, INFO)
C*****END precision > single
C
      IF (INFO .NE. 0) THEN
         WRITE (LOUT, *)
     1   'ERROR IN FACTORING THE JACOBIAN, INFO = ',INFO
         STOP
      ENDIF
C
C          COMPUTE THE RAW SENSITIVITY COEFFICIENTS WITH RESPECT TO
C             THE RATE CONSTANTS, D(MASS FRACTION)/D(RATE CONSTANT).
C
      WRITE (LSAVE) ISENSI
      DO 1000 I = 1, II
C
         CALL CKRDEX (I, RCKWRK, SAVEP)
         DP    = REL*SAVEP + ABS
         CALL CKRDEX (-I, RCKWRK, SAVEP+DP)
C
         CALL FUN (KK, JJ, NATJ, LBURNR, LENRGY, LMULTI, LVCOR, LTDIF,
     1             LVARMC, LTIME, JFIXT, TFIXT, FLRT, P, WT, EPS, DT,
     2             NTEMP, XGIVEN, TGIVEN, X, SN, S, WNDFAC, SCRCHK(1,1),
     3             YV, SCRCHK(1,2), SCRCHK(1,3), SCRCHK(1,4), COND, D,
     4             DKJ, TDR, ICKWRK, RCKWRK, IMCWRK, RMCWRK, F,
     5             SCRCHK(1,5), GFAC)
C
         CALL CKRDEX (-I, RCKWRK, SAVEP)
         DO 100 J = 1, JJ
            DO 90 N = 1, NATJ
               SN(N,J) =  - (F(N,J)-FN(N,J)) / DP
90          CONTINUE
100      CONTINUE
C
C*****precision > double
         CALL DGBSL (A, LDA, NATJ*JJ, ML, MU, IP, SN, 0)
C*****END precision > double
C
C*****precision > single
C         CALL SGBSL (A, LDA, NATJ*JJ, ML, MU, IP, SN, 0)
C*****END precision > single
C
C          SN(N,J) NOW CONTAINS THE RAW SENSIVITY MATRIX DS(N,J)/DAi
C
C          NORMALIZE THE SENSIVITY COEFFICIENTS
C
         DO 200 J = 1, JJ
            SN(NT,J) = SN(NT,J) * SAVEP / S(NT,J)
200      CONTINUE
C
         IF (LMOLE) THEN
C
            DO 500 J = 1, JJ
               SUM = 0.0
               DO 300 L = 1, KK
                  SUM = SUM + SN(NYS+L,J) / WT(L)
300            CONTINUE
               CALL CKMMWY (S(NY,J), ICKWRK, RCKWRK, WTM)
               WTMSUM = WTM * SUM
               DO 400 K = 1, KK
                  SN(NYS+K,J) = SAVEP * (SN(NYS+K,J) / S(NYS+K,J)
     2                          - WTMSUM)
400            CONTINUE
500         CONTINUE
C
         ELSE
C
            DO 700 J = 1, JJ
               DO 600 K = 1, KK
                  SN(NYS+K,J) = SN(NYS+K,J) * SAVEP / S(NYS+K,J)
600            CONTINUE
700         CONTINUE
C
         ENDIF
C
         DO 800 J = 1, JJ
            SN(NM,J) = SN(NM,J) * SAVEP / S(NM,J)
800      CONTINUE
C
         WRITE (LSAVE)  I, ((SN(N,J),N=1,NATJ), J=1,JJ)
         WRITE (LRCRVR) I, ((SN(N,J),N=1,NATJ), J=1,JJ)
C
1000  CONTINUE
C
      RETURN
      END
C
      SUBROUTINE REHSEN (II, KK, JJ, NATJ, LDA, LBURNR, LENRGY, LMULTI,
     1                   LTDIF, LMOLE, LSAVE, LRCRVR, LOUT, LVARMC,
     2                   LTIME, JFIXT, TFIXT, FLRT, P, WT, EPS, DT,
     3                   NTEMP, XGIVEN, TGIVEN, X, SN, S, WNDFAC, ABS,
     4                   REL, SCRCHK, YV, COND, D, DKJ, TDR, ICKWRK,
     5                   RCKWRK, IMCWRK, RMCWRK, F, FN, A, DS, SSAVE,
     6                   IP, GFAC, A6)
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION WT(*), EPS(*), XGIVEN(NTEMP), TGIVEN(NTEMP), X(*),
     1          SN(NATJ,*), S(NATJ, *), F(NATJ,*), FN(NATJ,*),
     2          SCRCHK(KK,*), YV(KK,*), COND(*), D(KK,*), TDR(KK,*),
     3          DS(*), SSAVE(*), IP(*), ICKWRK(*), RCKWRK(*),
     4          IMCWRK(*), RMCWRK(*), A(LDA, *), DKJ(KK,KK,*), A6(*)
C
      CHARACTER*16 IHSENS
      COMMON /LOCS/ NT, NM, NYS, NY, NTR
      LOGICAL LBURNR, LENRGY, LTDIF, LVARMC, LTIME, LMOLE, LMULTI,
     1        LVCOR
      DATA IHSENS/'HSENSITIVITY    '/
C
C/////////////////////////////////////////////////////////////////
C
      LTIME = .FALSE.
      LVARMC = .TRUE.
      ML = 2*NATJ - 1
      MU = 2*NATJ - 1
C
C          FORM THE JACOBIAN, AND AT THE SAME TIME EVALUATE THE
C               UNPERTURBED FUNCTION FN
C
      CALL JACOB  (KK, JJ, NATJ, LBURNR, LENRGY, LMULTI, LVCOR, LTDIF,
     1             LVARMC, LTIME, JFIXT, TFIXT, FLRT, P, WT, EPS, DT,
     2             NTEMP, XGIVEN, TGIVEN, X, SN, S, WNDFAC, ABS, REL,
     3             SCRCHK, YV, COND, D, DKJ, TDR, ICKWRK, RCKWRK,
     4             IMCWRK, RMCWRK, F, FN, A, DS, SSAVE, GFAC)
C
C          FACTOR THE JACOBIAN
C
C*****precision > double
      CALL DGBFA (A, LDA, NATJ*JJ, ML, MU, IP, INFO)
C*****END precision > double
C
C*****precision > single
C      CALL SGBFA (A, LDA, NATJ*JJ, ML, MU, IP, INFO)
C*****END precision > single
C
      IF (INFO .NE. 0) THEN
         WRITE (LOUT, *)
     1   'ERROR IN FACTORING THE JACOBIAN, INFO = ',INFO
         STOP
      ENDIF
C
C          COMPUTE THE RAW SENSITIVITY COEFFICIENTS WITH RESPECT TO
C             THE RATE CONSTANTS, D(MASS FRACTION)/D(ENTHALPY).
C
      WRITE (LSAVE) IHSENS
      DO 1000 K = 1, KK
C
         CALL CKRHEX (K, RCKWRK, A6)
         DP    = REL*A6(1) + ABS
C
         DO 90 L = 1, NTR
            A6(L) = A6(L) + DP
   90    CONTINUE
C
         CALL CKRHEX (-K, RCKWRK, A6)
C
         CALL FUN (KK, JJ, NATJ, LBURNR, LENRGY, LMULTI, LVCOR, LTDIF,
     1             LVARMC, LTIME, JFIXT, TFIXT, FLRT, P, WT, EPS, DT,
     2             NTEMP, XGIVEN, TGIVEN, X, SN, S, WNDFAC, SCRCHK(1,1),
     3             YV, SCRCHK(1,2), SCRCHK(1,3), SCRCHK(1,4), COND, D,
     4             DKJ, TDR, ICKWRK, RCKWRK, IMCWRK, RMCWRK, F,
     5             SCRCHK(1,5), GFAC)
C
         DO 95 L = 1, NTR
            A6(L) = A6(L) - DP
   95    CONTINUE
         CALL CKRHEX (-K, RCKWRK, A6)
C
         DO 100 J = 1, JJ
            DO 100 N = 1, NATJ
               SN(N,J) =  - (F(N,J)-FN(N,J)) / DP
100      CONTINUE
C
C*****precision > double
         CALL DGBSL (A, LDA, NATJ*JJ, ML, MU, IP, SN, 0)
C*****END precision > double
C
C*****precision > single
C         CALL SGBSL (A, LDA, NATJ*JJ, ML, MU, IP, SN, 0)
C*****END precision > single
C
C          SN(N,J) NOW CONTAINS THE RAW SENSIVITY MATRIX DS(N,J)/DAi
C
C          NORMALIZE THE SENSIVITY COEFFICIENTS
C
         DO 200 J = 1, JJ
            SN(NT,J) = SN(NT,J) * A6(1) / S(NT,J)
200      CONTINUE
C
         IF (LMOLE) THEN
C
            DO 500 J = 1, JJ
               SUM = 0.0
               DO 300 L = 1, KK
                  SUM = SUM + SN(NYS+L,J) / WT(L)
300            CONTINUE
               CALL CKMMWY (S(NY,J), ICKWRK, RCKWRK, WTM)
               WTMSUM = WTM * SUM
               DO 400 L = 1, KK
                  SN(NYS+L,J) = A6(1) * (SN(NYS+L,J) / S(NYS+L,J)
     2                          - WTMSUM)
400            CONTINUE
500         CONTINUE
C
         ELSE
C
            DO 700 J = 1, JJ
               DO 600 L = 1, KK
                  SN(NYS+L,J) = SN(NYS+L,J) * A6(1) /
     1                          S(NYS+L,J)
600            CONTINUE
700         CONTINUE
C
         ENDIF
C
         DO 800 J = 1, JJ
            SN(NM,J) = SN(NM,J) * A6(1) / S(NM,J)
800      CONTINUE
C
         WRITE (LSAVE)  K, ((SN(N,J),N=1,NATJ), J=1,JJ)
         WRITE (LRCRVR) K, ((SN(N,J),N=1,NATJ), J=1,JJ)
C
1000  CONTINUE
C
      RETURN
      END
C------------------------------------------------------------------
C
      SUBROUTINE AREA (X,A)
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      A = 1.0
      RETURN
      END
C
      SUBROUTINE UPCASE(ISTR)
      CHARACTER ISTR*(*), LCASE(26)*1, UCASE(26)*1
      DATA LCASE /'a','b','c','d','e','f','g','h','i','j','k','l','m',
     1            'n','o','p','q','r','s','t','u','v','w','x','y','z'/,
     2     UCASE /'A','B','C','D','E','F','G','H','I','J','K','L','M',
     3            'N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/
C
      DO 10 J = 1, LEN(ISTR)
         DO 10 N = 1, 26
            IF (ISTR(J:J) .EQ. LCASE(N)) ISTR(J:J) = UCASE(N)
   10 CONTINUE
      RETURN
      END
