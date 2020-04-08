      PROGRAM DRIVER
C*****precision > double
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
C          LENLWK allocates the logical working space
C          LENIWK allocates the integer working space
C          LENRWK allocates the real working space
C          LENCWK allocates the character working space
C          LENSYM is the length of a character string
C          LIN is the unit from which user input is read
C          LOUT is the unit to which printed output is written
C          LRIN is the unit from which the restart file is read
C          LROUT is the unit to which the save file is written
C          LRCRVR is the unit to which the scratch save file is written
C          LINKCK is unit from which the Chemkin Linking file is read
C          LINTP is unit from which the Transport Linking file is read
C
      PARAMETER (LENLWK=570, LENIWK=42000, LENRWK=17000000, LENCWK=140,
     1           LENSYM=16, LIN=5, LOUT=6, LRIN=14, LROUT=15,
     2           LRCRVR=16, LINKCK=25, LINKTP=35)
C
C          All storage needed by the flame code is allocated in the
C          following three arrays.  LWORK is for logical variables,
C          IWORK is for integer variables, RWORK is of real variables,
C          CWORK is for character variables.
C
      DIMENSION IWORK(LENIWK), RWORK(LENRWK)
      CHARACTER CWORK(LENCWK)*(LENSYM)
      LOGICAL LWORK(LENLWK)
C
      NMAX = 500
C            open the user input file
C
C*****VAX/vms
C      OPEN (LIN, STATUS='OLD', FORM='FORMATTED')
C      OPEN (LOUT, STATUS='UNKNOWN', FORM='FORMATTED')
C      OPEN (LRIN, STATUS='OLD', FORM='UNFORMATTED')
C      OPEN (LROUT, STATUS='UNKNOWN', FORM='UNFORMATTED')
C      OPEN (LRCRVR, STATUS='UNKNOWN', FORM='UNFORMATTED')
C      OPEN (LINKCK, STATUS='OLD', FORM='UNFORMATTED')
C      OPEN (LINKTP, STATUS='OLD', FORM='UNFORMATTED')
C*****END VAX/vms
C
C*****unix
      OPEN (LIN, FORM='FORMATTED', FILE='../input/inp')
      OPEN (LOUT, FORM='FORMATTED', FILE='../output/flout')
      OPEN (LRIN, FORM='UNFORMATTED', FILE='../output/restart')
      OPEN (LROUT, FORM='UNFORMATTED', FILE='../output/save')
      OPEN (LRCRVR, FORM='UNFORMATTED', FILE='../output/recover')
      OPEN (LINKCK, FORM='UNFORMATTED', FILE='../link/cklink')
      OPEN (LINKTP, FORM='UNFORMATTED', FILE='../link/tplink')
C*****END unix
C
      CALL PREMIX (NMAX, LIN, LOUT, LINKCK, LINKTP, LRIN, LROUT,
     1             LRCRVR, LENLWK, LWORK, LENIWK, IWORK, LENRWK,
     2             RWORK, LENCWK, CWORK)
C
      STOP
      END
