SUBROUTINE MXMAOPTR(PA,KA,KAD,PB,KB,KBD,PC,KC,KCA,KAR,KAC,KBC,LDACC)

!**** *MXMAOPTR - Prepares call to MXMAOP

!     Purpose.
!     --------
!        Transposition of PB and PC befor call to MXMAOP for bit reproducibility.
!

!**   Interface.
!     ----------
!        CALL MXMAOPTR(PA,KA,KAD,PB,KB,KBD,PC,KC,KCA,KAR,KAC,KBC)

!        Explicit arguments : See SGEMMX documentaion.
!        --------------------
!         PA     - input matrix PA                                     (in)
!         KA     - memory jump between two lines in PA (generally 1)   (in)
!         KAD    - memory jump between two columns in PA               (in)
!         PB     - input matrix PB                                     (in)
!         KB     - memory jump between two lines in PB (generally 1)   (in)
!         KBD    - memory jump between two columns in PB               (in)
!         PC     - output matrix PC                                    (out)
!         KC     - memory jump between two lines in PC (generally 1)   (in)
!         KCA    - memory jump between two columns in PC               (in)
!         KAR    - number of useful lines of PA                        (in)
!         KAC    - number of useful columns of PA (and lines of PB)    (in)
!         KBC    - number of useful columns of PB                      (in)

!        Implicit arguments :
!        --------------------

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE OML_MOD   ,ONLY : OML_MAX_THREADS, OML_IN_PARALLEL

!     ------------------------------------------------------------------

IMPLICIT NONE

REAL(KIND=JPRB)   ,INTENT(IN)    :: PA(KAR,KAC) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KAD 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PB(KBC,KAC) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KB 
INTEGER(KIND=JPIM),INTENT(IN)    :: KBD 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PC(KBC,KAR) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KC 
INTEGER(KIND=JPIM),INTENT(IN)    :: KCA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KAR 
INTEGER(KIND=JPIM),INTENT(IN)    :: KAC 
INTEGER(KIND=JPIM),INTENT(IN)    :: KBC 
LOGICAL, OPTIONAL, INTENT(IN)    :: LDACC

!     ------------------------------------------------------------------

REAL(KIND=JPRB)   :: PCTR(KAR,KBC)
REAL(KIND=JPRB)   :: PBTR(KAC,KBC)
INTEGER(KIND=JPIM):: JL,JS
LOGICAL           :: LLACC


REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
#include "mxmaop.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MXMAOPTR',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

LLACC=.FALSE.
IF (PRESENT(LDACC)) LLACC=LDACC

IF (LLACC) THEN
!  write (0,*) "branche mxmaoptr llacc" 
  CALL MXMAOP(PA,KA,KAD,PB,KB,KBD,PC,KC,KCA,KAR,KAC,KBC,LLACC)

ELSE

  !$OMP PARALLEL DO PRIVATE(JL,JS)
  DO JS=1,KBC
    DO JL=1,KAC
      PBTR(JL,JS)=PB(JS,JL)
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  
  CALL MXMAOP(PA,KA,KAD,PBTR,KB,KBD,PCTR,KC,KCA,KAR,KAC,KBC)
  
  !$OMP PARALLEL DO PRIVATE(JL,JS)
  DO JS=1,KBC
    DO JL=1,KAR
      PC(JS,JL)=PCTR(JL,JS)
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO

ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('MXMAOPTR',1,ZHOOK_HANDLE)
END SUBROUTINE MXMAOPTR

