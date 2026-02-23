#if defined(USE_OPENACC)
SUBROUTINE SGEMMXACC(PA,KA,KAD,PB,KB,KBD,PC,KC,KCA,KAR,KAC,KBC)

!**   Interface.
!     ----------
!        CALL SGEMMXACC(PA,KA,KAD,PB,KB,KBD,PC,KC,KCA,KAR,KAC,KBC)

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

USE OPENACC
USE CUBLAS

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

!     ------------------------------------------------------------------

REAL(KIND=JPRB)   :: PCTR(KAR,KBC)
REAL(KIND=JPRB)   :: PBTR(KAC,KBC)
INTEGER(KIND=JPIM):: JL,JS
LOGICAL           :: LLACC

REAL(KIND=JPRB), PARAMETER :: ALPHA=1.0
REAL(KIND=JPRB), PARAMETER :: BETA=0.0

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SGEMMXACC',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!$ACC DATA PRESENT(PA,PB,PC)
!$ACC HOST_DATA USE_DEVICE(PA,PB,PC)
#if defined(PARKIND1_SINGLE)

  CALL CUBLASSGEMM('N','T',KBC,KAR,KAC,ALPHA,&
        &PB(1,1),KBC,PA(1,1),KAR,BETA,PC(1,1),KBC)
#else

  CALL CUBLASDGEMM('N','T',KBC,KAR,KAC,ALPHA,&
        &PB(1,1),KBC,PA(1,1),KAR,BETA,PC(1,1),KBC)

#endif
!$ACC END HOST_DATA
!$ACC WAIT
!$ACC END DATA

IF (LHOOK) CALL DR_HOOK('SGEMMXACC',1,ZHOOK_HANDLE)
END SUBROUTINE SGEMMXACC
#elif defined(USE_OPENMP)
SUBROUTINE SGEMMXACC(PA,KA,KAD,PB,KB,KBD,PC,KC,KCA,KAR,KAC,KBC)

!**   Interface.
!     ----------
!        CALL SGEMMXACC(PA,KA,KAD,PB,KB,KBD,PC,KC,KCA,KAR,KAC,KBC)

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
!USE OML_MOD   ,ONLY : OML_MAX_THREADS, OML_IN_PARALLEL

USE HIPFORT
USE HIPFORT_ROCBLAS
USE ISO_C_BINDING

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

!     ------------------------------------------------------------------

REAL(KIND=JPRB)   :: PCTR(KAR,KBC)
REAL(KIND=JPRB)   :: PBTR(KAC,KBC)
INTEGER(KIND=JPIM):: JL,JS,IRESU
LOGICAL           :: LLACC

REAL(KIND=JPRB), PARAMETER :: ALPHA=1.0
REAL(KIND=JPRB), PARAMETER :: BETA=0.0

TYPE(C_PTR), SAVE :: YLCBH
TYPE(C_PTR) :: handle
integer ::status
LOGICAL, SAVE :: LLCBH_INITIALIZED = .FALSE.
TYPE(C_PTR), SAVE :: STREAM 

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SGEMMXACC',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------
IF (.NOT. LLCBH_INITIALIZED) THEN
  CALL CHECKROCBLAS(ROCBLAS_CREATE_HANDLE (YLCBH))
  CALL CHECKROCBLAS(ROCBLAS_GET_STREAM(YLCBH,STREAM))
  LLCBH_INITIALIZED = .TRUE.
ENDIF

!$OMP TARGET DATA MAP(TOFROM:PA,PB,PC)
!$OMP TARGET DATA USE_DEVICE_ADDR(PA,PB,PC)
#if defined(PARKIND1_SINGLE)

  CALL CHECKROCBLAS( ROCBLAS_SGEMM(YLCBH,ROCBLAS_OPERATION_NONE,ROCBLAS_OPERATION_TRANSPOSE,KBC,KAR,KAC,ALPHA,&
        &PB(1,1),KBC,PA(1,1),KAR,BETA,PC(1,1),KBC))
#else

  CALL CHECKROCBLAS(ROCBLAS_DGEMM(YLCBH,ROCBLAS_OPERATION_NONE,ROCBLAS_OPERATION_TRANSPOSE,KBC,KAR,KAC,ALPHA,&
        &PB(1,1),KBC,PA(1,1),KAR,BETA,PC(1,1),KBC))

#endif
!$OMP END TARGET DATA
!$OMP END TARGET DATA
CALL CHECKHIP (HIPSTREAMSYNCHRONIZE (STREAM))

IF (LHOOK) CALL DR_HOOK('SGEMMXACC',1,ZHOOK_HANDLE)
CONTAINS

SUBROUTINE CHECKROCBLAS (STATUS)

USE HIPFORT
USE HIPFORT_ROCBLAS

INTEGER :: STATUS

IF (STATUS /= ROCBLAS_STATUS_SUCCESS) THEN
  PRINT *, 'ROCBLAS ERROR: STATUS =', STATUS
  PRINT *, 'ROCBLAS_STATUS_INVALID_HANDLE=',ROCBLAS_STATUS_INVALID_HANDLE
  CALL ABOR1 ('FXTRAN_ACDC_GEMM')
END IF

END SUBROUTINE CHECKROCBLAS

SUBROUTINE CHECKHIP (STATUS)

USE HIPFORT
USE HIPFORT_ROCBLAS

INTEGER :: STATUS

IF (STATUS /= HIPSUCCESS) THEN
  PRINT *, 'HIP ERROR: STATUS = ', STATUS
  CALL ABOR1 ('FXTRAN_ACDC_GEMM')
END IF

END SUBROUTINE CHECKHIP

END SUBROUTINE SGEMMXACC

#else
SUBROUTINE SGEMMXACC(PA,KA,KAD,PB,KB,KBD,PC,KC,KCA,KAR,KAC,KBC)

!**   Interface.
!     ----------
!        CALL SGEMMXACC(PA,KA,KAD,PB,KB,KBD,PC,KC,KCA,KAR,KAC,KBC)

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

!     ------------------------------------------------------------------

REAL(KIND=JPRB)   :: PCTR(KAR,KBC)
REAL(KIND=JPRB)   :: PBTR(KAC,KBC)
INTEGER(KIND=JPIM):: JL,JS
LOGICAL           :: LLACC

REAL(KIND=JPRB), PARAMETER :: ALPHA=1.0
REAL(KIND=JPRB), PARAMETER :: BETA=0.0

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SGEMMXACC',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

#if defined(PARKIND1_SINGLE)

  CALL SGEMM('N','T',KBC,KAR,KAC,ALPHA,&
        &PB(1,1),KBC,PA(1,1),KAR,BETA,PC(1,1),KBC)
#else

  CALL DGEMM('N','T',KBC,KAR,KAC,ALPHA,&
        &PB(1,1),KBC,PA(1,1),KAR,BETA,PC(1,1),KBC)
#endif

IF (LHOOK) CALL DR_HOOK('SGEMMXACC',1,ZHOOK_HANDLE)
END SUBROUTINE SGEMMXACC

#endif
