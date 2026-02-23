SUBROUTINE MXTURE(KLX,KLXS,KVX,KVXS,KIX,KIXS,KT1,KT2,PA,PB,PC,PD,PE,PY,PX)

!**** *MXTURE*   - Resolution of a set of triangular tridiagonal systems.

!     Purpose.    Resolution of a set of triangular tridiagonal systems.
!     --------

!    Example, for KLX=4:

!    PY is known, PX is unknown. We want to find PX, solution of the
!     following system, for l=1 to KVX, i=1 to KIX.

!    * KT = -2:

!    I PA(l,1)    0       0       0    I   I PX(l,1,i) I   I PY(l,1,i) I
!    I PB(l,1) PA(l,2)    0       0    I   I PX(l,2,i) I   I PY(l,2,i) I
!    I PC(l,1) PB(l,2) PA(l,3)    0    I * I PX(l,3,i) I = I PY(l,3,i) I
!    I    0    PC(l,2) PB(l,3) PA(l,4) I   I PX(l,4,i) I   I PY(l,4,i) I

!    * KT = -3:

!    I    1       0       0       0    I   I PX(l,1,i) I   I PY(l,1,i) I
!    I PB(l,1)    1       0       0    I   I PX(l,2,i) I   I PY(l,2,i) I
!    I PC(l,1) PB(l,2)    1       0    I * I PX(l,3,i) I = I PY(l,3,i) I
!    I    0    PC(l,2) PB(l,3)    1    I   I PX(l,4,i) I   I PY(l,4,i) I

!    Dummy array PA is not used in this case.

!    * KT = 1:

!    I    1    ZB(l,1) ZC(l,1)    0    I   I PX(l,1,i) I   I PY(l,1,i) I
!    I    0       1    ZB(l,2) ZC(l,2) I   I PX(l,2,i) I   I PY(l,2,i) I
!    I    0       0       1    ZB(l,3) I * I PX(l,3,i) I = I PY(l,3,i) I
!    I    0       0       0       1    I   I PX(l,4,i) I   I PY(l,4,i) I

!    where ZB(l,j)=PB(l,j)/PA(l,j); ZC(l,j)=PC(l,j)/PA(l,j) .

!    * KT = 2:

!    I PA(l,1) PD(l,1) PE(l,1)    0    I   I PX(l,1,i) I   I PY(l,1,i) I
!    I    0    PA(l,2) PD(l,2) PE(l,2) I   I PX(l,2,i) I   I PY(l,2,i) I
!    I    0       0    PA(l,3) PD(l,3) I * I PX(l,3,i) I = I PY(l,3,i) I
!    I    0       0       0    PA(l,4) I   I PX(l,4,i) I   I PY(l,4,i) I

!    * KT = 3:

!    I    1    PD(l,1) PE(l,1)    0    I   I PX(l,1,i) I   I PY(l,1,i) I
!    I    0       1    PD(l,2) PE(l,2) I   I PX(l,2,i) I   I PY(l,2,i) I
!    I    0       0       1    PD(l,3) I * I PX(l,3,i) I = I PY(l,3,i) I
!    I    0       0       0       1    I   I PX(l,4,i) I   I PY(l,4,i) I

!    Dummy array PA is not used in this case.

!**   Interface.
!     ----------
!        *CALL* *MXTURE(KLX,KLXS,KVX,KVXS,KIX,KIXS,KT1,KT2,&
!                        &PA,PB,PC,PD,PE,PY,PX)

!        Explicit arguments :
!        --------------------
!         KLX:        - Dimension of the system.                     (input)
!         KLXS:       - Surdimension of matrices PX,PY              (input)
!         KVX:        - Number of variables (second                  (input)
!                       dimension of PX and PY).
!         KVXS:       - Surdimension corresponding to KVX.           (input)
!         KIX:        - Number of variables multiplied by the same
!                       matrix.                                      (input)
!         KIXS:       - Surdimension of KIX in PX,PY                 (input)
!         KT1 :       - Type of matrix, 1st part (see figures above).(input)
!         KT2 :       - Type of matrix, 2nd part (see figures above).(input)
!         PA,PB,PC:   - non-zero diagonals of the triangular         (input)
!                       matrixes for 1st part and KT2=1 
!                       (see figures above).
!         PA,PD,PE:   - non-zero diagonals of the triangular         (input)
!                       matrixes for 2nd part (see figures above) 
!                       Caution! All PA coefficients must be non 0.
!         PY:         - known vector.                                (input)
!         PX:         - unknown vector.                              (output)

!        Implicit arguments :
!        --------------------
!        none.

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        None.
!        Called by SPC2,SPC2AD,SPCSIDG_PART1,SPCSIDG_PART1AD.

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        K. YESSAD: OCTOBER 1993.

!     Modifications.
!     --------------
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLX
INTEGER(KIND=JPIM),INTENT(IN)    :: KLXS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KVX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KVXS
INTEGER(KIND=JPIM),INTENT(IN)    :: KIX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIXS
INTEGER(KIND=JPIM),INTENT(IN)    :: KT1
INTEGER(KIND=JPIM),INTENT(IN)    :: KT2
REAL(KIND=JPRB)   ,INTENT(IN)    :: PA(KLX,KVX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PB(KLX,KVX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PC(KLX,KVX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PD(KLX,KVX)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PE(KLX,KVX)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PY(KIXS,KLX,KVXS) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PX(KIXS,KLX,KVXS) 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IT1,IT2, JI, JL, JV,IIXL,IIX

REAL(KIND=JPRB) :: ZBB, ZCC
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MXTURE',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!*       1.    PRELIMINARY INITIALISATIONS.
!              ----------------------------

IT1=MIN(3,MAX(-3,KT1))
IF (KT1 == 0.OR.KT1 == -1) IT1=-2
IT2=MIN(3,MAX(-3,KT2))
IF (KT2 == 0.OR.KT2 == -1) IT2=-2

IIXL=MAX(2,KIX)
IIX=KIX-MOD(KIX,2)

!      ----------------------------------------------------------------

!*       2.    COMPUTATION OF PX.
!              ------------------

! solves the linear system L.PY=PYin
!!$OMP PARALLEL PRIVATE(JI,JV,JL,ZBB,ZCC)
IF (IT1 == -3) THEN
  DO JI=1,IIX,2
    !!$OMP DO 
    DO JV=1,KVX
      PY(JI,1,JV)=PY(JI,1,JV)
      PY(JI+1,1,JV)=PY(JI+1,1,JV)
    ENDDO
  ENDDO
  DO JI=IIX+1,KIX
    DO JV=1,KVX
      PY(JI,1,JV)=PY(JI,1,JV)
    ENDDO
  ENDDO

  IF (KLX >= 2) THEN
    DO JI=1,IIX,2
      !!$OMP DO
      DO JV=1,KVX
        PY(JI,2,JV)=PY(JI,2,JV)-PB(1,JV)*PY(JI,1,JV)
        PY(JI+1,2,JV)=PY(JI+1,2,JV)-PB(1,JV)*PY(JI+1,1,JV)
      ENDDO
    ENDDO
    DO JI=IIX+1,KIX
      DO JV=1,KVX
        PY(JI,2,JV)=PY(JI,2,JV)-PB(1,JV)*PY(JI,1,JV)
      ENDDO
    ENDDO
  ENDIF

  IF (KLX >= 3) THEN
    DO JI=1,IIX,2
      !!$OMP DO SIMD
      DO JV=1,KVX
        DO JL=3,KLX
          PY(JI,JL,JV)=PY(JI,JL,JV)-PB(JL-1,JV)*PY(JI,JL-1,JV)&
             &-PC(JL-2,JV)*PY(JI,JL-2,JV)
          PY(JI+1,JL,JV)=PY(JI+1,JL,JV)-PB(JL-1,JV)*PY(JI+1,JL-1,JV)&
             &-PC(JL-2,JV)*PY(JI+1,JL-2,JV)
        ENDDO
      ENDDO
    ENDDO
    DO JI=IIX+1,KIX
      !!$OMP DO SIMD
      DO JV=1,KVX
        DO JL=3,KLX
           PY(JI,JL,JV)=PY(JI,JL,JV)-PB(JL-1,JV)*PY(JI,JL-1,JV)&
             &-PC(JL-2,JV)*PY(JI,JL-2,JV)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

ELSEIF (IT1==-2) THEN
  DO JI=1,IIX,2
    !!$OMP DO
    DO JV=1,KVX
      PY(JI,1,JV)=PY(JI,1,JV)/PA(1,JV)
      PY(JI+1,1,JV)=PY(JI+1,1,JV)/PA(1,JV)    
    ENDDO
  ENDDO
  DO JI=IIX+1,KIX
    !!$OMP DO
    DO JV=1,KVX
      PY(JI,1,JV)=PY(JI,1,JV)/PA(1,JV)
    ENDDO
  ENDDO

  IF (KLX >= 2) THEN
    DO JI=1,IIX,2
      !!$OMP DO 
      DO JV=1,KVX
        PY(JI,2,JV)=(PY(JI,2,JV)-PB(1,JV)*PY(JI,1,JV))/PA(2,JV)
        PY(JI+1,2,JV)=(PY(JI+1,2,JV)-PB(1,JV)*PY(JI+1,1,JV))/PA(2,JV)
      ENDDO
    ENDDO
    DO JI=IIX+1,KIX
      !!$OMP DO
      DO JV=1,KVX
        PY(JI,2,JV)=(PY(JI,2,JV)-PB(1,JV)*PY(JI,1,JV))/PA(2,JV)
      ENDDO
    ENDDO
  ENDIF

  IF (KLX >= 3) THEN
    DO JI=1,IIX,2
      !!$OMP DO SIMD
      DO JV=1,KVX
        DO JL=3,KLX
          PY(JI,JL,JV)=(PY(JI,JL,JV)-PC(JL-2,JV)*PY(JI,JL-2,JV)&
           &-PB(JL-1,JV)*PY(JI,JL-1,JV))/PA(JL,JV)
          PY(JI+1,JL,JV)=(PY(JI+1,JL,JV)-PC(JL-2,JV)*PY(JI+1,JL-2,JV)&
           &-PB(JL-1,JV)*PY(JI+1,JL-1,JV))/PA(JL,JV)
        ENDDO
      ENDDO
    ENDDO
    DO JI=IIX+1,KIX
      !!$OMP DO SIMD
      DO JV=1,KVX
        DO JL=3,KLX
          PY(JI,JL,JV)=(PY(JI,JL,JV)-PC(JL-2,JV)*PY(JI,JL-2,JV)&
           &-PB(JL-1,JV)*PY(JI,JL-1,JV))/PA(JL,JV)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

ENDIF

! solves U.PX=Linv.PYin with Linv.PYin in PY
IF (IT2 == 1) THEN
  DO JI=1,IIX,2
    !!$OMP DO
    DO JV=1,KVX
      PX(JI,KLX,JV)=PY(JI,KLX,JV)
      PX(JI+1,KLX,JV)=PY(JI+1,KLX,JV)
    ENDDO
  ENDDO
  DO JI=IIX+1,KIX
    !!$OMP DO
    DO JV=1,KVX
      PX(JI,KLX,JV)=PY(JI,KLX,JV)
    ENDDO
  ENDDO

  IF (KLX >= 2) THEN
    DO JI=1,IIX,2
      !!$OMP DO
      DO JV=1,KVX
        ZBB=PB(KLX-1,JV)/PA(KLX-1,JV)
        PX(JI,KLX-1,JV)=PY(JI,KLX-1,JV)-ZBB*PX(JI,KLX,JV)
        PX(JI+1,KLX-1,JV)=PY(JI+1,KLX-1,JV)-ZBB*PX(JI+1,KLX,JV)
      ENDDO
    ENDDO
    DO JI=IIX+1,KIX
      !!$OMP DO
      DO JV=1,KVX
        ZBB=PB(KLX-1,JV)/PA(KLX-1,JV)
        PX(JI,KLX-1,JV)=PY(JI,KLX-1,JV)-ZBB*PX(JI,KLX,JV)
      ENDDO
    ENDDO
  ENDIF

  IF (KLX >= 3) THEN
    DO JI=1,IIX,2
      !!$OMP DO SIMD
      DO JV=1,KVX
        DO JL=KLX-2,1,-1
          ZBB=PB(JL,JV)/PA(JL,JV)
          ZCC=PC(JL,JV)/PA(JL,JV)
          PX(JI,JL,JV)=PY(JI,JL,JV)-ZBB*PX(JI,JL+1,JV)-ZCC*PX(JI,JL+2,JV)
          PX(JI+1,JL,JV)=PY(JI+1,JL,JV)-ZBB*PX(JI+1,JL+1,JV)-ZCC*PX(JI+1,JL+2,JV)
        ENDDO
      ENDDO
    ENDDO
    DO JI=1+IIX,KIX
      !!$OMP DO SIMD
      DO JV=1,KVX
        DO JL=KLX-2,1,-1
          ZBB=PB(JL,JV)/PA(JL,JV)
          ZCC=PC(JL,JV)/PA(JL,JV)
          PX(JI,JL,JV)=PY(JI,JL,JV)-ZBB*PX(JI,JL+1,JV)-ZCC*PX(JI,JL+2,JV)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

ELSEIF (IT2==2) THEN
  DO JI=1,IIX,2
    !!$OMP DO
    DO JV=1,KVX
      PX(JI,KLX,JV)=PY(JI,KLX,JV)/PA(KLX,JV)
      PX(JI+1,KLX,JV)=PY(JI+1,KLX,JV)/PA(KLX,JV)
    ENDDO
  ENDDO
  DO JI=IIX+1,KIX
    !!$OMP DO
    DO JV=1,KVX
      PX(JI,KLX,JV)=PY(JI,KLX,JV)/PA(KLX,JV)
    ENDDO
  ENDDO

  IF (KLX >= 2) THEN
    DO JI=1,IIX,2
      !!$OMP DO
      DO JV=1,KVX
        PX(JI,KLX-1,JV)=&
         & (PY(JI,KLX-1,JV)-PD(KLX-1,JV)*PX(JI,KLX,JV))/PA(KLX-1,JV)
        PX(JI+1,KLX-1,JV)=&
         & (PY(JI+1,KLX-1,JV)-PD(KLX-1,JV)*PX(JI+1,KLX,JV))/PA(KLX-1,JV)
      ENDDO
    ENDDO
    DO JI=IIX+1,KIX
      !!$OMP DO
      DO JV=1,KVX
        PX(JI,KLX-1,JV)=&
         & (PY(JI,KLX-1,JV)-PD(KLX-1,JV)*PX(JI,KLX,JV))/PA(KLX-1,JV)
      ENDDO
    ENDDO
  ENDIF

  IF (KLX >= 3) THEN
    DO JI=1,IIX,2
      !!$OMP DO SIMD
      DO JV=1,KVX
        DO JL=KLX-2,1,-1
          PX(JI,JL,JV)=(PY(JI,JL,JV)-PD(JL,JV)*PX(JI,JL+1,JV)&
           &-PE(JL,JV)*PX(JI,JL+2,JV))/PA(JL,JV)
          PX(JI+1,JL,JV)=(PY(JI+1,JL,JV)-PD(JL,JV)*PX(JI+1,JL+1,JV)&
           &-PE(JL,JV)*PX(JI+1,JL+2,JV))/PA(JL,JV)
        ENDDO
      ENDDO
    ENDDO
    DO JI=IIX+1,KIX
      !!$OMP DO SIMD
      DO JV=1,KVX
        DO JL=KLX-2,1,-1
          PX(JI,JL,JV)=(PY(JI,JL,JV)-PD(JL,JV)*PX(JI,JL+1,JV)&
           &-PE(JL,JV)*PX(JI,JL+2,JV))/PA(JL,JV)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

ELSEIF (IT2 == 3) THEN
  DO JI=1,IIX,2
    !!$OMP DO
    DO JV=1,KVX
      PX(JI,KLX,JV)=PY(JI,KLX,JV)
      PX(JI+1,KLX,JV)=PY(JI+1,KLX,JV)
    ENDDO
  ENDDO
  DO JI=IIX+1,KIX
    !!$OMP DO
    DO JV=1,KVX
      PX(JI,KLX,JV)=PY(JI,KLX,JV)
    ENDDO
  ENDDO

  IF (KLX >= 2) THEN
    DO JI=1,IIX,2
      !!$OMP DO
      DO JV=1,KVX
        ZBB=PD(KLX-1,JV)
        PX(JI,KLX-1,JV)=PY(JI,KLX-1,JV)-ZBB*PX(JI,KLX,JV)
        PX(JI+1,KLX-1,JV)=PY(JI+1,KLX-1,JV)-ZBB*PX(JI+1,KLX,JV)
      ENDDO
    ENDDO
    DO JI=IIX+1,KIX
      !!$OMP DO
      DO JV=1,KVX
        ZBB=PD(KLX-1,JV)
        PX(JI,KLX-1,JV)=PY(JI,KLX-1,JV)-ZBB*PX(JI,KLX,JV)
      ENDDO
    ENDDO
  ENDIF

  IF (KLX >= 3) THEN
    DO JI=1,IIX,2
      !!$OMP DO SIMD
      DO JV=1,KVX
        DO JL=KLX-2,1,-1
          ZBB=PD(JL,JV)
          ZCC=PE(JL,JV)
          PX(JI,JL,JV)=PY(JI,JL,JV)-ZBB*PX(JI,JL+1,JV)-ZCC*PX(JI,JL+2,JV)
          PX(JI+1,JL,JV)=PY(JI,JL,JV)-ZBB*PX(JI+1,JL+1,JV)-ZCC*PX(JI+1,JL+2,JV)
        ENDDO
      ENDDO
    ENDDO
    DO JI=IIX+1,KIX
      !!$OMP DO SIMD
      DO JV=1,KVX
        DO JL=KLX-2,1,-1
          ZBB=PD(JL,JV)
          ZCC=PE(JL,JV)
          PX(JI,JL,JV)=PY(JI,JL,JV)-ZBB*PX(JI,JL+1,JV)-ZCC*PX(JI,JL+2,JV)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

ENDIF
!!$OMP END PARALLEL
!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('MXTURE',1,ZHOOK_HANDLE)
END SUBROUTINE MXTURE

