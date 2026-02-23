SUBROUTINE SPCSIDG_PART1ACC (YDGEOMETRY, YDDYN, KSTA,KEND,PSDIVP,&
          &PSPDIVP,KMLOCSTA,KMLOCEND,LDSIDG,KSZNISNAX,NISNAX,SIHEG,SIHEG2)
!SPCSIDG_PART1ACC : SPECTRAL SPACE SEMI-IMPLICIT CALCULATION, FIRST HORIZONTAL
!OPERATOR. INVERSION OF A PENTADIAGONAL SYSTEM IN STRETCHED GRIDS.
!
!     YDGEOMETRY, YDDYN (IN) : SET-UP VARIABLES
!     KSTA                   : First line (spectral number) processed
!     KEND                   : Last line (spectral number) processed
!     PSDIVP (INOUT)         : Divergence, before inversion
!     PSPDIVP (INOUT)        : Divergence, after inversion
!     KMLOCSTA, KMLOCEND     : Range of JMLOC computed on this task
!     LSIDG, KSZNISNAX,NISNAX: Not used yet (stretched-grid for local model). 
!     SIHEG, SIHEG2 (IN)     : Parameter arrays used in the inversion, LU
!     factorized in setup
!
!
!
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMDYN       , ONLY : TDYN

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTA
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND
INTEGER(KIND=JPIM),INTENT(IN)    :: KMLOCSTA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KMLOCEND
LOGICAl           ,INTENT(IN)    :: LDSIDG
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSDIVP(KSTA:KEND,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPDIVP(KSTA:KEND,YDGEOMETRY%YRDIMV%NFLEVG)
INTEGER(KIND=JPIM),INTENT(IN)    :: KSZNISNAX
INTEGER(KIND=JPIM),INTENT(IN)    :: NISNAX(0:KSZNISNAX)
REAL(KIND=JPRB)   ,INTENT(IN)    :: SIHEG(YDGEOMETRY%YRDIMV%NFLEVG*&
                                          &YDGEOMETRY%YRMP%NSPEC2VF/2,3)
REAL(KIND=JPRB)   ,INTENT(IN)    :: SIHEG2(YDGEOMETRY%YRDIM%NSMAX+1,&
                                          &YDGEOMETRY%YRDIMV%NFLEVG,2:3)
!logical , intent(in)  :: llacc
!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JV,IIS,II0

INTEGER(KIND=JPIM) :: SIHEGSTA,SIHEGEND

INTEGER(KIND=JPIM), PARAMETER    :: ISIZESP=54!!62
INTEGER(KIND=JPIM), PARAMETER    :: ISIZELEV=1 !0!!8!!54-10 meilleur pour l'instant
REAL(KIND=JPRB) :: ZPAS(ISIZESP+3,ISIZELEV)
REAL(KIND=JPRB) :: ZPBS(ISIZESP+3,ISIZELEV)
REAL(KIND=JPRB) :: ZPCS(ISIZESP+3,ISIZELEV)
REAL(KIND=JPRB) :: ZPIN(ISIZESP+3,ISIZELEV,2)
INTEGER(KIND=JPIM) :: II, IS0, IS02, ISE, JN,JLEVTOT,JMLOC,JI,IOFFSET1,KLX
INTEGER(KIND=JPIM) :: IM, ISTA, IEND,JSP,JSPTOT,JLB,ISPBLOC,IRMN,JLEV,JLEVBLOC

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!     ------------------------------------------------------------------

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SPCSIDG_PART1ACC',0,ZHOOK_HANDLE)

!             Inversion of two tridiagonal systems (Helmholtz equation)
!                --> (SIMI*DIVprIM(t+dt)).
#ifdef USE_OPENACC
!$ACC DATA PRESENT(SIHEG,SIHEG2,YDDYN%SIHEGIND,&
!$ACC& PSDIVP,PSPDIVP,YDGEOMETRY%YRDIM%NSMAX,YDGEOMETRY%YRDIMV%NFLEVG,&
!$ACC& YDGEOMETRY%YRLAP%MYMS,YDGEOMETRY%YRMP%NSPSTAF) CREATE(ZPAS,ZPBS,ZPCS,ZPIN)
!$ACC PARALLEL PRIVATE(IM,ISTA,KLX,ISE,ZPAS,ZPBS,ZPCS,ZPIN,JLB,ISPBLOC,IRMN,SIHEGSTA,SIHEGEND) DEFAULT(NONE)
!$ACC CACHE(ZPAS(1:ISIZESP+3,ISIZELEV),ZPBS(1:ISIZESP+3,ISIZELEV),ZPCS(1:ISIZESP+3,ISIZELEV),ZPIN(1:ISIZESP+3,1:ISIZELEV,1:2))
!$ACC LOOP GANG COLLAPSE(2) 
#endif
#ifdef USE_OPENMP
!$OMP TARGET DATA MAP(PRESENT:SIHEG,SIHEG2,YDDYN%SIHEGIND,&
!$OMP& PSDIVP,PSPDIVP,YDGEOMETRY%YRDIM%NSMAX,YDGEOMETRY%YRDIMV%NFLEVG,&
!$OMP& YDGEOMETRY%YRLAP%MYMS,YDGEOMETRY%YRMP%NSPSTAF) MAP(ALLOC:ZPAS,ZPBS,ZPCS,ZPIN) !! if(llacc)
!$OMP TARGET TEAMS DISTRIBUTE PRIVATE(IM,ISTA,KLX,ISE,ZPAS,ZPBS,ZPCS,ZPIN,JLB,ISPBLOC,IRMN,SIHEGSTA,SIHEGEND) COLLAPSE(2) THREAD_LIMIT(64) !! if(llacc)
#endif
DO JMLOC=KMLOCSTA,KMLOCEND                               !!each block inverts one JMLOC
  DO JLEVBLOC=1,(YDGEOMETRY%YRDIMV%NFLEVG-1)/ISIZELEV+1  !!and ISIZELEV vertical levels
    IM=YDGEOMETRY%YRLAP%MYMS(JMLOC)
    ISTA=YDGEOMETRY%YRMP%NSPSTAF(IM)
    SIHEGSTA=YDDYN%SIHEGIND(JMLOC)
    SIHEGEND=YDDYN%SIHEGIND(JMLOC+1)-1
    KLX=YDGEOMETRY%YRDIM%NSMAX+1-IM
    
    IF (IM > 0) THEN

      !               Inversion of a symmetric matrix.

      !!solves L.(U.X)=Y for U.X
      IF (KLX >= 3) THEN
#ifdef USE_OPENACC        
        !$ACC LOOP SEQ
#endif
        DO JLB=1,(KLX-3)/ISIZESP+1  !!processes KLX-3+1 elements, from 3 to KLX
          ISPBLOC=(JLB-1)*ISIZESP
#ifdef USE_OPENACC
          !$ACC LOOP VECTOR PRIVATE(JLEV,JLEVTOT,IOFFSET1,ISE)
#endif

          DO JSPTOT=ISPBLOC+1,MIN(ISPBLOC+ISIZESP+2,KLX)      !!pre-loads parts of the data in shared memory
            ISE=ISTA+2*(JSPTOT-1)
            DO JLEV=1,MIN(ISIZELEV,YDGEOMETRY%YRDIMV%NFLEVG-(JLEVBLOC-1)*ISIZELEV)
              JLEVTOT=JLEV+(JLEVBLOC-1)*ISIZELEV
              IOFFSET1=SIHEGSTA-1+(JLEVTOT-1)*(YDGEOMETRY%YRDIM%NSMAX+1-IM)
              ZPAS(JSPTOT-ISPBLOC,JLEV)=SIHEG(IOFFSET1+JSPTOT,1)
              ZPBS(JSPTOT-ISPBLOC,JLEV)=SIHEG(IOFFSET1+JSPTOT,2)
              ZPCS(JSPTOT-ISPBLOC,JLEV)=SIHEG(IOFFSET1+JSPTOT,3)
              DO JI=1,2
                ZPIN(JSPTOT-ISPBLOC,JLEV,JI)=PSDIVP(ISE+JI-1,JLEVTOT)
              ENDDO
            ENDDO
          ENDDO
          IF (JLB==1) THEN                                     !!processes the pre-loaded data, if first block
#ifdef USE_OPENACC

            !$ACC LOOP VECTOR PRIVATE(JLEVTOT,JSP) COLLAPSE(2)
#endif

            DO JLEV=1,MIN(ISIZELEV,YDGEOMETRY%YRDIMV%NFLEVG-(JLEVBLOC-1)*ISIZELEV)
              DO JI=1,2
                ZPIN(1,JLEV,JI)=ZPIN(1,JLEV,JI)/ZPAS(1,JLEV)
                ZPIN(2,JLEV,JI)=(ZPIN(2,JLEV,JI)-ZPBS(1,JLEV)*ZPIN(1,JLEV,JI))/ZPAS(2,JLEV)
                DO JSP=3,MIN(ISPBLOC+ISIZESP+2,KLX)-ISPBLOC
                  ZPIN(JSP,JLEV,JI)=(ZPIN(JSP,JLEV,JI)-ZPCS(JSP-2,JLEV)*ZPIN(JSP-2,JLEV,JI)&
                   & -ZPBS(JSP-1,JLEV)*ZPIN(JSP-1,JLEV,JI))/ZPAS(JSP,JLEV)  
                ENDDO
              ENDDO
            ENDDO
          ELSE                                                !!processes the pre-loaded data, if following blocks
            ISE=ISTA+2*(ISPBLOC+1-1)
#ifdef USE_OPENACC

            !$ACC LOOP VECTOR PRIVATE(JLEVTOT,JSP) COLLAPSE(2)
#endif

            DO JLEV=1,MIN(ISIZELEV,YDGEOMETRY%YRDIMV%NFLEVG-(JLEVBLOC-1)*ISIZELEV)
              DO JI=1,2
                JLEVTOT=JLEV+(JLEVBLOC-1)*ISIZELEV
                ZPIN(1,JLEV,JI)=PSPDIVP(ISE+JI-1,JLEVTOT)
                ZPIN(2,JLEV,JI)=PSPDIVP(ISE+JI+1,JLEVTOT)
                DO JSP=3,MIN(ISPBLOC+ISIZESP+2,KLX)-ISPBLOC
                  ZPIN(JSP,JLEV,JI)=(ZPIN(JSP,JLEV,JI)-ZPCS(JSP-2,JLEV)*ZPIN(JSP-2,JLEV,JI)&
                   & -ZPBS(JSP-1,JLEV)*ZPIN(JSP-1,JLEV,JI))/ZPAS(JSP,JLEV)  
                ENDDO
              ENDDO
            ENDDO
          ENDIF
#ifdef USE_OPENACC

          !$ACC LOOP VECTOR PRIVATE(JLEV,JI,ISE)
#endif

          DO JSPTOT=ISPBLOC+1,MIN(ISPBLOC+ISIZESP+2,KLX)     !!moves the processed data from shared memory to main variable
            ISE=ISTA+2*(JSPTOT-1)
            DO JLEV=1,MIN(ISIZELEV,YDGEOMETRY%YRDIMV%NFLEVG-(JLEVBLOC-1)*ISIZELEV)
              DO JI=1,2
                PSPDIVP(ISE+JI-1,JLEV+(JLEVBLOC-1)*ISIZELEV)=ZPIN(JSPTOT-ISPBLOC,JLEV,JI)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ELSE                                                   !!direct processing, without shared memory, if small-size data
        ISE=ISTA+2*(1-1)
#ifdef USE_OPENACC

        !$ACC LOOP VECTOR PRIVATE(JLEVTOT,IOFFSET1) COLLAPSE(2)
#endif

        DO JLEV=1,MIN(ISIZELEV,YDGEOMETRY%YRDIMV%NFLEVG-(JLEVBLOC-1)*ISIZELEV)
          DO JI=1,2
            JLEVTOT=JLEV+(JLEVBLOC-1)*ISIZELEV
            IOFFSET1=SIHEGSTA-1+(JLEVTOT-1)*(YDGEOMETRY%YRDIM%NSMAX+1-IM)
            PSPDIVP(ISE+JI-1,JLEVTOT)=PSDIVP(ISE+JI-1,JLEVTOT)/SIHEG(IOFFSET1+1,1)
            IF (KLX >= 2) THEN
              PSPDIVP(ISE+JI+1,JLEVTOT)=(PSDIVP(ISE+JI+1,JLEVTOT)&
               &-SIHEG(IOFFSET1+1,2)*PSPDIVP(ISE+JI-1,JLEVTOT))/SIHEG(IOFFSET1+2,1)
            ENDIF
          ENDDO
        ENDDO
      ENDIF

      !!solves U.X=Linv.Y for X
      IF (KLX >= 3) THEN
        IRMN=MOD(KLX-3,ISIZESP)+1
#ifdef USE_OPENACC
        !$ACC LOOP SEQ
#endif
        DO JLB=(KLX-3)/ISIZESP+1,1,-1  !!processes KLX-3+1 elements, from KLX-2 to 1
          ISPBLOC=(JLB-2)*ISIZESP+IRMN
#ifdef USE_OPENACC
          !$ACC LOOP VECTOR PRIVATE(JLEVTOT,JLEV,IOFFSET1,JI,ISE)
#endif
          DO JSPTOT=ISPBLOC+ISIZESP+2,MAX(ISPBLOC+1,1),-1   !!pre-loads parts of the data in shared memory
            ISE=ISTA+2*(JSPTOT-1)
            DO JLEV=1,MIN(ISIZELEV,YDGEOMETRY%YRDIMV%NFLEVG-(JLEVBLOC-1)*ISIZELEV)
              JLEVTOT=JLEV+(JLEVBLOC-1)*ISIZELEV
              IOFFSET1=SIHEGSTA-1+(JLEVTOT-1)*(YDGEOMETRY%YRDIM%NSMAX+1-IM)
              ZPAS(JSPTOT-ISPBLOC,JLEV)=SIHEG(IOFFSET1+JSPTOT,1)
              ZPBS(JSPTOT-ISPBLOC,JLEV)=SIHEG(IOFFSET1+JSPTOT,2)
              ZPCS(JSPTOT-ISPBLOC,JLEV)=SIHEG(IOFFSET1+JSPTOT,3)
              DO JI=1,2
                ZPIN(JSPTOT-ISPBLOC,JLEV,JI)=PSPDIVP(ISE+JI-1,JLEVTOT)
              ENDDO
            ENDDO
          ENDDO
          IF (JLB==(KLX-3)/ISIZESP+1) THEN   !!processes the pre-loaded data, if last block
#ifdef USE_OPENACC

             !$ACC LOOP VECTOR PRIVATE(JSP) COLLAPSE(2)
#endif
             DO JLEV=1,MIN(ISIZELEV,YDGEOMETRY%YRDIMV%NFLEVG-(JLEVBLOC-1)*ISIZELEV)
               DO JI=1,2
                 ZPIN(KLX-ISPBLOC,JLEV,JI)=ZPIN(KLX-ISPBLOC,JLEV,JI)
                 ZPIN(KLX-1-ISPBLOC,JLEV,JI)=ZPIN(KLX-1-ISPBLOC,JLEV,JI)&
                  &-ZPBS(KLX-1-ISPBLOC,JLEV)/ZPAS(KLX-1-ISPBLOC,JLEV)*ZPIN(KLX-ISPBLOC,JLEV,JI)
                 DO JSP=ISIZESP,MAX(ISPBLOC+1,1)-ISPBLOC,-1
                   ZPIN(JSP,JLEV,JI)=ZPIN(JSP,JLEV,JI)-ZPBS(JSP,JLEV)/ZPAS(JSP,JLEV)*ZPIN(JSP+1,JLEV,JI)&
                    &-ZPCS(JSP,JLEV)/ZPAS(JSP,JLEV)*ZPIN(JSP+2,JLEV,JI)     
                 ENDDO
               ENDDO
             ENDDO
          ELSE                              !!processes the pre-loaded data, if preceding blocks
             ISE=ISTA+2*(ISPBLOC+ISIZESP+1-1)
#ifdef USE_OPENACC

             !$ACC LOOP VECTOR PRIVATE(JLEVTOT,JSP) COLLAPSE(2)
#endif
             DO JLEV=1,MIN(ISIZELEV,YDGEOMETRY%YRDIMV%NFLEVG-(JLEVBLOC-1)*ISIZELEV)
               DO JI=1,2
                 JLEVTOT=JLEV+(JLEVBLOC-1)*ISIZELEV
                 ZPIN(ISIZESP+1,JLEV,JI)=PSPDIVP(ISE+JI-1,JLEVTOT)
                 ZPIN(ISIZESP+2,JLEV,JI)=PSPDIVP(ISE+JI+1,JLEVTOT)
                 DO JSP=ISIZESP,MAX(ISPBLOC+1,1)-ISPBLOC,-1
                   ZPIN(JSP,JLEV,JI)=ZPIN(JSP,JLEV,JI)-ZPBS(JSP,JLEV)/ZPAS(JSP,JLEV)*ZPIN(JSP+1,JLEV,JI)&
                   &-ZPCS(JSP,JLEV)/ZPAS(JSP,JLEV)*ZPIN(JSP+2,JLEV,JI)     
                 ENDDO
               ENDDO
             ENDDO
          ENDIF
#ifdef USE_OPENACC
          !$ACC LOOP VECTOR PRIVATE(JLEV,JI,ISE)
#endif
          DO JSPTOT=ISPBLOC+ISIZESP+2,MAX(ISPBLOC+1,1),-1   !!moves the processed data from shared memory to main variable
            ISE=ISTA+2*(JSPTOT-1)
            DO JLEV=1,MIN(ISIZELEV,YDGEOMETRY%YRDIMV%NFLEVG-(JLEVBLOC-1)*ISIZELEV)
              DO JI=1,2
                PSPDIVP(ISE+JI-1,JLEV+(JLEVBLOC-1)*ISIZELEV)=ZPIN(JSPTOT-ISPBLOC,JLEV,JI)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ELSE                                                 !!direct processing, without shared memory, if small-size data
        ISE=ISTA+2*(KLX-1-1)
#ifdef USE_OPENACC
        !$ACC LOOP VECTOR PRIVATE(JLEVTOT,IOFFSET1) COLLAPSE(2)
#endif
        DO JLEV=1,MIN(ISIZELEV,YDGEOMETRY%YRDIMV%NFLEVG-(JLEVBLOC-1)*ISIZELEV)
          DO JI=1,2
            JLEVTOT=JLEV+(JLEVBLOC-1)*ISIZELEV
            IOFFSET1=SIHEGSTA-1+(JLEVTOT-1)*(YDGEOMETRY%YRDIM%NSMAX+1-IM)
            IF (KLX >= 2) THEN
              PSPDIVP(ISE+JI-1,JLEVTOT)=PSPDIVP(ISE+JI-1,JLEVTOT)&
               &-SIHEG(IOFFSET1+KLX-1,2)/SIHEG(IOFFSET1+KLX-1,1)*PSPDIVP(ISE+JI+1,JLEVTOT)
            ENDIF
          ENDDO
        ENDDO
      ENDIF


    ELSE               !!Inversion of a non-symmetrix matrix, case IM==0
      DO JI=1,2
 
        IF (JI==1) THEN
  
          !!solves L.(U.X)=Y for U.X
          IF (KLX >= 3) THEN
#ifdef USE_OPENACC

            !$ACC LOOP SEQ
#endif
            DO JLB=1,(KLX-3)/ISIZESP+1  !!processes KLX-3+1 elements, from 3 to KLX
              ISPBLOC=(JLB-1)*ISIZESP
#ifdef USE_OPENACC

              !$ACC LOOP VECTOR PRIVATE(JLEV,JLEVTOT,IOFFSET1,ISE)
#endif
              DO JSPTOT=ISPBLOC+1,MIN(ISPBLOC+ISIZESP+2,KLX)      !!pre-loads parts of the data in shared memory
                ISE=ISTA+2*(JSPTOT-1)
                DO JLEV=1,MIN(ISIZELEV,YDGEOMETRY%YRDIMV%NFLEVG-(JLEVBLOC-1)*ISIZELEV)
                  JLEVTOT=JLEV+(JLEVBLOC-1)*ISIZELEV
                  IOFFSET1=SIHEGSTA-1+(JLEVTOT-1)*(YDGEOMETRY%YRDIM%NSMAX+1-IM)
                  ZPAS(JSPTOT-ISPBLOC,JLEV)=SIHEG(IOFFSET1+JSPTOT,1)
                  ZPBS(JSPTOT-ISPBLOC,JLEV)=SIHEG(IOFFSET1+JSPTOT,2)
                  ZPCS(JSPTOT-ISPBLOC,JLEV)=SIHEG(IOFFSET1+JSPTOT,3)
                  ZPIN(JSPTOT-ISPBLOC,JLEV,JI)=PSDIVP(ISE+JI-1,JLEVTOT)
                ENDDO
              ENDDO


              IF (JLB==1) THEN        !!processes the pre-loaded data, if first block
#ifdef USE_OPENACC

                !$ACC LOOP VECTOR PRIVATE(JSP)
#endif
                DO JLEV=1,MIN(ISIZELEV,YDGEOMETRY%YRDIMV%NFLEVG-(JLEVBLOC-1)*ISIZELEV)
                  ZPIN(1,JLEV,JI)=ZPIN(1,JLEV,JI)/ZPAS(1,JLEV)
                  ZPIN(2,JLEV,JI)=(ZPIN(2,JLEV,JI)-ZPBS(1,JLEV)*ZPIN(1,JLEV,JI))/ZPAS(2,JLEV)
                  DO JSP=3,MIN(ISPBLOC+ISIZESP+2,KLX)-ISPBLOC
                    ZPIN(JSP,JLEV,JI)=(ZPIN(JSP,JLEV,JI)-ZPCS(JSP-2,JLEV)*ZPIN(JSP-2,JLEV,JI)&
                     & -ZPBS(JSP-1,JLEV)*ZPIN(JSP-1,JLEV,JI))/ZPAS(JSP,JLEV)  
                  ENDDO
                ENDDO
              ELSE                    !!processes the pre-loaded data, if following blocks
                ISE=ISTA+2*(ISPBLOC+1-1)
#ifdef USE_OPENACC

                !$ACC LOOP VECTOR PRIVATE(JLEVTOT,JSP)
#endif
                DO JLEV=1,MIN(ISIZELEV,YDGEOMETRY%YRDIMV%NFLEVG-(JLEVBLOC-1)*ISIZELEV)
                  JLEVTOT=JLEV+(JLEVBLOC-1)*ISIZELEV
                  ZPIN(1,JLEV,JI)=PSPDIVP(ISE+JI-1,JLEVTOT)
                  ZPIN(2,JLEV,JI)=PSPDIVP(ISE+JI+1,JLEVTOT)
                  DO JSP=3,MIN(ISPBLOC+ISIZESP+2,KLX)-ISPBLOC
                    ZPIN(JSP,JLEV,JI)=(ZPIN(JSP,JLEV,JI)-ZPCS(JSP-2,JLEV)*ZPIN(JSP-2,JLEV,JI)&
                     & -ZPBS(JSP-1,JLEV)*ZPIN(JSP-1,JLEV,JI))/ZPAS(JSP,JLEV)  
                  ENDDO
                ENDDO
              ENDIF
#ifdef USE_OPENACC

              !$ACC LOOP VECTOR PRIVATE(JLEV,ISE)
#endif

              DO JSPTOT=ISPBLOC+1,MIN(ISPBLOC+ISIZESP+2,KLX)   !!moves the processed data from shared memory to main variable 
                ISE=ISTA+2*(JSPTOT-1)
                DO JLEV=1,MIN(ISIZELEV,YDGEOMETRY%YRDIMV%NFLEVG-(JLEVBLOC-1)*ISIZELEV)
                  PSPDIVP(ISE+JI-1,JLEV+(JLEVBLOC-1)*ISIZELEV)=ZPIN(JSPTOT-ISPBLOC,JLEV,JI)
                ENDDO
              ENDDO
            ENDDO
          ELSE                   !!direct processing, without shared memory, if small-size data
            ISE=ISTA+2*(1-1)
#ifdef USE_OPENACC

            !$ACC LOOP VECTOR PRIVATE(JLEVTOT,IOFFSET1)
#endif
            DO JLEV=1,MIN(ISIZELEV,YDGEOMETRY%YRDIMV%NFLEVG-(JLEVBLOC-1)*ISIZELEV)
              JLEVTOT=JLEV+(JLEVBLOC-1)*ISIZELEV
              IOFFSET1=SIHEGSTA-1+(JLEVTOT-1)*(YDGEOMETRY%YRDIM%NSMAX+1-IM)
              PSPDIVP(ISE+JI-1,JLEVTOT)=PSDIVP(ISE+JI-1,JLEVTOT)/SIHEG(IOFFSET1+1,1)
              IF (KLX >= 2) THEN
                PSPDIVP(ISE+JI+1,JLEVTOT)=(PSDIVP(ISE+JI+1,JLEVTOT)&
                 &-SIHEG(IOFFSET1+1,2)*PSPDIVP(ISE+JI-1,JLEVTOT))/SIHEG(IOFFSET1+2,1)
              ENDIF
            ENDDO
          ENDIF
 
 
          !!solves U.X=Linv.Y for X
          IF (KLX >= 3) THEN
            IRMN=MOD(KLX-3,ISIZESP)+1
#ifdef USE_OPENACC

            !$ACC LOOP SEQ
#endif
            DO JLB=(KLX-3)/ISIZESP+1,1,-1  !!processes KLX-3+1 elements, from KLX-2 to 1
              ISPBLOC=(JLB-2)*ISIZESP+IRMN
#ifdef USE_OPENACC

              !$ACC LOOP VECTOR PRIVATE(JLEV,JLEVTOT,IOFFSET1,ISE)
#endif
              DO JSPTOT=ISPBLOC+ISIZESP+2,MAX(ISPBLOC+1,1),-1   !!pre-loads parts of the data in shared memory
                ISE=ISTA+2*(JSPTOT-1)
                DO JLEV=1,MIN(ISIZELEV,YDGEOMETRY%YRDIMV%NFLEVG-(JLEVBLOC-1)*ISIZELEV)
                  JLEVTOT=JLEV+(JLEVBLOC-1)*ISIZELEV
                  ZPBS(JSPTOT-ISPBLOC,JLEV)=SIHEG2(JSPTOT,JLEVTOT,2)
                  ZPCS(JSPTOT-ISPBLOC,JLEV)=SIHEG2(JSPTOT,JLEVTOT,3)
                  ZPIN(JSPTOT-ISPBLOC,JLEV,JI)=PSPDIVP(ISE+JI-1,JLEVTOT)
                ENDDO
              ENDDO
              IF (JLB==(KLX-3)/ISIZESP+1) THEN    !!processes the pre-loaded data, if last block
#ifdef USE_OPENACC

                !$ACC LOOP VECTOR PRIVATE(JSP)
#endif
                DO JLEV=1,MIN(ISIZELEV,YDGEOMETRY%YRDIMV%NFLEVG-(JLEVBLOC-1)*ISIZELEV)
                  ZPIN(KLX-ISPBLOC,JLEV,JI)=ZPIN(KLX-ISPBLOC,JLEV,JI)
                  ZPIN(KLX-1-ISPBLOC,JLEV,JI)=ZPIN(KLX-1-ISPBLOC,JLEV,JI)&
                    &-ZPBS(KLX-1-ISPBLOC,JLEV)*ZPIN(KLX-ISPBLOC,JLEV,JI)
                  DO JSP=ISIZESP,MAX(ISPBLOC+1,1)-ISPBLOC,-1
                    ZPIN(JSP,JLEV,JI)=ZPIN(JSP,JLEV,JI)-ZPBS(JSP,JLEV)*ZPIN(JSP+1,JLEV,JI)&
                      &-ZPCS(JSP,JLEV)*ZPIN(JSP+2,JLEV,JI)
                  ENDDO
                ENDDO
              ELSE                                !!processes the pre-loaded data, if preceding blocks
                ISE=ISTA+2*(ISPBLOC+ISIZESP+1-1)
#ifdef USE_OPENACC

                !$ACC LOOP VECTOR PRIVATE(JLEVTOT,JSP)
#endif
                DO JLEV=1,MIN(ISIZELEV,YDGEOMETRY%YRDIMV%NFLEVG-(JLEVBLOC-1)*ISIZELEV)
                  JLEVTOT=JLEV+(JLEVBLOC-1)*ISIZELEV
                  ZPIN(ISIZESP+1,JLEV,JI)=PSPDIVP(ISE+JI-1,JLEVTOT)
                  ZPIN(ISIZESP+2,JLEV,JI)=PSPDIVP(ISE+JI+1,JLEVTOT)
                  DO JSP=ISIZESP,MAX(ISPBLOC+1,1)-ISPBLOC,-1
                    ZPIN(JSP,JLEV,JI)=ZPIN(JSP,JLEV,JI)-ZPBS(JSP,JLEV)*ZPIN(JSP+1,JLEV,JI)&
                      &-ZPCS(JSP,JLEV)*ZPIN(JSP+2,JLEV,JI)
                  ENDDO
                ENDDO
              ENDIF
#ifdef USE_OPENACC

              !$ACC LOOP VECTOR PRIVATE(JLEV)
#endif
              DO JSPTOT=ISPBLOC+ISIZESP+2,MAX(ISPBLOC+1,1),-1   !!moves the processed data from shared memory to main variable
                ISE=ISTA+2*(JSPTOT-1)
                DO JLEV=1,MIN(ISIZELEV,YDGEOMETRY%YRDIMV%NFLEVG-(JLEVBLOC-1)*ISIZELEV)
                  PSPDIVP(ISE+JI-1,JLEV+(JLEVBLOC-1)*ISIZELEV)=ZPIN(JSPTOT-ISPBLOC,JLEV,JI)
                ENDDO
              ENDDO
            ENDDO
          ELSE                                                 !!direct processing, without shared memory, if small-size data
            ISE=ISTA+2*(KLX-1-1)
#ifdef USE_OPENACC

            !$ACC LOOP VECTOR PRIVATE(JLEVTOT,IOFFSET1)
#endif
            DO JLEV=1,MIN(ISIZELEV,YDGEOMETRY%YRDIMV%NFLEVG-(JLEVBLOC-1)*ISIZELEV)
              JLEVTOT=JLEV+(JLEVBLOC-1)*ISIZELEV
              IOFFSET1=(JLEVTOT-1)*(YDGEOMETRY%YRDIM%NSMAX+1-IM)
              IF (KLX >= 2) THEN
                PSPDIVP(ISE+JI-1,JLEVTOT)=PSPDIVP(ISE+JI-1,JLEVTOT)&
                   &-SIHEG2(KLX-1,JLEVTOT,2)*PSPDIVP(ISE+JI+1,JLEVTOT)
              ENDIF
            ENDDO
          ENDIF
  
        ELSE           !!case IM==0, JI==2

      !!for KM=0 values with JI=2 are set to 0, as was the case in the original code
      !!(in the 49t0 version of SPCSI, array ZSPDIVG was initialized to zero, and the
      !!values with JI=2 for KM=0 were not computed)
#ifdef USE_OPENACC

          !$ACC LOOP VECTOR PRIVATE(JLEVTOT,JLEV)
#endif
          DO JSP=1,YDGEOMETRY%YRDIM%NSMAX+1-IM
            ISE=ISTA+2*(JSP-1)
            DO JLEV=1,MIN(ISIZELEV,YDGEOMETRY%YRDIMV%NFLEVG-(JLEVBLOC-1)*ISIZELEV)
              JLEVTOT=JLEV+(JLEVBLOC-1)*ISIZELEV
              PSPDIVP(ISE+JI-1,JLEVTOT)=0.0_JPRB
            ENDDO
          ENDDO
        ENDIF
      ENDDO !!JI loop, case IM==0

    ENDIF   !!IM>0 or IM==0

  ENDDO  !!JLEVTOT loop
ENDDO    !!JMLOC loop
#ifdef USE_OPENACC
!$ACC END PARALLEL
#endif
#ifdef USE_OPENMP
!$OMP END TARGET TEAMS DISTRIBUTE
#endif

#ifdef USE_OPENACC
!$ACC END DATA
#endif

#ifdef USE_OPENMP
!$OMP END TARGET DATA
#endif

IF (LHOOK) CALL DR_HOOK('SPCSIDG_PART1ACC',1,ZHOOK_HANDLE)
END SUBROUTINE SPCSIDG_PART1ACC

