PROGRAM TESTBATCHED

USE FXTRAN_ACDC_GEMM_MOD_MANYBLOCKS, ONLY : FXTRAN_ACDC_GEMM_MANYBLOCKS

IMPLICIT NONE

INTEGER(4), PARAMETER :: KPROMA=40 !4
INTEGER(4), PARAMETER :: KLEVIN=12 !3
INTEGER(4), PARAMETER :: KLEVOUT=11 !3
INTEGER(4), PARAMETER :: KGPBLKS=50 !2

INTEGER(4) :: JBLK,JLEVIN,JLEVOUT,JROF

REAL(8) :: ZIN (KPROMA,KLEVIN,KGPBLKS)
REAL(8) :: ZOUT(KPROMA,KLEVOUT-1,KGPBLKS)
REAL(8) :: ZOUT_MAIN(KPROMA,KLEVOUT-1,KGPBLKS)
REAL(8) :: ZINTE(KLEVOUT,KLEVIN)

INTEGER(4)        :: IRESU
REAL(8),PARAMETER :: ALPHA=1.0
REAL(8),PARAMETER :: BETA =0.0

REAL(8)   :: DELTA

LOGICAL :: VERBOSE,LDONE

CALL RANDOM_NUMBER(ZIN)
CALL RANDOM_NUMBER(ZINTE)
ZOUT(:,:,:)=0.
ZOUT_MAIN(:,:,:)=0.

VERBOSE=.FALSE.

IF (VERBOSE) THEN
  write (0,*) "initial data - zinte"
    do jlevout=1,klevout
      do jlevin=1,klevin
        write (0,*) "i j zinte",jlevout,jlevin,zinte(jlevout,jlevin)
      enddo
    enddo
  write (0,*) "initial data - zin"
    DO JBLK = 1, KGPBLKS
      DO JLEVIN = 1, KLEVIN
        DO JROF = 1, KPROMA
          write (0,*) "i j block zin ",jrof,jlevin,jblk, ZIN (JROF,JLEVIN,JBLK)
        ENDDO
      ENDDO
    ENDDO
ENDIF

DO JBLK = 1, KGPBLKS
  ZOUT_MAIN (:,:,JBLK) = 0.0
  DO JLEVIN = 1, KLEVIN
    DO JLEVOUT = 1, KLEVOUT-1
      DO JROF = 1, KPROMA
        ZOUT_MAIN (JROF,JLEVOUT,JBLK) = ZOUT_MAIN (JROF,JLEVOUT,JBLK) + ZINTE(JLEVOUT,JLEVIN) * ZIN (JROF,JLEVIN,JBLK)
      ENDDO
    ENDDO
  ENDDO
ENDDO 

IF (VERBOSE) THEN
  WRITE (0,*) "COMPUTATION ON THE CPU"
  do jblk=1,kgpblks
    write (0,*) "block number ",jblk
    do jlevout=1,klevout
      do jrof=1,kproma
        write (0,*) "i j value",jrof,jlevout,zout_main(jrof,jlevout,jblk)
      enddo
    enddo
  enddo
  
  write (0,*) "=========="
ENDIF

#if (__NVCOMPILER==1)
!$ACC DATA COPYIN(ZIN,ZINTE) COPYOUT(ZOUT)
#else
!$OMP TARGET DATA MAP(TO:ZIN,ZINTE) MAP(FROM:ZOUT)
#endif

CALL FXTRAN_ACDC_GEMM_MANYBLOCKS(1,KPROMA,'N','T',KPROMA,KLEVOUT-1,KLEVIN,1.0_8,ZIN,&
                                  &KPROMA,ZINTE,KLEVOUT,0.0_8,ZOUT,KPROMA,LDONE,.TRUE.,KGPBLKS)
#if (__NVCOMPILER==1)
!$ACC END DATA
#else
!$OMP END TARGET DATA
#endif

WRITE(0,*) "LDONE : ",LDONE

IF (VERBOSE) THEN
  WRITE(0,*) "multiplication result : ",iresu
  WRITE (0,*) "COMPUTATION ON THE GPU, BATCHED"
  do jblk=1,kgpblks
    write (0,*) "block number ",jblk
    do jlevout=1,klevout-1
      do jrof=1,kproma
        write (0,*) "i j value",jrof,jlevout,zout(jrof,jlevout,jblk)
      enddo
    enddo
  enddo
ENDIF

WRITE(0,*) "Difference initial multiplication / batched multiplication :"
  do jblk=1,kgpblks
    write (0,*) "block number ",jblk
    do jlevout=1,klevout-1
      do jrof=1,kproma
        delta=abs(zout(jrof,jlevout,jblk)-zout_main(jrof,jlevout,jblk))/(abs(zout_main(jrof,jlevout,jblk)))
        IF (delta .GE. 0.00000000000001) write (0,*) "i  delta ",jrof,delta
        IF (delta .GE. 0.00000000000001) write (0,*) "i j value value_main",jrof,jlevout,zout(jrof,jlevout,jblk),zout_main(jrof,jlevout,jblk)
      enddo
    enddo
  enddo


END PROGRAM TESTBATCHED
