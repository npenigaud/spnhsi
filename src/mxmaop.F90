SUBROUTINE MXMAOP(PA,KA,KAD,PB,KB,KBD,PC,KC,KCA,KAR,KAC,KBC,LDACC)

!**** *MXMAOP - Optimize call to SGEMMX

!     Purpose.
!     --------
!        Make sure SGEMMX is called in a way to insure maximum optimization.
!        Does matrix product PC=PA*PB

!**   Interface.
!     ----------
!        CALL MXMAOP(PA,KA,KAD,PB,KB,KBD,PC,KC,KCA,KAR,KAC,KBC)

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

!     Method.
!     -------

!     Externals.   SGEMMX in Cray library.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 88-01-28
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        K. Yessad (Apr 2009): add comments.
!        J. Hague  (Oct 2012): Parallelise call to SGEMMX
!        R. El Khatib 03-Aug-2023 Workaround for bitwise reproducibility in single precision on NEC
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE OML_MOD   ,ONLY : OML_MAX_THREADS, OML_IN_PARALLEL

!     ------------------------------------------------------------------

IMPLICIT NONE

REAL(KIND=JPRB)   ,INTENT(IN)    :: PA(*) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KAD 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PB(*) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KB 
INTEGER(KIND=JPIM),INTENT(IN)    :: KBD 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PC(*) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KC 
INTEGER(KIND=JPIM),INTENT(IN)    :: KCA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KAR 
INTEGER(KIND=JPIM),INTENT(IN)    :: KAC 
INTEGER(KIND=JPIM),INTENT(IN)    :: KBC
LOGICAL, OPTIONAL, INTENT(IN)    :: LDACC
character :: clenv2
integer(kind=jpim) :: ilenv2
!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JJJ,JLEN,JT,JCHUNK
LOGICAL            :: LLACC

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MXMAOP',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------
call get_environment_variable("LLSIMPLE_DGEMM",clenv2,ilenv2)
if (ilenv2==0) clenv2='0'
LLACC=.FALSE.
IF (PRESENT(LDACC)) LLACC=LDACC
if (clenv2=='1') then
  if (llacc) then
    call simple_dgemm3(pc,pa,kar,kac,pb,kbc,kac,llacc)
  else
    call simple_dgemm2(pc,pa,kar,kac,pb,kac,kbc,llacc)
!    call simple_dgemm2(pc,pa,kad,kbd,pb,kbd,kbd,llacc)
  endif
else
IF (LLACC) THEN
!   write (0,*) "branche mxmaop llacc"
   CALL SGEMMXACC(PA,KA,KAD,PB,KB,KBD,PC,KC,KCA,KAR,KAC,KBC)

ELSE

!*       1.       PERFORM LEGENDRE TRANFORM.
!                 --------------------------

#ifdef __NEC__
! workaround for bitwise reproducibility : sequential call in single precision
  IF(OML_IN_PARALLEL() .OR. .NOT.(JPRB==JPRD)) THEN
#else
  IF(OML_IN_PARALLEL()) THEN
#endif

    IF (KAR >= KBC) THEN
      CALL SGEMMX(KAR,KBC,KAC,1.0_JPRB,PA,KA,KAD,PB,KB,KBD,0.0_JPRB,PC,KC,KCA)
    ELSE
      CALL SGEMMX(KBC,KAR,KAC,1.0_JPRB,PB,KBD,KB,PA,KAD,KA,0.0_JPRB,PC,KCA,KC)
    ENDIF
  
  ELSE
  
    JT=OML_MAX_THREADS()
  
    IF (KAR >= KBC) THEN
      JCHUNK=(KAR-1)/JT+1
  !$OMP PARALLEL DO PRIVATE(JJJ,JLEN)
      DO JJJ=1,KAR,JCHUNK
        JLEN=MIN(JCHUNK,KAR-JJJ+1)
         IF(JLEN>0)  CALL SGEMMX(JLEN,KBC,KAC,1.0_JPRB,PA((JJJ-1)*KA+1),KA,KAD,PB,KB,KBD,0.0_JPRB,PC((JJJ-1)*KC+1),KC,KCA)
      ENDDO
  !$OMP END PARALLEL DO 
  
    ELSE
  
      JCHUNK=(KBC-1)/JT+1
  !$OMP PARALLEL DO PRIVATE(JJJ,JLEN)
      DO JJJ=1,KBC,JCHUNK
        JLEN=MIN(JCHUNK,KBC-JJJ+1)
        IF(JLEN>0) CALL SGEMMX(JLEN,KAR,KAC,1.0_JPRB,PB((JJJ-1)*KBD+1),KBD,KB,PA,KAD,KA,0.0_JPRB,PC((JJJ-1)*KCA+1),KCA,KC)
      ENDDO
  !$OMP END PARALLEL DO 
    ENDIF
  
  ENDIF

ENDIF
endif
!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('MXMAOP',1,ZHOOK_HANDLE)
contains

subroutine simple_dgemm2(pab,pa,ligneA,colA,PB,ligneB,colB,ldacc)
integer(kind=jpim), intent(in) :: ligneA,colA,ligneB,colB
logical, intent(in) :: ldacc
real(kind=jprb), intent(out) :: pab(ligneA,colB)
real(kind=jprb), intent(in) :: pa(ligneA,colA)
real(kind=jprb), intent(in) :: pb(ligneB,colB)

integer(kind=jpim) :: ii,ij,ik
integer(kind=jpim) :: ni,nj,nk

do ik=1,ligneA
  do ii=1,colB
    pab(ik,ii)=0._jprb
    do ij=1,colA
      pab(ik,ii)=pab(ik,ii)+pa(ik,ij)*pb(ij,ii)
    end do
  end do
end do

end subroutine simple_dgemm2

subroutine simple_dgemm3(pab,pa,ligneA,colA,PB,ligneB,colB,ldacc)
integer(kind=jpim), intent(in) :: ligneA,colA,ligneB,colB
logical, intent(in) :: ldacc
real(kind=jprb), intent(out) :: pab(ligneB,ligneA)
real(kind=jprb), intent(in) :: pa(ligneA,colA)
real(kind=jprb), intent(in) :: pb(ligneB,colB)

integer(kind=jpim) :: ii,ij,ik

#ifdef USE_OPENACC
!$acc parallel loop gang vector collapse(2) private(ij) present(pab,pa,pb) !if(ldacc)
#endif
#ifdef USE_OPENMP
!$omp target teams distribute parallel do simd private(ij) collapse(2) map(present:pab,pa,pb) if(ldacc)
#endif
do ik=1,ligneB
  do ii=1,ligneA
    pab(ik,ii)=0._jprb
    do ij=1,colB
      pab(ik,ii)=pab(ik,ii)+pb(ik,ij)*pa(ii,ij)
    end do
  end do
end do
#ifdef USE_OPENMP
!$omp end target teams distribute parallel do simd
#endif

end subroutine simple_dgemm3

END SUBROUTINE MXMAOP

