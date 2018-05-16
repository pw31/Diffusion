!======================================================================
      LOGICAL FUNCTION IS_NAN16(X)
!======================================================================
#ifdef IFORT
      USE,INTRINSIC :: IEEE_ARITHMETIC
#endif
      IMPLICIT NONE
      INTEGER,PARAMETER :: qp = selected_real_kind ( 33, 4931 )
      REAL(KIND=qp) :: X
#ifdef IFORT
      IS_NAN16 = (IEEE_IS_NAN(X).ne.0)
#else
      IS_NAN16 = (.not.(x.eq.x))
#endif
      RETURN
      END FUNCTION IS_NAN16
