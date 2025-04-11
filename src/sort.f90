      !-----------------------------------------------------------------
      !
      ! AUTHOR: Z. NIKOLAOU
      !
      !-----------------------------------------------------------------
      SUBROUTINE SORT(N,NC,R_ARRAY,C_ARRAY,SR_ARRAY,SC_ARRAY)
      IMPLICIT NONE

      INTEGER :: N,NC,I,J
      DOUBLE PRECISION :: R_ARRAY(N),SR_ARRAY(N),TEMPR
      CHARACTER(LEN=NC) :: C_ARRAY(N),SC_ARRAY(N),TEMPC
      !     
      SR_ARRAY(1:N)=R_ARRAY(1:N)
      SC_ARRAY(1:N)=C_ARRAY(1:N)

      DO J=2,N 
       TEMPR=SR_ARRAY(J) 
       TEMPC=SC_ARRAY(J)
       I=J-1       
       DO WHILE( (I.GE.1).AND.(SR_ARRAY(I).LT.TEMPR) ) 
        SR_ARRAY(I+1)=SR_ARRAY(I)
        SC_ARRAY(I+1)=SC_ARRAY(I)
        I=I-1
        IF(I.EQ.0) EXIT
       ENDDO
        SR_ARRAY(I+1)=TEMPR
        SC_ARRAY(I+1)=TEMPC
      ENDDO

      END
      !-----------------------------------------------------------------
