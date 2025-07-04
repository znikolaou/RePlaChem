      !-----------------------------------------------------------------
      !
      ! AUTHOR: Z. NIKOLAOU
      !
      !-----------------------------------------------------------------
      SUBROUTINE GET_STATS(OIC,STATS)
      USE GLOBAL, ONLY : NTRG,NSPEC,NDATA,NCASE
      IMPLICIT NONE
      INTEGER :: I,J
      DOUBLE PRECISION :: OIC(NTRG,NSPEC),STATS(NTRG,NSPEC,3)

      DO J=1,NSPEC
       DO I=1,NTRG
        STATS(I,J,1)=OIC(I,J)+STATS(I,J,1)
        STATS(I,J,2)=MIN(OIC(I,J),STATS(I,J,2))
        STATS(I,J,3)=MAX(OIC(I,J),STATS(I,J,3))
       ENDDO
      ENDDO
      STATS(1:NTRG,1:NSPEC,1)=STATS(1:NTRG,1:NSPEC,1)/(NDATA*NCASE)

      RETURN
      END
      !-----------------------------------------------------------------
