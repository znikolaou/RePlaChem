      PROGRAM TEST_GET_CHARGED_SPECIES
      USE GLOBAL 
      IMPLICIT NONE
      INTEGER :: I
      CHARACTER(LEN=*), PARAMETER :: FL='../input_test/species.txt'

      CALL SET_SPECIES(FL)
      CALL SET_IS_NEUTRAL()

      DO I=1,NSPEC
       WRITE(*,*) I,TRIM(ADJUSTL(SPEC(I))),'   IS NEUTRAL:', &
               IS_SPEC_NEUTRAL(I)
      ENDDO

      CLOSE(1)
      STOP
      END
