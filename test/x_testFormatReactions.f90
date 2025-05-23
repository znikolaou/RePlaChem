      PROGRAM TEST_FORMAT_REACTIONS
      USE GLOBAL 
      IMPLICIT NONE 
      INTEGER :: I
      CHARACTER(LEN=*), PARAMETER :: FL='../input_test/reactions.txt'
      
      CALL SET_REACTIONS(FL)
      CALL SET_FORMATTED_REACTIONS()

      DO I=1,NREAC
       WRITE(*,'(I4,X,A)') I,TRIM(ADJUSTL(REAC(I)))
       WRITE(*,'(X,A)') TRIM(ADJUSTL(REACF(I)))
      ENDDO

      CLOSE(1)
      STOP
      END
