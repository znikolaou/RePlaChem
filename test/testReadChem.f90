      PROGRAM TEST_READ_CHEM
      USE GLOBAL
      IMPLICIT NONE 
      CHARACTER(LEN=*), PARAMETER :: DIR='../input_test/'
      INTEGER :: I,J

      CALL READ_CHEM(DIR,'reactions.txt','species.txt')

      CLOSE(1)
      STOP
      END
