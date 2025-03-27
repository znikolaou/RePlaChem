      PROGRAM TEST_READ_CHEM
      USE GLOBAL
      IMPLICIT NONE 
      CHARACTER(LEN=*), PARAMETER :: DIR='./input_chemistry/'// &
                                         'methaneAntwerp/'
      INTEGER :: I,J

      CALL READ_CHEM(DIR,'kinetAntwerpMethaneNoTabs.inp','')

      CLOSE(1)
      STOP
      END
