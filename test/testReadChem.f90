      PROGRAM TEST_READ_CHEM
      USE GLOBAL
      IMPLICIT NONE 
      CHARACTER(LEN=*), PARAMETER :: DIR='./input_chemistry/'// &
                                         'methaneAntwerp/'
      CHARACTER(LEN=*), PARAMETER :: DIRB='./input_chemistry/'// &
                                         'chemkinFormatGri3//'


      INTEGER :: I,J

      CALL READ_CHEM(DIR,'kinetAntwerpMethaneNoTabs.inp')
      !CALL READ_CHEM(DIRB,'grimech3.txt')

      CLOSE(1)
      STOP
      END
