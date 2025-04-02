      PROGRAM TEST_ZDPLASKIN
      USE ZDPLASKIN_PARSE
      IMPLICIT NONE
      INTEGER :: I
      CHARACTER(LEN=*), PARAMETER :: &
      FLA='./input_chemistry/argon2step/kinet2stepArgonModified.inp', &
      FLB='./input_chemistry/methaneAntwerp/'//&
              'kinetAntwerpMethaneNoTabs.inp', &
      FLC='./input_chemistry/chemkinFormatGri3/grimech3.txt'

      
      CALL ZDP_INIT(FLB)

      STOP
      END
