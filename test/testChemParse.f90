      PROGRAM TEST_CHEMPARSE
      USE CHEM_PARSE
      IMPLICIT NONE
      INTEGER :: I
      CHARACTER(LEN=*), PARAMETER :: &
      FLA='./input_chemistry/argon2step/kinet2stepArgon.inp'
       
      CALL CM_INIT(FLA)

      STOP
      END
