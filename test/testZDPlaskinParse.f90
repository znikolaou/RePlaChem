      PROGRAM TEST_ZDPLASKIN
      USE ZDPLASKIN_PARSE
      IMPLICIT NONE
      INTEGER :: I
      CHARACTER(LEN=*), PARAMETER :: &
       FL='../input_chemistry/argon2step/kinet2stepArgon.inp', &
       FLA='../input_chemistry/methaneAntwerp/kinetAntwerpMethane.inp'

      CALL ZDP_INIT(FLA)

      STOP
      END
