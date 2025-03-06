      PROGRAM TEST_READ_LINES
      USE GLOBAL, ONLY : NSMX_LINE,NLINEMX
      IMPLICIT NONE 
      CHARACTER(LEN=*), PARAMETER :: DIR='../input_chemistry/'
      CHARACTER(LEN=*), PARAMETER :: FL1= &
              DIR//'argon2step/kinet2stepArgon.inp'
      CHARACTER(LEN=*), PARAMETER :: FL2= &
              DIR//'methaneAntwerp/kinetAntwerpMethane.inp'

      CHARACTER(LEN=NSMX_LINE) :: LINES(NLINEMX)
      INTEGER NL,I

      CALL READ_LINES(FL2,LINES,NL)
      DO I=1,NL
       WRITE(*,'(I5XA)') I,TRIM(ADJUSTL(LINES(I)))
      ENDDO

      STOP
      END
