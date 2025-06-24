      !-----------------------------------------------------------------
      !
      ! AUTHOR: Z. NIKOLAOU
      !
      !-----------------------------------------------------------------
      PROGRAM TEST_PCMO
      USE PCMO
      IMPLICIT NONE
      INTEGER, PARAMETER :: N=100
      INTEGER :: I,IC,NR,DATA_COUNT
      DOUBLE PRECISION :: TIME,A(N),AR(N)

      DO I=1,N
       A(I)=I
      ENDDO

      CALL PCMO_WRITE_BIN('./output/','reactionRates',1,10,20.0d0,N,A)

      OPEN(UNIT=1,FILE='./output/reactionRates_00000001.dat', &
           STATUS='OLD',FORM='UNFORMATTED')
      READ(1) IC,DATA_COUNT,TIME,NR,AR

      WRITE(*,*) IC,DATA_COUNT,TIME,NR
      DO I=1,NR
       WRITE(*,*) AR(I)
      ENDDO
      
      STOP
      END
      !-----------------------------------------------------------------

