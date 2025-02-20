      PROGRAM TEST_SPLIT_STRING
      IMPLICIT NONE 
      INTEGER :: N,NA,NB,I
      PARAMETER(N=100)
      CHARACTER(LEN=*) :: STRA,STRB
      CHARACTER(LEN=100) :: COLSA(N),COLSB(N)
      PARAMETER(STRA='A + B + C => D + E + F')
      PARAMETER(STRB=' A*B*C*-----123***F')

      CALL SPLIT_STRING(STRA,' ',N,COLSA,NA)
      CALL SPLIT_STRING(STRB,'*',N,COLSB,NB)

      WRITE(*,*) STRA,NA
      DO I=1,NA
       WRITE(*,*) COLSA(I)
      ENDDO
      WRITE(*,*) STRB,NB
      DO I=1,NB
       WRITE(*,*) COLSB(I)
      ENDDO

      STOP
      END
