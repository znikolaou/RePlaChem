      PROGRAM TEST_SPLIT_STRING
      IMPLICIT NONE 
      INTEGER :: N,NA,NB,NC,I
      PARAMETER(N=100)
      CHARACTER(LEN=*) :: STRA,STRB,STRC
      CHARACTER(LEN=100) :: COLSA(N),COLSB(N),COLSC(N)
      PARAMETER(STRA='   A + B + C => D + E + F')
      PARAMETER(STRB='  A*B*C*-----123***F   ')
      PARAMETER(STRC='SOME!!!!TEXT!!!!WITH SEP')

      CALL SPLIT_STRING(STRA,' ',N,COLSA,NA)
      CALL SPLIT_STRING(STRB,'*',N,COLSB,NB)
      CALL SPLIT_STRING(STRC,'!!!!',N,COLSC,NC)

      WRITE(*,*) STRA,NA
      DO I=1,NA
       WRITE(*,*) COLSA(I)
      ENDDO
      WRITE(*,*) STRB,NB
      DO I=1,NB
       WRITE(*,*) COLSB(I)
      ENDDO
      WRITE(*,*) STRC,NC
      DO I=1,NC
       WRITE(*,*) COLSC(I)
      ENDDO

      STOP
      END
