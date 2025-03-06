      !-----------------------------------------------------------------
      SUBROUTINE READ_INT_AND_STR(FL,N,STRL,NL)
      USE GLOBAL, ONLY : NSMX
      CHARACTER(LEN=*) :: FL 
      INTEGER :: N,ID,STAT,I,NL,NA
      CHARACTER(LEN=NSMX) :: STRL(N),CL(2),C
      LOGICAL ISEMPTY,ISCOMMENT
      PARAMETER(ID=1)
      
      OPEN(UNIT=ID,FILE=TRIM(ADJUSTL(FL)),STATUS='OLD',FORM='FORMATTED')

      I=0
      STAT=0
      DO WHILE(STAT.EQ.0)
       READ(ID,'(A)',IOSTAT=STAT) C 
       IF(.NOT.ISEMPTY(C)) THEN
        IF(.NOT.ISCOMMENT(C)) THEN
         CALL SPLIT_STRING(C,' ',2,CL,NA)
          I=I+1
          STRL(I)=CL(2)
        ENDIF
       ENDIF
      ENDDO      
      NL=I

      CLOSE(ID)

      END SUBROUTINE 
      !-----------------------------------------------------------------
      SUBROUTINE REMOVE_TABS_FROM_FILE(FL)
      IMPLICIT NONE
      CHARACTER(LEN=*) :: FL
            
      CALL SYSTEM("sed -i 's/\t/\ /g'"//' '//FL) 
      CALL SYSTEM("sed -i 's/\t/\ /g'"//' '//FL)
       
      END
      !-----------------------------------------------------------------
      SUBROUTINE READ_LINES(FL,LINES,NL)
      USE GLOBAL, ONLY : NSMX_LINE,NLINEMX
      IMPLICIT NONE
      CHARACTER(LEN=*) FL
      CHARACTER(NSMX_LINE) :: L,LINES(NLINEMX)
      INTEGER, PARAMETER :: ID=1
      INTEGER :: IOS,I,NL

      WRITE(*,*) '***READ_LINES:***'
      WRITE(*,*) 'FILE:'
      WRITE(*,*) FL
      OPEN(UNIT=ID,FILE=FL,STATUS='OLD',FORM='FORMATTED',IOSTAT=IOS)

      IF(IOS.LT.0) THEN
       WRITE(*,*) 'ERROR READING FILE'
       STOP
      ENDIF

      NL=0
      DO
       READ(ID,'(A)',IOSTAT=IOS) L
       IF(IOS.NE.0) EXIT
       NL=NL+1
       IF(NL.GT.NLINEMX) THEN
        WRITE(*,*) 'MAX LINES EXCEEDED, TERMINATING ...'
        STOP
       ENDIF
      ENDDO
      WRITE(*,*) 'NO OF LINES:',NL

      REWIND(ID)
      DO I=1,NL
       READ(ID,'(A)') LINES(I)
      ENDDO

      CLOSE(ID)
            
      WRITE(*,*) '***READ_LINES:***'

      END
      !-----------------------------------------------------------------
