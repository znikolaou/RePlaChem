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
      USE GLOBAL, ONLY : NSMX,NLINEMX
      IMPLICIT NONE
      CHARACTER(LEN=*) FL
      CHARACTER(NSMX) :: L,LINES(NLINEMX)
      INTEGER, PARAMETER :: ID=1
      INTEGER :: IOS,I,NL

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
      WRITE(*,'(AXI6)') 'NO OF LINES:',NL
      
      REWIND(ID)
      DO I=1,NL
       READ(ID,'(A)') LINES(I)
      ENDDO

      CLOSE(ID)
            
      END
      !-----------------------------------------------------------------
      SUBROUTINE READ_CONTROL()
      USE GLOBAL, ONLY: NTRG,NCASE,NDATA,NSPMX,ZERO,SPEC,CONTROL_FL, &
       INDX_TRG,ETOL,CHEMFL
      IMPLICIT NONE
      INTEGER :: I
      INTEGER, PARAMETER:: IU=1
      
      WRITE(*,*) '***READ_CONTROL***'
      
      OPEN(UNIT=IU,FILE=CONTROL_FL,STATUS='OLD',FORM='FORMATTED')
      
      READ(IU,*) 
      READ(IU,*) CHEMFL
      READ(IU,*) 
      READ(IU,*) NCASE,NDATA
      READ(IU,*)
      READ(IU,*) NTRG
      READ(IU,*)
      INDX_TRG(1:NSPMX)=0
      ETOL(1:NSPMX)=ZERO
      DO I=1,NTRG
       READ(IU,*) INDX_TRG(I),ETOL(I)      
      ENDDO
      WRITE(*,*) 'CHEMICAL MECHANISM FILE:',TRIM(ADJUSTL(CHEMFL))
      WRITE(*,*) 'NO TARGETS SPECIES, NO CASES, DATA SIZE'
      WRITE(*,*) NTRG,NCASE,NDATA
      WRITE(*,*) 'TARGET SPECIES INDEX, TARGET SPECIES, TOL:'
      DO I=1,NTRG
       WRITE(*,*) I,INDX_TRG(I),ETOL(I)       
      ENDDO

      CLOSE(IU)
 
      WRITE(*,*) '***READ_CONTROL***'
                 
      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE READ_RATES(IDAT,DIR,RNM,NSPEC,NREAC,WIJ,RR,JIJ)
      USE GLOBAL, ONLY : NSMX
      IMPLICIT NONE
      INTEGER I,J,IDAT,NSPEC,NREAC
      CHARACTER(LEN=*) :: DIR,RNM
      DOUBLE PRECISION :: WIJ(NREAC,NSPEC),RR(NREAC),JIJ(NSPEC,NSPEC)
      CHARACTER(LEN=4) :: CDAT
      CHARACTER(LEN=NSMX) :: FLNM
      INTEGER :: NSP,NRE
      DOUBLE PRECISION :: TIME,THETA
 
      !BUILD FILENAMES
      WRITE(CDAT,'(I4.4)') IDAT     
      FLNM=DIR//RNM//'_'//CDAT//'.dat'

      OPEN(UNIT=1,FILE=FLNM,STATUS='OLD',FORM='UNFORMATTED')
      READ(1) NSP,NRE
      CLOSE(1)
      
      IF(NSP.NE.NSPEC.OR.NRE.NE.NREAC) THEN
       WRITE(*,*) 'READ_RATES:ERROR, MISMATCH',NSP,NSPEC,NRE,NREAC
       WRITE(*,*) 'TERMINATING ...'
       STOP
      ELSE
       OPEN(UNIT=1,FILE=FLNM,STATUS='OLD',FORM='UNFORMATTED')
       !TODO READ SPECIES RATES MATRIX DIRECTLY
       READ(1) NSP,NRE,TIME,THETA,WIJ,RR,JIJ
       CLOSE(1)
      ENDIF
           
      END SUBROUTINE
      !-----------------------------------------------------------------
      SUBROUTINE READ_REACTION_RATES(IDAT,DIR,NREAC,RR)
      USE GLOBAL, ONLY : NSMX,REAC_RATE_FL
      IMPLICIT NONE
      INTEGER I,J,IDAT,NREAC,TSTEP
      CHARACTER(LEN=*) :: DIR
      DOUBLE PRECISION :: RR(NREAC)
      CHARACTER(LEN=8) :: CDAT
      CHARACTER(LEN=NSMX) :: FLNM
      INTEGER :: NSP,NRE,COUNT
      DOUBLE PRECISION :: TIME
 
      !BUILD FILENAMES
      WRITE(CDAT,'(I8.8)') IDAT     
      FLNM=TRIM(ADJUSTL(DIR))//REAC_RATE_FL//'_'//CDAT//'.dat'

      OPEN(UNIT=1,FILE=FLNM,STATUS='OLD',FORM='UNFORMATTED')
      READ(1) COUNT,TSTEP,TIME,NRE
      CLOSE(1)
      
      IF(NRE.NE.NREAC) THEN
       WRITE(*,*) 'READ_RATES:ERROR, MISMATCH',NRE,NREAC
       WRITE(*,*) 'TERMINATING ...'
       STOP
      ELSE
       OPEN(UNIT=1,FILE=FLNM,STATUS='OLD',FORM='UNFORMATTED')
       READ(1) COUNT,TSTEP,TIME,NRE,RR
       CLOSE(1)
      ENDIF
           
      END SUBROUTINE
      !-----------------------------------------------------------------
      SUBROUTINE READ_SPECIES_RATES_MATRIX(IDAT,DIR,RNM,NSPEC,NREAC,WIJ)
      USE GLOBAL, ONLY : NSMX
      IMPLICIT NONE
      INTEGER I,J,IDAT,NSPEC,NREAC
      CHARACTER(LEN=*) :: DIR,RNM
      DOUBLE PRECISION :: WIJ(NREAC,NSPEC)
      CHARACTER(LEN=4) :: CDAT
      CHARACTER(LEN=NSMX) :: FLNM
      INTEGER :: NSP,NRE
      DOUBLE PRECISION :: TIME
 
      !BUILD FILENAMES
      WRITE(CDAT,'(I4.4)') IDAT     
      FLNM=TRIM(ADJUSTL(DIR))//RNM
      
      WRITE(*,*) 'READING RATES FILE:',TRIM(ADJUSTL(FLNM))
      OPEN(UNIT=1,FILE=FLNM,STATUS='OLD',FORM='UNFORMATTED')
      READ(1) NRE,NSP
      CLOSE(1)
      
      IF(NSP.NE.NSPEC.OR.NRE.NE.NREAC) THEN
       WRITE(*,*) 'READ_RATES:ERROR, MISMATCH',NSP,NSPEC,NRE,NREAC
       WRITE(*,*) 'TERMINATING ...'
       STOP
      ELSE
       OPEN(UNIT=1,FILE=FLNM,STATUS='OLD',FORM='UNFORMATTED')
       READ(1) NRE,NSP,WIJ
       CLOSE(1)
      ENDIF
           
      END SUBROUTINE
      !-----------------------------------------------------------------
      FUNCTION BUILD_CASE_DIR(ICASE)
      USE GLOBAL, ONLY: INDIR,RATE_DIR,NSFLMX
      IMPLICIT NONE
      INTEGER :: ICASE
      CHARACTER(LEN=NSFLMX) :: BUILD_CASE_DIR
      CHARACTER(LEN=6) :: CHI

      WRITE(CHI,'(I6.6)') ICASE
      BUILD_CASE_DIR=INDIR//RATE_DIR//'case_'//CHI//'/'

      END FUNCTION
      !----------------------------------------------------------------- 
