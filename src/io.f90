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

      WRITE(*,*) '***READ_LINES:***'
      WRITE(*,'(AXA)') 'FILE:',TRIM(ADJUSTL(FL))

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
            
      WRITE(*,*) '***READ_LINES:***'

      END
      !-----------------------------------------------------------------
      SUBROUTINE READ_CONTROL(NSPMX,NTRG,NCASE,NDATA,IRDMETH, &
                              INDX_TRG,ETOL,CHEMFL,SPECFL)
      USE PRECIS, ONLY : DBL_P
      USE GLOBAL, ONLY: SPEC,CONTROL_FL
      IMPLICIT NONE
      CHARACTER(LEN=*) :: CHEMFL,SPECFL
      INTEGER :: I,NSPMX,NTRG,NCASE,NDATA,IRDMETH,INDX_TRG(NSPMX)
      REAL(KIND=DBL_P) :: ETOL(NSPMX)
      INTEGER, PARAMETER:: IU=1
     
      WRITE(*,*) '***READ_CONTROL***'
      
      OPEN(UNIT=IU,FILE=CONTROL_FL,STATUS='OLD',FORM='FORMATTED')

      READ(IU,*) 
      READ(IU,*) CHEMFL
      READ(IU,*) 
      READ(IU,*) SPECFL
      READ(IU,*) 
      READ(IU,*) NTRG,NCASE,NDATA,IRDMETH
      READ(IU,*)
      DO I=1,NTRG
       READ(IU,*) INDX_TRG(I),ETOL(I)      
      ENDDO
      
      WRITE(*,*) 'CHEMICAL MECHANISM FILE:',TRIM(ADJUSTL(CHEMFL))
      WRITE(*,*) 'SPECIES FILE:',TRIM(ADJUSTL(SPECFL))
      WRITE(*,*) 'NO TARGETS SPECIES, NO CASES, DATA SIZE, RED. METHOD'
      WRITE(*,*) NTRG,NCASE,NDATA,IRDMETH
      WRITE(*,*) 'TARGET SPECIES INDEX, TARGET SPECIES, TOL:'
      DO I=1,NTRG
       WRITE(*,*) I,INDX_TRG(I),ETOL(I)       
      ENDDO

      CLOSE(IU)

      IF(IRDMETH.EQ.1) THEN
       WRITE(*,*) 'REDUCTION USING DRG WITH DFS'
      ELSEIF(IRDMETH.EQ.2) THEN
       WRITE(*,*) 'REDUCTION USING DRGEP WITH DIJIKSTRAS ALGORITHM'
      ELSEIF(IRDMETH.EQ.3) THEN
       WRITE(*,*) 'REDUCTION USING JAC-DRGEP WITH DIJIKSTRAS ALGORITHM'
      ENDIF

      WRITE(*,*) '***READ_CONTROL***'
                 
      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE READ_RATES(IDAT,DIR,RNM,NSPEC,NREAC,WIJ,RR,JIJ)
      USE GLOBAL, ONLY : NSMX
      USE PRECIS, ONLY : DBL_P
      IMPLICIT NONE
      INTEGER I,J,IDAT,NSPEC,NREAC
      CHARACTER(LEN=*) :: DIR,RNM
      REAL(KIND=DBL_P) :: WIJ(NREAC,NSPEC),RR(NREAC),JIJ(NSPEC,NSPEC)
      CHARACTER(LEN=4) :: CDAT
      CHARACTER(LEN=NSMX) :: FLNM
      INTEGER :: NSP,NRE
      REAL(KIND=DBL_P) :: TIME,THETA
 
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
      SUBROUTINE READ_REACTION_RATES(IDAT,DIR,RNM,NREAC,RR)
      USE GLOBAL, ONLY : NSMX
      USE PRECIS, ONLY : DBL_P
      IMPLICIT NONE
      INTEGER I,J,IDAT,NREAC,TSTEP
      CHARACTER(LEN=*) :: DIR,RNM
      REAL(KIND=DBL_P) :: RR(NREAC)
      CHARACTER(LEN=8) :: CDAT
      CHARACTER(LEN=NSMX) :: FLNM
      INTEGER :: NSP,NRE
      REAL(KIND=DBL_P) :: TIME
 
      !BUILD FILENAMES
      WRITE(CDAT,'(I8.8)') IDAT     
      FLNM=TRIM(ADJUSTL(DIR))//RNM//'_'//CDAT//'.dat'

      OPEN(UNIT=1,FILE=FLNM,STATUS='OLD',FORM='UNFORMATTED')
      READ(1) TSTEP,TIME,NRE
      CLOSE(1)
      
      IF(NRE.NE.NREAC) THEN
       WRITE(*,*) 'READ_RATES:ERROR, MISMATCH',NRE,NREAC
       WRITE(*,*) 'TERMINATING ...'
       STOP
      ELSE
       OPEN(UNIT=1,FILE=FLNM,STATUS='OLD',FORM='UNFORMATTED')
       READ(1) TSTEP,TIME,NRE,RR
       CLOSE(1)
      ENDIF
           
      END SUBROUTINE
      !-----------------------------------------------------------------
      SUBROUTINE READ_SPECIES_RATES_MATRIX(IDAT,DIR,RNM,NSPEC,NREAC,WIJ)
      USE GLOBAL, ONLY : NSMX
      USE PRECIS, ONLY : DBL_P
      IMPLICIT NONE
      INTEGER I,J,IDAT,NSPEC,NREAC
      CHARACTER(LEN=*) :: DIR,RNM
      REAL(KIND=DBL_P) :: WIJ(NREAC,NSPEC)
      CHARACTER(LEN=4) :: CDAT
      CHARACTER(LEN=NSMX) :: FLNM
      INTEGER :: NSP,NRE
      REAL(KIND=DBL_P) :: TIME
 
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
      
