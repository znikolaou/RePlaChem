      !-----------------------------------------------------------------
      !
      ! AUTHOR: Z. NIKOLAOU
      !
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
      SUBROUTINE READ_REACTION_RATES(IDAT,DIR,RR)
      USE GLOBAL, ONLY : NSMX,NREAC,REAC_RATE_FL
      IMPLICIT NONE
      INTEGER I,J,IDAT,TSTEP
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
      SUBROUTINE WRITE_CHEM_MECH_FMT_ZDP(SETE,SETSP,SETSP_BOLS,SETRE)
      USE GLOBAL
      IMPLICIT NONE
      INTEGER :: I
      INTEGER, PARAMETER :: IO=7
      INTEGER :: SETE(NELEM),SETSP(NSPEC),SETSP_BOLS(NSPEC),SETRE(NREAC)
      CHARACTER(LEN=LEN(OUTDIR)+LEN(CHEMRED_FL)+1) :: FL
      CHARACTER(LEN=8)  :: DATE
      CHARACTER(LEN=10) :: TIME
      CHARACTER(LEN=5)  :: ZONE

      CALL DATETIME(DATE,TIME,ZONE)
    
      FL=OUTDIR//'/'//CHEMRED_FL
      OPEN(UNIT=IO,FILE=FL,STATUS='REPLACE',FORM='FORMATTED')

      WRITE(IO,'(A)') LOGO
      WRITE(IO,'(A)') AUTHOR
      WRITE(IO,'(A)') DASH
      WRITE(IO,'(A)') '#'
      WRITE(IO,'(A)') '# THIS IS AN AUTO-GENERATED FILE.'
      WRITE(IO,'(4(AX))') '# CREATED: ',TRIM(DATE),TRIM(TIME),TRIM(ZONE)
      WRITE(IO,'(A)') '#'
      WRITE(IO,'(A1X2(A10X))') '#','TARGETS', 'ACC. THR.'
      DO I=1,NTRG
       WRITE(IO,'(A1XA10XF10.6)') '#', &
               TRIM(ADJUSTL(SPEC(INDX_TRG(I)))),ETOL(I)
      ENDDO
      WRITE(IO,'(A)') '#'
      WRITE(IO,'(A1(8X)AXA)') '#','DET. MECH.','RED. MECH.'
      WRITE(IO,'(A1XA6X2(I8X))') '#','NREAC=',NREAC,SUM(SETRE)
      WRITE(IO,'(A1XA6X2(I8X))') '#','NSPEC=',NSPEC,SUM(SETSP)
     
      CALL WRITE_SECTION_FMT_ZDP(IO,'ELEMENTS',NELEM,ELEM,SETE)
      CALL WRITE_SECTION_FMT_ZDP(IO,'SPECIES',NSPEC,SPEC,SETSP)
      CALL WRITE_SECTION_FMT_ZDP(IO,'BOLSIG',NSPEC,SPEC,SETSP_BOLS)
      CALL WRITE_SECTION_REAC_FMT_ZDP(IO,'REACTIONS',NREAC, &
       NREAC_DOLLAR,REAC,REAC_CONST,SETRE,REAC_SEC_DOLLAR_LIST)
 
      CLOSE(IO)

      RETURN
      END 
      !-----------------------------------------------------------------
      SUBROUTINE WRITE_SECTION_FMT_ZDP(IU,SECNAME,N,CLIST,SET)
      IMPLICIT NONE
      INTEGER :: IU,N,I,SET(N)
      CHARACTER(LEN=*) :: SECNAME,CLIST(N)
 
      WRITE(IU,'(A)') '#'
      WRITE(IU,'(A)') SECNAME
      DO I=1,N
       IF(SET(I).EQ.1) THEN !KEEP
        WRITE(IU,'(A)') TRIM(ADJUSTL(CLIST(I)))
       ENDIF
      ENDDO
      WRITE(IU,'(A)') 'END'
      WRITE(IU,'(A)') '#' 
      RETURN

      END
      !-----------------------------------------------------------------
      SUBROUTINE WRITE_SECTION_REAC_FMT_ZDP(IU,SECNAME,N,NDOLLAR,CLIST, &
                                            CONST,SET,REAC_DOLLAR_LIST)
      IMPLICIT NONE
      INTEGER :: IU,N,I,NDOLLAR,SET(N)
      CHARACTER(LEN=*) :: SECNAME,CLIST(N),CONST(N), &
                          REAC_DOLLAR_LIST(NDOLLAR)
 
      WRITE(IU,'(A)') '#'
      WRITE(IU,'(A)') SECNAME
      WRITE(IU,'(A)') '# USER-PROVIDED EXPRESSIONS:'
      DO I=1,NDOLLAR
       WRITE(IU,'(A)') TRIM(ADJUSTL(REAC_DOLLAR_LIST(I)))
      ENDDO
      WRITE(IU,'(A)') '#'
      DO I=1,N
       IF(SET(I).EQ.1) THEN !KEEP
       WRITE(IU,'(A)') TRIM(ADJUSTL(CLIST(I)))//' !'// & 
                       TRIM(ADJUSTL(CONST(I))) 
       ENDIF
      ENDDO
      WRITE(IU,'(A)') 'END'
      WRITE(IU,'(A)') '#'
 
      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE WRITE_OICS(ICASE,IDATA,NTRG,INDX_TRG,OIC)
      USE GLOBAL, ONLY: NSFLMX,NSPEC,SPEC,OUTDIR
      IMPLICIT NONE
      INTEGER, PARAMETER :: IU=1
      CHARACTER(LEN=NSFLMX) :: FL
      CHARACTER(LEN=3) :: IC
      CHARACTER(LEN=8) :: ID
      INTEGER :: ICASE,IDATA,NTRG,INDX_TRG(NTRG)
      DOUBLE PRECISION :: OIC(NTRG,NSPEC)

      WRITE(IC,'(I3.3)') ICASE
      WRITE(ID,'(I8.8)') IDATA
      FL=TRIM(ADJUSTL(OUTDIR))//'case'//IC//'data'//ID//'.dat'
      OPEN(UNIT=IU,FILE=FL,STATUS='REPLACE',FORM='UNFORMATTED')
      WRITE(IU) NTRG,NSPEC,OIC
      CLOSE(IU)
      
      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE WRITE_STATS(TRGD_USET,STATS)
      USE GLOBAL, ONLY : NTRG,NSPEC,NSMX,SPEC,INDX_TRG,OUTDIR, &
                         STATS_FL,LOGO,AUTHOR,DASH,ETOL
      IMPLICIT NONE
      CHARACTER(LEN=1) :: IS_ACCEPT
      CHARACTER(LEN=NSMX) :: SRT_SPEC(NSPEC)
      CHARACTER(LEN=30), PARAMETER :: FMTA='(2(I10X)A10XF10.6XA10)', &
                                      FMTB='(3(AX)I4XAXI4XAXF10.6)' 
      INTEGER, PARAMETER :: IU=1
      INTEGER :: I,J,IT,IC,TRGD_USET(NTRG,NSPEC)
      DOUBLE PRECISION :: STATS(NTRG,NSPEC,3),SRT_STATS(NSPEC)

      OPEN(UNIT=IU,FILE=OUTDIR//STATS_FL,STATUS='REPLACE', &
           FORM='FORMATTED')

      WRITE(IU,'(A/A/A)') LOGO,AUTHOR,DASH
      WRITE(IU,*)
      WRITE(IU,'(AXI4)') 'NO TARGETS:',NTRG
      WRITE(IU,*)

      DO I=1,NTRG
       IT=INDX_TRG(I)
       WRITE(IU,'(4(A10X))') 'TARGET','INDEX','NO DEP.','ACC. THR.'
       WRITE(IU,'(A10XI10XI10XF10.6)') TRIM(ADJUSTL(SPEC(IT))),IT, &
                  SUM(TRGD_USET(I,:))-1,ETOL(I)    
       WRITE(IU,'(6(A10X))') 'COUNT', 'INDEX','SPEC','MAX OIC', &
                             'ACCEPT(A)'
       CALL SORT(NSPEC,NSMX,STATS(I,:,3),SPEC(1:NSPEC),SRT_STATS, &
                 SRT_SPEC)
       IC=0
       DO J=1,NSPEC
        IF(J.NE.IT) THEN
         IC=IC+1
         WRITE(IU,FMTA) IC,J,TRIM(ADJUSTL(SRT_SPEC(J))),SRT_STATS(J), &
                 IS_ACCEPT(SRT_STATS(J),ETOL(I))
        ENDIF 
       ENDDO
       
      WRITE(IU,*) 
      ENDDO

      WRITE(IU,'(A)') '!END'
   
      CLOSE(IU)

      RETURN
      END
      !-----------------------------------------------------------------
      FUNCTION IS_ACCEPT(V,E)
      IMPLICIT NONE
      DOUBLE PRECISION :: V,E
      CHARACTER(LEN=1) :: IS_ACCEPT

      IF(V.GE.E) THEN
       IS_ACCEPT='A'
      ELSE
       IS_ACCEPT='X'
      ENDIF

      END FUNCTION
      !-----------------------------------------------------------------
      SUBROUTINE WRITE_REDUCTION_INFO(SP_USET,RE_USET)
      USE GLOBAL, ONLY: NSPEC,NREAC,SPEC,REAC
      IMPLICIT NONE
      INTEGER :: I,NSP_RED,NRE_RED,SP_USET(NSPEC),RE_USET(NREAC)

      WRITE(*,*) 'ELIMINATED SPECIES:'
      DO I=1,NSPEC
       IF(SP_USET(I).EQ.0) THEN
        WRITE(*,*) TRIM(ADJUSTL(SPEC(I)))
       ENDIF
      ENDDO

      WRITE(*,*) 'KEPT SPECIES:'
      DO I=1,NSPEC
       IF(SP_USET(I).EQ.1) THEN
        WRITE(*,*) TRIM(ADJUSTL(SPEC(I)))
       ENDIF
      ENDDO

      WRITE(*,*) 'KEPT REACTIONS:'
      DO I=1,NREAC
       IF(RE_USET(I).EQ.1) THEN
        WRITE(*,*) TRIM(ADJUSTL(REAC(I)))
       ENDIF
      ENDDO

      NSP_RED=SUM(SP_USET)
      NRE_RED=SUM(RE_USET)
      IF(NSP_RED.LE.0.OR.NRE_RED.LE.0) THEN
       WRITE(*,*) 'ERROR: NO RED. MECH. CREATED TERMINATING ...'
       STOP
      ENDIF
      IF(NSP_RED.EQ.NSPEC) THEN
       WRITE(*,*) 'WARNING: NO REDUCTION!'
      ENDIF

      RETURN
      END
      !-----------------------------------------------------------------
      FUNCTION GET_FORMAT(ID)
      IMPLICIT NONE
      INTEGER :: ID
      CHARACTER(LEN=50) :: GET_FORMAT

      GET_FORMAT='*'
      IF(ID.EQ.1) GET_FORMAT='(A10)'
     

      END FUNCTION
      !-----------------------------------------------------------------
