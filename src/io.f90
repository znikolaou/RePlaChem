      !-----------------------------------------------------------------
      !
      ! AUTHOR: Z. NIKOLAOU
      !
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
      SUBROUTINE READ_REACTION_RATES(IDAT,DIR,RR)
      USE GLOBAL, ONLY : NSMX,NREAC,REAC_RATE_FL
      IMPLICIT NONE
      INTEGER I,J,IDAT,TSTEP,NSP,NRE,COUNT
      CHARACTER(LEN=*) :: DIR
      DOUBLE PRECISION :: TIME,RR(NREAC)
      CHARACTER(LEN=8) :: CDAT
      CHARACTER(LEN=NSMX) :: FLNM
 
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
      SUBROUTINE WRITE_SUMMARY(TRGD_USET,STATS,SP_USET,RE_USET)
      USE GLOBAL, ONLY : NTRG,NSPEC,NREAC,NSMX,SPEC,INDX_TRG,OUTDIR, &
                         REAC,STATS_FL,LOGO,AUTHOR,DASH,ETOL
      IMPLICIT NONE
      CHARACTER(LEN=72) :: BUILD_LINE
      CHARACTER(LEN=1) :: IS_ACCEPT
      CHARACTER(LEN=NSMX) :: SRT_SPEC(NSPEC)
      CHARACTER(LEN=30), PARAMETER :: FMTA='(2(I10X)A10XF10.6XA10)', &
                                      FMTB='(3(AX)I4XAXI4XAXF10.6)' 
      INTEGER, PARAMETER :: IU=1
      INTEGER :: I,J,IT,IC,TRGD_USET(NTRG,NSPEC),SP_USET(NSPEC), &
                 RE_USET(NREAC),NSP_RED,NRE_RED
      DOUBLE PRECISION :: STATS(NTRG,NSPEC,3),SRT_STATS(NSPEC)

      OPEN(UNIT=IU,FILE=OUTDIR//STATS_FL,STATUS='REPLACE', &
           FORM='FORMATTED')

      WRITE(IU,'(A/A/A)') LOGO,AUTHOR,DASH
      WRITE(IU,'(AXI4)') 'NO TARGETS:',NTRG

      DO I=1,NTRG
       IT=INDX_TRG(I)
       WRITE(IU,'(A)') BUILD_LINE(1)
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
      ENDDO
      
      WRITE(IU,'(A)') BUILD_LINE(1)
      
      NSP_RED=SUM(SP_USET)
      NRE_RED=SUM(RE_USET)
      IF(NSP_RED.EQ.NSPEC) THEN
       WRITE(*,*) '*WARNING: NO REDUCTION!'
      ENDIF
      
      WRITE(IU,'(A)') BUILD_LINE(1)
      
      WRITE(IU,'(AXI4)') 'ELIMINATED SPECIES:',NSPEC-NSP_RED
      IC=0
      DO I=1,NSPEC
       IF(SP_USET(I).EQ.0) THEN
        IC=IC+1
        WRITE(IU,'(I4XA)') IC,TRIM(ADJUSTL(SPEC(I)))
       ENDIF
      ENDDO
      
      WRITE(IU,'(A)') BUILD_LINE(1)

      WRITE(IU,'(AXI4)') 'KEPT SPECIES:',NSP_RED
      IC=0
      DO I=1,NSPEC
       IF(SP_USET(I).EQ.1) THEN
        IC=IC+1
        WRITE(IU,'(I4XA)') IC,TRIM(ADJUSTL(SPEC(I)))
       ENDIF
      ENDDO

      WRITE(IU,'(A)') BUILD_LINE(1)
      
      WRITE(IU,'(AXI4)') 'ELIMINATED REACTIONS:',NREAC-NRE_RED
      IC=0
      DO I=1,NREAC
       IF(RE_USET(I).EQ.0) THEN
        IC=IC+1
        WRITE(IU,'(I4XA)') IC,TRIM(ADJUSTL(REAC(I)))
       ENDIF
      ENDDO
      
      WRITE(IU,'(A)') BUILD_LINE(1)

      WRITE(IU,'(AXI4)') 'KEPT REACTIONS:',NRE_RED
      IC=0
      DO I=1,NREAC
       IF(RE_USET(I).EQ.1) THEN
        IC=IC+1
        WRITE(IU,'(I4XA)') IC,TRIM(ADJUSTL(REAC(I)))
       ENDIF
      ENDDO
      
      WRITE(IU,'(A)') BUILD_LINE(1)

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
      
