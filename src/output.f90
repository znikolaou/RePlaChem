      !-----------------------------------------------------------------
      !
      ! AUTHOR: Z. NIKOLAOU
      !
      !-----------------------------------------------------------------
      SUBROUTINE OUTPUT(NSPEC,NREAC,SETSP,SETRE,CSPNM,CSPNMF, &
                        CRENMF,FLSP_KPP,FLRE_KPP)
      USE GLOBAL, ONLY :NSMX
      IMPLICIT NONE 
      INTEGER :: I,J,K,INDX,NSPEC,NREAC,LS
      INTEGER :: SETSP(NSPEC),SETRE(NREAC)
      CHARACTER(LEN=*) :: CSPNM(NSPEC),CSPNMF(NSPEC),CRENMF(NREAC), &
                          FLSP_KPP,FLRE_KPP
      CHARACTER(LEN=NSMX) :: C,CF
      INTEGER :: GET_TEXT_WITH_SPACES_COUNT,NSPEC_SKEL,NREAC_SKEL
      CHARACTER(LEN=8)  :: DATE
      CHARACTER(LEN=10) :: TIME
      CHARACTER(LEN=5)  :: ZONE

      CALL DATETIME(DATE,TIME,ZONE)

      NSPEC_SKEL=SUM(SETSP)  
      NREAC_SKEL=SUM(SETRE)  

      !REACTION FILE
      OPEN(UNIT=1,FILE=FLRE_KPP,STATUS='NEW',FORM='FORMATTED')

      WRITE(1,'(A)') '// AUTHOR: Z.M. NIKOLAOU'
      WRITE(1,'(A)') '// SKELETAL MECHANISM-OUTPUT FROM DRG.'
      WRITE(1,'(4(AX))') '// CREATED: ',TRIM(DATE),TRIM(TIME),TRIM(ZONE)
      WRITE(1,'(A)') '//' 
      WRITE(1,'(A)') '#EQUATIONS' 

      J=0
      DO I=1,NREAC
       IF(SETRE(I).EQ.1) THEN
        J=J+1
        C=CRENMF(I)
        !WRITE(*,*) C
        !C(1:INDEX(C,'}'))=' '
        C=TRIM(ADJUSTL(C))
        WRITE(1,'(A,I4.4,A,X,A)') '{',J,'}',TRIM(ADJUSTL(C))
       ENDIF 
      ENDDO
      !
      GOTO 1

      !OUTPUT DUMMY REACTIONS FOR ELIMINATED SPECIES
      WRITE(1,'(A)') '{DUMMY REACTIONS-ELIMINATED SPEC}'
      DO I=1,NSPEC
       IF(SETSP(I).EQ.0) THEN
        C=CSPNM(I)
        WRITE(1,'(2(A,X,A,X))') TRIM(ADJUSTL(C)),'=',TRIM(ADJUSTL(C)), &
                           ': 0.0D0 ;' 
       ENDIF
      ENDDO
      !
      WRITE(1,'(A)') '//' 
      CLOSE(1)

1     CONTINUE

      !SPECIES FILE
      OPEN(UNIT=1,FILE=FLSP_KPP,STATUS='NEW',FORM='FORMATTED')

      WRITE(1,'(A)') '// AUTHOR: Z.M. NIKOLAOU'
      WRITE(1,'(A)') '// SKELETAL MECHANISM-OUTPUT FROM DRG.'
      WRITE(1,'(4(AX))') '// CREATED: ',TRIM(DATE),TRIM(TIME),TRIM(ZONE)
      WRITE(1,'(A)') '//' 
      WRITE(1,'(A)') '#DEFVAR' 
      DO J=1,NSPEC
       C=CSPNM(J)
       CALL REMOVE_SPACES(LEN(C),C,LS)
       
       IF(SETSP(J).EQ.1) THEN
        !WRITE(*,'(A,X,I3)') C
        K=0
        I=0         
        DO WHILE(K.EQ.0.AND.I.LE.NSPEC-1)
         I=I+1
         CF=CSPNMF(I)
         CF=CF(1:INDEX(CF,'='))
         K=GET_TEXT_WITH_SPACES_COUNT(CF,C)         
         !K=GET_TEXT_WITH_SPACES_COUNT(CSPNMF(I),C)
         !write(*,*) K
        ENDDO            
        !WRITE(*,*) CSPNMF(I)
        WRITE(1,'(A)') TRIM(ADJUSTL(CSPNMF(I)))
       ENDIF
      ENDDO
      WRITE(1,'(A)') '//' 
      CLOSE(1)

      END
      !-----------------------------------------------------------------
      SUBROUTINE WRITE_CHEM_MECH_FMT_ZDP(DIR,FLSKEL,NELEM,NSPEC, &
                                         NREAC,ELEM,SPEC,REAC, &
                                         REAC_CONST, &
                                         SETE,SETSP,SETSP_BOLS,SETRE)
      IMPLICIT NONE
      INTEGER :: NELEM,NSPEC,NBOLS,NREAC
      INTEGER, PARAMETER :: IO=7
      INTEGER :: SETE(NELEM),SETSP(NSPEC),SETSP_BOLS(NSPEC),SETRE(NREAC)
      CHARACTER(LEN=*) :: ELEM(NELEM),SPEC(NSPEC), &
                          REAC(NREAC),REAC_CONST(NREAC),DIR,FLSKEL 
      CHARACTER(LEN=LEN(DIR)+LEN(FLSKEL)+1) :: FL
      CHARACTER(LEN=8)  :: DATE
      CHARACTER(LEN=10) :: TIME
      CHARACTER(LEN=5)  :: ZONE

      CALL DATETIME(DATE,TIME,ZONE)
    
      FL=DIR//'/'//FLSKEL
      OPEN(UNIT=IO,FILE=FL,STATUS='REPLACE',FORM='FORMATTED')

      WRITE(IO,'(A)') '# PLASEREDCHEM.V01'
      WRITE(IO,'(A)') '#'
      WRITE(IO,'(A)') '# AUTHOR: Z. NIKOLAOU'
      WRITE(IO,'(A)') '#'
      WRITE(IO,'(A)') '# AUTO-GENERATED FILE: SKELETAL CHEM. MECH.'
      WRITE(IO,'(A)') '# FORMAT=ZDPlasKin'
      WRITE(IO,'(4(AX))') '# CREATED: ',TRIM(DATE),TRIM(TIME),TRIM(ZONE)
      WRITE(IO,'(A)') '#' 
     
      CALL WRITE_SECTION_FMT_ZDP(IO,'ELEMENTS',NELEM,ELEM,SETE)
      CALL WRITE_SECTION_FMT_ZDP(IO,'SPECIES',NSPEC,SPEC,SETSP)
      CALL WRITE_SECTION_FMT_ZDP(IO,'BOLSIG',NSPEC,SPEC,SETSP_BOLS)
      !CALL WRITE_SECTION_FMT_ZDP(IO,'BOLSIG',NBOLS,BOLS)
      CALL WRITE_SECTION_REAC_FMT_ZDP(IO,'REACTIONS',NREAC,REAC, &
                                      REAC_CONST,SETRE)
 
      CLOSE(IO)

      RETURN
      END 
      !-----------------------------------------------------------------
      SUBROUTINE WRITE_SECTION_FMT_ZDP(IU,SECNAME,N,CLIST,SET)
      IMPLICIT NONE
      INTEGER :: IU,N,I,SET(N)
      CHARACTER(LEN=*) :: SECNAME,CLIST(N)

      WRITE(IU,'(A)') ' ' 
      WRITE(IU,'(A)') '#'
      WRITE(IU,'(A)') SECNAME
      DO I=1,N
       IF(SET(I).EQ.1) THEN !KEEP
        WRITE(IU,'(A)') TRIM(ADJUSTL(CLIST(I)))
       ENDIF
      ENDDO
      WRITE(IU,'(A)') 'END'
      WRITE(IU,'(A)') '#'
      WRITE(IU,'(A)') ' '
 
      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE WRITE_SECTION_REAC_FMT_ZDP(IU,SECNAME,N,CLIST,CONST, &
                                            SET)
      IMPLICIT NONE
      INTEGER :: IU,N,I,SET(N)
      CHARACTER(LEN=*) :: SECNAME,CLIST(N),CONST(N)

      !TODO: WRITE ALSO EXPRESSIONS WITH '$' SIGN WHICH HOLD RATE EXPRESSIONS
      WRITE(IU,'(A)') ' '
      WRITE(IU,'(A)') '#'
      WRITE(IU,'(A)') SECNAME
      DO I=1,N
       IF(SET(I).EQ.1) THEN !KEEP
       WRITE(IU,'(A)') TRIM(ADJUSTL(CLIST(I)))//' !'// & 
                       TRIM(ADJUSTL(CONST(I))) 
       ENDIF
      ENDDO
      WRITE(IU,'(A)') 'END'
      WRITE(IU,'(A)') '#'
      WRITE(IU,'(A)') ' '
 
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
      SUBROUTINE WRITE_SPECIES(NCASE,NDATA,NTRG,INDX_TRG)
      USE GLOBAL, ONLY : NSFLMX,NSPEC,SPEC,OUTDIR
      IMPLICIT NONE
      CHARACTER(LEN=NSFLMX) :: FL
      INTEGER :: NTRG,INDX_TRG(NTRG),I,NCASE,NDATA
      INTEGER, PARAMETER :: IU=1

      FL=OUTDIR//'speciesAndTargets.txt'
      OPEN(UNIT=IU,FILE=FL,STATUS='REPLACE',FORM='FORMATTED')
      WRITE(IU,*) 'No of cases:'
      WRITE(IU,*) NCASE
      WRITE(IU,*) 'No of data for each case:'
      WRITE(IU,*) NDATA
      WRITE(IU,*) 'No species:'
      WRITE(IU,*) NSPEC
      WRITE(IU,*) 'Species:'
      DO I=1,NSPEC
       WRITE(IU,*) TRIM(ADJUSTL(SPEC(I)))
      ENDDO
      WRITE(IU,*) 'No of targets:'
      WRITE(IU,*) NTRG
      WRITE(IU,*) 'Target species:'
      DO I=1,NTRG
       WRITE(IU,*) TRIM(ADJUSTL(SPEC(INDX_TRG(I))))
      ENDDO
      CLOSE(IU)

      RETURN
      END
      !-----------------------------------------------------------------
