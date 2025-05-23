      !-----------------------------------------------------------------
      !
      ! AUTHOR: Z. NIKOLAOU
      !
      !-----------------------------------------------------------------
      SUBROUTINE WRITE_CHEM_MECH_FMT_ZDP(SETE,SETSP,SETSP_BOLS,SETRE)
      USE GLOBAL, ONLY : OUTDIR,CHEMRED_FL,NTRG,NELEM,NSPEC,NREAC, &
       NREAC_DOLLAR,INDX_TRG,ELEM,SPEC,REAC,REAC_SEC_DOLLAR_LIST, &
       REAC_CONST,ETOL
      IMPLICIT NONE
      INTEGER :: NBOLS,I
      INTEGER, PARAMETER :: IO=7
      INTEGER :: SETE(NELEM),SETSP(NSPEC),SETSP_BOLS(NSPEC), &
                 SETRE(NREAC)
      !CHARACTER(LEN=*) :: ELEM(NELEM),SPEC(NSPEC), &
      !                    REAC_DOLLAR_LIST(NREAC_DOLLAR), &
      !                    REAC(NREAC),REAC_CONST(NREAC)
      CHARACTER(LEN=LEN(OUTDIR)+LEN(CHEMRED_FL)+1) :: FL
      CHARACTER(LEN=8)  :: DATE
      CHARACTER(LEN=10) :: TIME
      CHARACTER(LEN=5)  :: ZONE

      CALL DATETIME(DATE,TIME,ZONE)
    
      FL=OUTDIR//'/'//CHEMRED_FL
      OPEN(UNIT=IO,FILE=FL,STATUS='REPLACE',FORM='FORMATTED')

      WRITE(IO,'(A)') '# PLASEREDCHEM.V01'
      WRITE(IO,'(A)') '#'
      WRITE(IO,'(A)') '# AUTHOR: ZACHARIAS NIKOLAOU'
      WRITE(IO,'(A)') '#'
      WRITE(IO,'(A)') '# THIS IS AN AUTO-GENERATED FILE.'
      WRITE(IO,'(A)') '#'
      WRITE(IO,'(4(AX))') '# CREATED: ',TRIM(DATE),TRIM(TIME),TRIM(ZONE)
      WRITE(IO,'(A)') '#'
      WRITE(IO,'(A)') '# REDUCTION PARAMS (TARGET SPECIES AND ETOL):'
      DO I=1,NTRG
       WRITE(IO,'(AXAXE12.5)') '# '//TRIM(ADJUSTL(SPEC(INDX_TRG(I)))), &
                              'ETOL=',ETOL(I)
      ENDDO
      WRITE(IO,'(A)') '#'
      WRITE(IO,'(A)') '# DETAILED MECHANISM:'
      WRITE(IO,'(AXI4)') '# NREAC=',NREAC
      WRITE(IO,'(AXI4)') '# NSPEC=',NSPEC
      WRITE(IO,'(A)') '# REDUCED MECHANISM:'
      WRITE(IO,'(AXI4)') '# NREAC=',SUM(SETRE)
      WRITE(IO,'(AXI4)') '# NSPEC=',SUM(SETSP)
      WRITE(IO,'(A)') '#' 
     
      CALL WRITE_SECTION_FMT_ZDP(IO,'ELEMENTS',NELEM,ELEM,SETE)
      CALL WRITE_SECTION_FMT_ZDP(IO,'SPECIES',NSPEC,SPEC,SETSP)
      CALL WRITE_SECTION_FMT_ZDP(IO,'BOLSIG',NSPEC,SPEC,SETSP_BOLS)
      CALL WRITE_SECTION_REAC_FMT_ZDP(IO,'REACTIONS',NREAC, &
       NREAC_DOLLAR,REAC,REAC_CONST,SETRE,REAC_SEC_DOLLAR_LIST)
 
      CLOSE(IO)

      RETURN
      END 
      !-----------------------------------------------------------------
      SUBROUTINE WRITE_REAC_CONST(IU,N,DOLLARLIST)
      IMPLICIT NONE
      INTEGER :: IU,N
      CHARACTER(LEN=*) :: DOLLARLIST(N)

      WRITE(IU,'(A)') '$'

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

