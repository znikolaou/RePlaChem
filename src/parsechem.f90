      !-----------------------------------------------------------------
      SUBROUTINE READ_CHEM(DIR,CHEMFL,SPECFL)
      !      
      ! AUTHOR: Z. NIKOLAOU
      !
      USE GLOBAL, ONLY : NSPMX,NREMX
      USE ZDPLASKIN_PARSE, ONLY : NELEM,NSPEC,NREAC,ELEM,SPEC,REAC, &
                                  REAC_F,IS_SPEC_CHARGED,REAC_SPEC, &
                                  ZDP_INIT
      IMPLICIT NONE
      CHARACTER(LEN=*) :: DIR,CHEMFL,SPECFL
    
      WRITE(*,*) '***READ_CHEM***'

      CALL REMOVE_TABS_FROM_FILE(DIR//SPECFL)
      CALL REMOVE_TABS_FROM_FILE(DIR//CHEMFL)
      
      !READ CHEMISTRY- CHEM. MECH. FORMAT DEPENDENT
      CALL ZDP_INIT(DIR//CHEMFL) 
     
      !TODO:
      CALL SET_STOICH_COEFFS()
      !CALL CHECK_STOICHIOMETRY()
       
      WRITE(*,*) '***READ_CHEM***'

      RETURN
      END
      !------------------------------------------------------------------
      SUBROUTINE SET_STOICH_COEFFS()
      USE GLOBAL, ONLY : NUR,NUP
      USE ZDPLASKIN_PARSE, ONLY: NSPEC,NREAC,SPEC,REAC
      IMPLICIT NONE
      INTEGER :: I

      NUR(1:NREAC,1:NSPEC)=0.0e0
      NUP(1:NREAC,1:NSPEC)=0.0e0
      WRITE(*,*) 'REACTION STOICH. COEFFS:'
      DO I=1,NREAC
       CALL SET_REAC_STOICH_COEFFS(I,NSPEC,NREAC,REAC(I),SPEC)
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE SET_REAC_STOICH_COEFFS(IR,NSPEC,NREAC,REAC,SPEC)
      !TODO: FIX FOR ANY_NEUTRAL REACTIONS 
      USE GLOBAL, ONLY : NSMX,NSPMX,NUR,NUP,EQ_SEP,RSPEC
      IMPLICIT NONE
      INTEGER :: IR,NSPEC,NREAC,NA,NR,NP,I,J
      CHARACTER(LEN=*) :: REAC,SPEC(NSPEC)
      CHARACTER(LEN=NSMX) :: RANDP(2),RLIST(NSPMX),PLIST(NSPMX), &
                             REACF,CH
      DOUBLE PRECISION :: NUM

      WRITE(*,*) IR,TRIM(ADJUSTL(REAC))

      CALL FORMAT_REACTION(NSMX,REAC,REACF)
      CALL SPLIT_STRING(REACF,EQ_SEP,2,RANDP,NA)
      CALL SPLIT_STRING_WITH_SPACES(RANDP(1),NSPMX,RLIST,NR)
      CALL SPLIT_STRING_WITH_SPACES(RANDP(2),NSPMX,PLIST,NP)
      WRITE(*,*) TRIM(ADJUSTL(REACF))
     
      WRITE(*,*) 'REACTANTS: ',TRIM(ADJUSTL(RANDP(1))) 
      DO I=1,NR
       CALL SPLIT_NUM_AND_CHAR(RLIST(I),NSMX,NUM,CH)
       WRITE(*,*) TRIM(ADJUSTL(CH)),': ',NUM
       DO J=1,NSPEC
        IF(TRIM(ADJUSTL(CH)).EQ.TRIM(ADJUSTL(SPEC(J)))) THEN
         NUR(IR,J)=NUR(IR,J)+NUM
        ENDIF
       ENDDO
      ENDDO
      WRITE(*,*) 'PRODUCTS: ',TRIM(ADJUSTL(RANDP(2)))
      DO I=1,NP
       CALL SPLIT_NUM_AND_CHAR(PLIST(I),NSMX,NUM,CH)
       WRITE(*,*) TRIM(ADJUSTL(CH)),': ',NUM
       DO J=1,NSPEC
        IF(TRIM(ADJUSTL(CH)).EQ.TRIM(ADJUSTL(SPEC(J)))) THEN
         NUP(IR,J)=NUP(IR,J)+NUM
        ENDIF
       ENDDO
      ENDDO

      DO I=1,NSPEC
       IF(RSPEC(IR,I).EQ.1) THEN
        WRITE(*,*) TRIM(ADJUSTL(SPEC(I))),': ','NUR=',NUR(IR,I), &
                   'NUP=',NUP(IR,I)
       ENDIF
      ENDDO


      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE FORMAT_REACTION(NSMX,REAC,REACF)
      IMPLICIT NONE
      INTEGER :: NSMX
      CHARACTER(LEN=NSMX) :: REAC,REACF,CWRK,C
    
      CALL REPLACE_TEXT(REAC,'^+','^POS',CWRK,NSMX)
      C=CWRK
      CALL REPLACE_TEXT(C,'^-','^NEG',CWRK,NSMX)
      C=CWRK
      CALL REPLACE_TEXT(C,'=>',' => ',CWRK,NSMX)
      C=CWRK
      CALL REPLACE_TEXT(C,'+',' ',CWRK,NSMX)
      C=CWRK
      CALL REPLACE_TEXT(C,'NEG','-',CWRK,NSMX)
      C=CWRK
      CALL REPLACE_TEXT(C,'POS','+',CWRK,NSMX)
      C=' '//TRIM(ADJUSTL(CWRK))//' '
      REACF=C
      
      RETURN
      END
      !-----------------------------------------------------------------
      FUNCTION GET_SPECIES_INDEX(SP)
      USE GLOBAL, ONLY : NSPEC,SPEC 
      IMPLICIT NONE
      CHARACTER(LEN=*) :: SP
      INTEGER :: I,GET_SPECIES_INDEX

      GET_SPECIES_INDEX=0
      DO I=1,NSPEC
       IF(TRIM(ADJUSTL(SP)).EQ.TRIM(ADJUSTL(SPEC(I)))) THEN
         GET_SPECIES_INDEX=I
         EXIT
       ENDIF
      ENDDO
      
      END FUNCTION
      !-----------------------------------------------------------------
      SUBROUTINE REMOVE_DUPLICATE_REACTIONS(NREAC,CREAC,CREAC_F,NF)
      IMPLICIT NONE
      INTEGER :: NREAC,NF
      CHARACTER(LEN=*) :: CREAC(NREAC),CREAC_F(NREAC)

      CALL REMOVE_DUPLICATE_STRINGS(NREAC,CREAC,CREAC_F,NF)
      WRITE(*,*) 'REMOVED',NREAC-NF,'DUPLICATE REACTIONS'

      RETURN
      END  
      !-----------------------------------------------------------------
