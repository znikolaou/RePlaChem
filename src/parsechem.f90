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
      USE GLOBAL, ONLY : NSMX,EQ_SEP
      USE ZDPLASKIN_PARSE, ONLY: NSPEC,NREAC,SPEC,REAC,RSPEC
      IMPLICIT NONE
      INTEGER, PARAMETER :: NC=2
      INTEGER :: I,J,NA,NUR,NUP,GET_TEXT_WITH_SPACES_COUNT
      CHARACTER(LEN=NSMX) :: COLMS(NC),REACF(NREAC)

      WRITE(*,*) 'STOICH. COEFFS:'
      CALL FORMAT_REACTIONS(NREAC,REAC(1:NREAC),REACF)
      DO I=1,NREAC
       !TODO:USE REAC RATHER THAN REAC_F
       CALL SPLIT_STRING(TRIM(ADJUSTL(REACF(I))), &
                         EQ_SEP,NC,COLMS,NA)
       WRITE(*,*) I,TRIM(ADJUSTL(REACF(I)))
       !WRITE(*,*) 'REAC:',TRIM(ADJUSTL(COLMS(1)))
       !WRITE(*,*) 'PROD:',TRIM(ADJUSTL(COLMS(2)))
       DO J=1,NSPEC
        IF(RSPEC(I,J).EQ.1) THEN
         NUR=GET_TEXT_WITH_SPACES_COUNT(COLMS(1)//' *', &
                            TRIM(ADJUSTL(SPEC(J))))
         NUP=GET_TEXT_WITH_SPACES_COUNT('* '//COLMS(2), &
                            TRIM(ADJUSTL(SPEC(J))))
         WRITE(*,*) TRIM(ADJUSTL(SPEC(J))),NUR,NUP
        ENDIF
       ENDDO
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE FORMAT_REACTIONS(NREAC,REAC,REACF)
      USE GLOBAL, ONLY: NSMX
      IMPLICIT NONE
      INTEGER :: I,NREAC
      CHARACTER(LEN=*) :: REAC(NREAC),REACF(NREAC)
      CHARACTER(LEN=NSMX) :: CWRK,C
    
      DO I=1,NREAC
       CALL REPLACE_TEXT(REAC(I),'^+','^POS',CWRK,NSMX)
       C=CWRK
       CALL REPLACE_TEXT(C,'^-','^NEG',CWRK,NSMX)
       C=CWRK
       CALL REPLACE_TEXT(C,'=>',' => ',CWRK,NSMX)
       C=CWRK
       CALL REPLACE_TEXT(C,'+',' + ',CWRK,NSMX)
       C=CWRK
       CALL REPLACE_TEXT(C,'NEG','-',CWRK,NSMX)
       C=CWRK
       CALL REPLACE_TEXT(C,'POS','+',CWRK,NSMX)
       C='* '//TRIM(ADJUSTL(CWRK))//' *'
       REACF(I)=C
      ENDDO

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
