      !-----------------------------------------------------------------
      SUBROUTINE READ_CHEM(DIR,CHEMFL,SPECFL)
      !      
      ! AUTHOR: Z. NIKOLAOU
      !
      USE GLOBAL, ONLY : NSPMX,NREMX
      USE ZDPLASKIN_PARSE, ONLY : NELEM,NSPEC,NREAC,ELEM,SPEC,REAC, &
                                  IS_SPEC_CHARGED,REAC_SPEC,ZDP_INIT
      IMPLICIT NONE
      CHARACTER(LEN=*) :: DIR,CHEMFL,SPECFL
    
      WRITE(*,*) '***READ_CHEM***'

      CALL REMOVE_TABS_FROM_FILE(DIR//SPECFL)
      CALL REMOVE_TABS_FROM_FILE(DIR//CHEMFL)
      
      !READ CHEMISTRY- CHEM. MECH. FORMAT DEPENDENT
      CALL ZDP_INIT(DIR//CHEMFL) 
     
      !TODO:
      !CALL CHECK_STOICHIOMETRY()
       
      WRITE(*,*) '***READ_CHEM***'

      RETURN
      END
      !------------------------------------------------------------------
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
