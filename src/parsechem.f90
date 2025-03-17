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
      CALL ZDP_INIT(DIR//CHEMFL) !NELEM,NSPEC,NREAC,ELEM,SPEC,REAC, &
                                 !IS_SPEC_CHARGED,REAC_SPEC
     
      !TODO:
      !CALL CHECK_STOICHIOMETRY()
       
      WRITE(*,*) '***READ_CHEM***'

      RETURN
      END
      !------------------------------------------------------------------
      !TODO: REFORMAT
      SUBROUTINE SEPARATE_REACTIONS(N,CREAC,REAC,PROD)
      USE GLOBAL, ONLY : NSMX
      IMPLICIT NONE
      INTEGER :: N,I
      CHARACTER(LEN=*) :: CREAC(N),REAC(N),PROD(N)
      CHARACTER(LEN=NSMX) :: CLIST(2)

      DO I=1,N
       IF(INDEX(CREAC(I),'bolsig').EQ.0) THEN 
        CALL SEPARATE_STRING(CREAC(I),'=>',CLIST)
       ELSE
        CALL SEPARATE_STRING(CREAC(I),'->',CLIST)
       ENDIF
        REAC(I)=CLIST(1)
        PROD(I)=CLIST(2)
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
      !TODO: CHECK WITH TEST
      SUBROUTINE GET_REAC_SPEC_FROM_REAC(NSPEC,NREAC,CREAC,CSPEC, &
                      CREAC_SPEC,NC)
      USE GLOBAL, ONLY : NSMX        
      IMPLICIT NONE
      !LOGICAL :: IS_ANY_NEUTRAL_REACTION,IS_CHARGED_SPECIES, &
      !           IS_STRING_PRESENT
      INTEGER :: NSPEC,NREAC,I,J,K,NC(NREAC),NA
      CHARACTER(LEN=*) :: CREAC(NREAC),CSPEC(NSPEC), &
                          CREAC_SPEC(NREAC,NSPEC)
      CHARACTER(LEN=NSMX) :: C(1),CF(1)
      
      DO I=1,NREAC
       C=CREAC(I)
       !CALL SET_FORMATTED_REACTIONS()
       CALL SPLIT_STRING(CF,' ',NSPEC,CREAC_SPEC(I,:),NA)
       NC(I)=NA
       !IF(IS_ANY_NEUTRAL_REACTION(CREAC(I))) THEN
       ! DO J=1,NSPEC
       !  IF(.NOT.IS_CHARGED_SPECIES(CSPEC(J))) THEN
       !   IF(.NOT.IS_STRING_PRESENT(NSPEC, &
       !           CREAC_SPEC(I,:),CSPEC(J))) THEN
       !    NC(I)=NC(I)+1
       !    CREAC_SPEC(I,NC(I))=CSPEC(J)
       !   ENDIF
       !  ENDIF
       ! ENDDO
       !ENDIF
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE SET_REACTION_SPECIES()
      USE GLOBAL
      IMPLICIT NONE
    
      !TODO: REFORMAT
      !CALL SET_SPECIES_FROM_LIST()
      !CALL SET_E_FOR_BOLSIG_IF_ANY()
      !CALL SET_NEUTRAL_IF_ANY()

      RETURN
      END
      !-----------------------------------------------------------------
      FUNCTION GET_SPECIES_INDEX(SP)
      USE GLOBAL 
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
      LOGICAL FUNCTION IS_LE_TO(N,NMAX)
      IMPLICIT NONE
      INTEGER :: N,NMAX
      
      IS_LE_TO=N.LE.NMAX

      END FUNCTION
      !------------------------------------------------------------------
      SUBROUTINE REMOVE_DUPLICATE_REACTIONS(NREAC,CREAC,CREAC_F,NF)
      IMPLICIT NONE
      INTEGER :: NREAC,NF
      CHARACTER(LEN=*) :: CREAC(NREAC),CREAC_F(NREAC)

      CALL REMOVE_DUPLICATE_STRINGS(NREAC,CREAC,CREAC_F,NF)
      WRITE(*,*) 'REMOVED',NREAC-NF,'DUPLICATE REACTIONS'

      RETURN
      END  
      !-----------------------------------------------------------------
      SUBROUTINE GET_CHARGED_SPECIES(NSPEC,CSPEC,CCHAR,NC)
      IMPLICIT NONE
      !LOGICAL :: IS_CHARGED_SPECIES
      INTEGER :: NSPEC,NC,I,J
      CHARACTER(LEN=*) :: CSPEC(NSPEC),CCHAR(NSPEC)

      J=0
      DO I=1,NSPEC
       !IF(IS_CHARGED_SPECIES(CSPEC(I))) THEN
       ! J=J+1
       ! CCHAR(J)=CSPEC(I)
       !ENDIF
      ENDDO
      NC=J

      RETURN
      END 
      !-----------------------------------------------------------------
      
