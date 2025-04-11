      !-----------------------------------------------------------------
      SUBROUTINE READ_CHEM(DIR,CHEMFL,SPECFL)
      !      
      ! AUTHOR: Z. NIKOLAOU
      !
      USE GLOBAL
      IMPLICIT NONE
      CHARACTER(LEN=*) :: DIR,CHEMFL,SPECFL
      INTEGER I,J
      INTEGER, PARAMETER :: I1=1,I2=2
      LOGICAL :: IS_LE_TO

      CALL INIT_GLOBALS()

      CALL REMOVE_TABS_FROM_FILE(DIR//SPECFL)
      CALL REMOVE_TABS_FROM_FILE(DIR//CHEMFL)
      
      CALL SET_SPECIES(DIR//SPECFL)
      CALL SET_IS_NEUTRAL()
      CALL SET_REACTIONS(DIR//CHEMFL)
      CALL SET_FORMATTED_REACTIONS()
      CALL SET_REACTION_SPECIES()
      
      WRITE(*,*) '***READ_CHEM***'

      IF(.NOT.IS_LE_TO(NREAC,NREMX)) THEN
       WRITE(*,*) '*ERROR: NREAC>NREMX:',NREAC,NREMX
       STOP
      ENDIF
      IF(.NOT.IS_LE_TO(NSPEC,NSPMX)) THEN
       WRITE(*,*) '*ERROR: NSPEC>NSPMX:',NSPEC,NSPMX
       STOP
      ENDIF
      
      WRITE(*,*) 'SPECIES:',NSPEC 
      DO I=1,NSPEC
       WRITE(*,'(I5,X,A)') I,TRIM(ADJUSTL(SPEC(I)))
      ENDDO 
      WRITE(*,*) 'REACTIONS:',NREAC 
      DO I=1,NREAC
       WRITE(*,'(I5,X,A)') I,TRIM(ADJUSTL(REAC(I)))
      ENDDO 
     
      WRITE(*,*) 'REACTION SPECIES:'
      DO I=1,NREAC
       WRITE(*,*) I,REAC(I)
       DO J=1,NSPEC
        IF(IS_SPEC_IN_REAC(J,I)) THEN
         WRITE(*,*) J,SPEC(J)
        ENDIF
       ENDDO
      ENDDO
     
      WRITE(*,*) '***READ_CHEM***'

      RETURN
      END
      !------------------------------------------------------------------
      SUBROUTINE INIT_GLOBALS()
      USE GLOBAL
      IMPLICIT NONE
      
      REAC(1:NREMX)=''
      REACF(1:NREMX)=''
      SPEC(1:NSPMX)=''
      IS_SPEC_IN_REAC(1:NSPMX,1:NREMX)=.FALSE.
      IS_SPEC_NEUTRAL(1:NSPMX)=.TRUE.

      RETURN
      END
      !------------------------------------------------------------------
      SUBROUTINE SET_REACTIONS(FL)
      USE GLOBAL, ONLY : REAC,NREAC,NREMX
      IMPLICIT NONE
      CHARACTER(LEN=*) :: FL

      CALL READ_INT_AND_STR(FL,NREMX,REAC,NREAC)

      RETURN
      END
      !------------------------------------------------------------------
      SUBROUTINE SET_SPECIES(FL)
      USE GLOBAL, ONLY : SPEC,NSPEC,NSPMX
      IMPLICIT NONE
      CHARACTER(LEN=*) :: FL

      CALL READ_INT_AND_STR(FL,NSPMX,SPEC,NSPEC)

      RETURN
      END
      !----------------------------------------------------------------- 
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
      LOGICAL :: IS_ANY_NEUTRAL_REACTION,IS_CHARGED_SPECIES, &
                 IS_STRING_PRESENT
      INTEGER :: NSPEC,NREAC,I,J,K,NC(NREAC),NA
      CHARACTER(LEN=*) :: CREAC(NREAC),CSPEC(NSPEC), &
                          CREAC_SPEC(NREAC,NSPEC)
      CHARACTER(LEN=NSMX) :: C(1),CF(1)
      
      DO I=1,NREAC
       C=CREAC(I)
       CALL SET_FORMATTED_REACTIONS()
       CALL SPLIT_STRING(CF,' ',NSPEC,CREAC_SPEC(I,:),NA)
       NC(I)=NA
       IF(IS_ANY_NEUTRAL_REACTION(CREAC(I))) THEN
        DO J=1,NSPEC
         IF(.NOT.IS_CHARGED_SPECIES(CSPEC(J))) THEN
          IF(.NOT.IS_STRING_PRESENT(NSPEC, &
                  CREAC_SPEC(I,:),CSPEC(J))) THEN
           NC(I)=NC(I)+1
           CREAC_SPEC(I,NC(I))=CSPEC(J)
          ENDIF
         ENDIF
        ENDDO
       ENDIF
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE SET_REACTION_SPECIES()
      USE GLOBAL
      IMPLICIT NONE
    
      CALL SET_SPECIES_FROM_LIST()
      CALL SET_E_FOR_BOLSIG_IF_ANY()
      CALL SET_NEUTRAL_IF_ANY()

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
      SUBROUTINE SET_SPECIES_FROM_LIST()
      USE GLOBAL
      IMPLICIT NONE 
      INTEGER :: I,J
      LOGICAL :: IS_SPEC_IN_REACTION

      IS_SPEC_IN_REAC(1:NSPMX,1:NREMX)=.FALSE.
      DO J=1,NREAC
       DO I=1,NSPEC
        IF(IS_SPEC_IN_REACTION(REACF(J),SPEC(I))) THEN 
         IS_SPEC_IN_REAC(I,J)=.TRUE.
        ENDIF
       ENDDO
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE SET_E_FOR_BOLSIG_IF_ANY()
      USE GLOBAL 
      IMPLICIT NONE
      INTEGER :: I,IE,GET_SPECIES_INDEX
      LOGICAL :: IS_BOLSIG_REACTION

      IE=GET_SPECIES_INDEX('E') 
      DO I=1,NREAC
       IF(IS_BOLSIG_REACTION(REACF(I))) THEN
        IS_SPEC_IN_REAC(IE,I)=.TRUE.
       ENDIF
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE SET_NEUTRAL_IF_ANY()
      USE GLOBAL 
      INTEGER :: I,J
      LOGICAL :: IS_ANY_NEUTRAL_REACTION,IS_CHARGED_SPECIES
      
      DO I=1,NREAC
       IF(IS_ANY_NEUTRAL_REACTION(REACF(I))) THEN
        DO J=1,NSPEC
         IF(.NOT.IS_CHARGED_SPECIES(SPEC(J))) THEN
          IS_SPEC_IN_REAC(J,I)=.TRUE.
         ENDIF
        ENDDO
       ENDIF
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
      LOGICAL FUNCTION IS_LE_TO(N,NMAX)
      IMPLICIT NONE
      INTEGER :: N,NMAX
      
      IS_LE_TO=N.LE.NMAX

      END FUNCTION
      !------------------------------------------------------------------
      FUNCTION IS_ANY_NEUTRAL_REACTION(REAC)
      USE GLOBAL, ONLY : NEUTRALID
      IMPLICIT NONE
      CHARACTER(LEN=*) :: REAC
      LOGICAL :: IS_ANY_NEUTRAL_REACTION

      IS_ANY_NEUTRAL_REACTION=INDEX(REAC,NEUTRALID).GT.0

      END FUNCTION
      !-----------------------------------------------------------------
      FUNCTION IS_BOLSIG_REACTION(REAC)
      USE GLOBAL, ONLY : BOLSIGID
      IMPLICIT NONE
      LOGICAL :: IS_BOLSIG_REACTION
      CHARACTER(LEN=*) :: REAC
      
      IF(INDEX(REAC,BOLSIGID).GT.0) THEN
       IS_BOLSIG_REACTION=.TRUE.
      ELSE
       IS_BOLSIG_REACTION=.FALSE.
      ENDIF

      END FUNCTION
      !-----------------------------------------------------------------
      FUNCTION IS_SPEC_IN_REACTION(REAC,SPEC)
      IMPLICIT NONE
      LOGICAL :: CASEA,CASEB,CASEC,IS_SPEC_IN_REACTION
      CHARACTER(LEN=*) :: REAC,SPEC
      CHARACTER(LEN=LEN(TRIM(ADJUSTL(SPEC)))+2) :: CA
      
      CA=' '//TRIM(ADJUSTL(SPEC))//' ' 
      IF(INDEX(REAC,CA).GT.0) THEN
       IS_SPEC_IN_REACTION=.TRUE.
      ELSE
       IS_SPEC_IN_REACTION=.FALSE.
      ENDIF

      END FUNCTION
      !-----------------------------------------------------------------
      FUNCTION IS_CHARGED_SPECIES(SPEC)
      USE GLOBAL, ONLY : EID
      IMPLICIT NONE
      CHARACTER(LEN=*) :: SPEC
      CHARACTER(LEN=LEN(TRIM(ADJUSTL(SPEC)))) :: C
      LOGICAL :: IS_CHARGED_SPECIES

      C=TRIM(ADJUSTL(SPEC))
      IF(C.EQ.EID) THEN
       IS_CHARGED_SPECIES=.TRUE.
      ELSE
       IS_CHARGED_SPECIES=(INDEX(C,'^+').GT.0).OR.(INDEX(C,'^-').GT.0) 
      ENDIF
             
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
      SUBROUTINE SET_IS_NEUTRAL()
      USE GLOBAL, ONLY : NSPEC,SPEC,IS_SPEC_NEUTRAL
      IMPLICIT NONE
      LOGICAL :: IS_CHARGED_SPECIES
      INTEGER :: NC,I
      
      DO I=1,NSPEC
       IF(IS_CHARGED_SPECIES(SPEC(I))) THEN
        IS_SPEC_NEUTRAL(I)=.FALSE.
       ELSE
        IS_SPEC_NEUTRAL(I)=.TRUE.
       ENDIF
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE GET_CHARGED_SPECIES(NSPEC,CSPEC,CCHAR,NC)
      IMPLICIT NONE
      LOGICAL :: IS_CHARGED_SPECIES
      INTEGER :: NSPEC,NC,I,J
      CHARACTER(LEN=*) :: CSPEC(NSPEC),CCHAR(NSPEC)

      J=0
      DO I=1,NSPEC
       IF(IS_CHARGED_SPECIES(CSPEC(I))) THEN
        J=J+1
        CCHAR(J)=CSPEC(I)
       ENDIF
      ENDDO
      NC=J

      RETURN
      END 
      !-----------------------------------------------------------------
      SUBROUTINE SET_FORMATTED_REACTIONS()
      !
      ! AUTHOR: Z. NIKOLAOU  
      !  
      USE GLOBAL
      IMPLICIT NONE
      INTEGER :: I
      CHARACTER(LEN=NSMX) :: C,CWRK

      DO I=1,NREAC
       C=REAC(I)
       CALL REPLACE_TEXT(C,'^+','^POS',CWRK,NSMX)
       C=CWRK
       CALL REPLACE_TEXT(C,':',' ',CWRK,NSMX)
       C=CWRK
       CALL REPLACE_TEXT(C,'^-','^NEG',CWRK,NSMX)
       C=CWRK
       CALL REPLACE_TEXT(C,'+',' ',CWRK,NSMX)
       C=CWRK
       CALL REPLACE_TEXT(C,'=>',' ',CWRK,NSMX)
       C=CWRK
       CALL REPLACE_TEXT(C,'->',' ',CWRK,NSMX)
       C=CWRK
       CALL REPLACE_TEXT(C,'(E-V)',' ',CWRK,NSMX)
       C=CWRK
       CALL REPLACE_TEXT(C,'ANY_NEUTRAL','  ANY_NEUTRAL ',CWRK,NSMX)
       C=CWRK
       CALL REPLACE_TEXT(C,'POS','+',CWRK,NSMX)
       C=CWRK
       CALL REPLACE_TEXT(C,'NEG','-',CWRK,NSMX)
       C='* '//TRIM(ADJUSTL(CWRK))//' *'
       REACF(I)=TRIM(ADJUSTL(C))
      ENDDO
      
      RETURN
      END
      !-----------------------------------------------------------------

