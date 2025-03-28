      !-----------------------------------------------------------------
      SUBROUTINE READ_CHEM(DIR,CHEMFL,SPECFL)
      !      
      ! AUTHOR: Z. NIKOLAOU
      !
      USE ZDPLASKIN_PARSE, ONLY : ZDP_INIT
      IMPLICIT NONE
      CHARACTER(LEN=*) :: DIR,CHEMFL,SPECFL
    
      WRITE(*,*) '***READ_CHEM***'

      CALL REMOVE_TABS_FROM_FILE(DIR//SPECFL)
      CALL REMOVE_TABS_FROM_FILE(DIR//CHEMFL)
      CALL ZDP_INIT(DIR//CHEMFL)
       !SET VARS:NELEM,NSPEC,NREAC,ELEM,SPEC,REAC,REAC_SPEC,RSPEC
      CALL SET_STOICH_COEFFS()
       !SET VARS:NUR,NUP,DELTANU
      CALL CHECK_CHARGE() 
      !TODO:CALL CHECK_STOICHIOMETRY()
      
      WRITE(*,*) '***READ_CHEM***'

      RETURN
      END
      !------------------------------------------------------------------
      SUBROUTINE SET_STOICH_COEFFS()
      USE GLOBAL, ONLY : NUR,NUP,DELTANU,ZERO
      USE ZDPLASKIN_PARSE, ONLY: NSPEC,NREAC,SPEC,REAC
      IMPLICIT NONE
      INTEGER :: I

      NUR(1:NREAC,1:NSPEC)=ZERO
      NUP(1:NREAC,1:NSPEC)=ZERO
      DELTANU(1:NREAC,1:NSPEC)=ZERO
      WRITE(*,*) 'REACTION STOICH. COEFFS:'
      DO I=1,NREAC
       CALL SET_REAC_STOICH_COEFFS(I,NSPEC,NREAC,REAC(I),SPEC)
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE SET_REAC_STOICH_COEFFS(IR,NSPEC,NREAC,REAC,SPEC)
      !TODO: FIX FOR ANY_NEUTRAL REACTIONS 
      USE GLOBAL, ONLY : NSMX,NSPMX,NUR,NUP,DELTANU,EQ_SEP,RSPEC, &
                         IS_ANY_NEUTRAL_REAC
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
     
      WRITE(*,*) '  REACTANTS: ',TRIM(ADJUSTL(RANDP(1))) 
      CALL SET_NU_FROM_LIST(IR,NR,RLIST(1:NR),NUR(1:NREAC,1:NSPEC))
      
      WRITE(*,*) '  PRODUCTS: ',TRIM(ADJUSTL(RANDP(2)))
      CALL SET_NU_FROM_LIST(IR,NP,PLIST(1:NP),NUP(1:NREAC,1:NSPEC))

      !TODO: MOVE TO ZDPLASKIN PARSE?     
      CALL SET_NU_FOR_ANY_NEUTRAL_REAC()
      CALL SET_NU_FOR_ANY_ION_POS_REAC()
      CALL SET_NU_FOR_ANY_ION_NEG_REAC()
      CALL SET_NU_FOR_ANY_SPEC_REAC()

      DELTANU(1:NREAC,1:NSPEC)=NUP(1:NREAC,1:NSPEC)-NUR(1:NREAC,1:NSPEC)

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE SET_NU_FOR_ANY_NEUTRAL_REAC()
      USE GLOBAL, ONLY : NSPEC,NREAC,NUR,NUP,IS_SPEC_CHARGED, &
                         IS_ANY_NEUTRAL_REAC,ONE
      IMPLICIT NONE
      INTEGER :: I,J

      DO I=1,NREAC
       IF(IS_ANY_NEUTRAL_REAC(I)) THEN
        DO J=1,NSPEC
         IF(.NOT.IS_SPEC_CHARGED(J)) THEN
          NUR(I,J)=ONE
          NUP(I,J)=ONE
         ENDIF
        ENDDO
       ENDIF        
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE SET_NU_FOR_ANY_ION_POS_REAC()
      USE GLOBAL, ONLY : NSPEC,NREAC,NUR,NUP,SPEC_CHARGE, &
                         IS_ANY_ION_POS_REAC,ZERO,ONE
      IMPLICIT NONE
      INTEGER :: I,J

      DO I=1,NREAC
       IF(IS_ANY_ION_POS_REAC(I)) THEN
        DO J=1,NSPEC
         IF(SPEC_CHARGE(I).GT.ZERO) THEN
          NUR(I,J)=ONE
          NUP(I,J)=ONE
         ENDIF
        ENDDO
       ENDIF        
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE SET_NU_FOR_ANY_ION_NEG_REAC()
      USE GLOBAL, ONLY : NSPEC,NREAC,NUR,NUP,SPEC_CHARGE, &
                         IS_ANY_ION_NEG_REAC,ZERO,ONE
      IMPLICIT NONE
      INTEGER :: I,J

      DO I=1,NREAC
       IF(IS_ANY_ION_NEG_REAC(I)) THEN
        DO J=1,NSPEC
         IF(SPEC_CHARGE(I).LT.ZERO) THEN
          NUR(I,J)=ONE
          NUP(I,J)=ONE
         ENDIF
        ENDDO
       ENDIF        
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE SET_NU_FOR_ANY_SPEC_REAC()
      USE GLOBAL, ONLY : NSPEC,NREAC,NUR,NUP,SPEC_CHARGE, &
                         IS_ANY_SPEC_REAC,ONE
      IMPLICIT NONE
      INTEGER :: I,J

      DO I=1,NREAC
       IF(IS_ANY_SPEC_REAC(I)) THEN
        DO J=1,NSPEC
         NUR(I,J)=ONE
         NUP(I,J)=ONE
        ENDDO
       ENDIF        
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE SET_NU_FROM_LIST(IR,NL,NO_AND_CHAR_LIST,NU)
      USE GLOBAL, ONLY : NSPEC,NREAC,NSMX,SPEC
      IMPLICIT NONE
      INTEGER :: I,J,IR,NL
      DOUBLE PRECISION :: NUM,NU(NREAC,NSPEC)
      CHARACTER(LEN=NSMX) :: NO_AND_CHAR_LIST(NL),CH
      
      DO I=1,NL
       CALL SPLIT_NUM_AND_CHAR(NO_AND_CHAR_LIST(I),NSMX,NUM,CH)
       DO J=1,NSPEC
        IF(TRIM(ADJUSTL(CH)).EQ.TRIM(ADJUSTL(SPEC(J)))) THEN
         NU(IR,J)=NU(IR,J)+NUM
        ENDIF
       ENDDO
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
      SUBROUTINE CHECK_CHARGE()
      USE GLOBAL, ONLY: NREAC,REAC,ZERO
      IMPLICIT NONE
      DOUBLE PRECISION :: GET_REACTION_CHARGE,CHARGE
      INTEGER :: I

      DO I=1,NREAC
       CHARGE=GET_REACTION_CHARGE(I)
       IF(CHARGE.NE.ZERO) THEN
        WRITE(*,*) '*CHECK_CHARGE: ERROR, CHARGE NOT ZERO FOR REAC:'
        WRITE(*,*) TRIM(ADJUSTL(REAC(I)))
        WRITE(*,*) 'CHARGE=',CHARGE
        WRITE(*,*) 'CHECK CHEM. MECH. TERMINATING ...'
        STOP
       ENDIF
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
      FUNCTION GET_REACTION_CHARGE(I)
      USE GLOBAL, ONLY: NSPEC,NUR,NUP,SPEC,SPEC_CHARGE,ZERO
      IMPLICIT NONE
      INTEGER :: I,J
      DOUBLE PRECISION :: GET_REACTION_CHARGE

      GET_REACTION_CHARGE=ZERO
      DO J=1,NSPEC
       GET_REACTION_CHARGE=(NUR(I,J)-NUP(I,J))*SPEC_CHARGE(J)+ &
                           GET_REACTION_CHARGE
      ENDDO

      END FUNCTION
      !-----------------------------------------------------------------
