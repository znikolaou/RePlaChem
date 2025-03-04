      !-----------------------------------------------------------------
      SUBROUTINE READ_CHEM(DIR,CHEMFL,SPECFL)
      !      
      ! AUTHOR: Z. NIKOLAOU
      !
      USE GLOBAL
      USE PRECIS, ONLY : DBL_P
      IMPLICIT NONE
      CHARACTER(LEN=*) :: DIR,CHEMFL,SPECFL
      CHARACTER(LEN=LEN(DIR)+6) :: RDIR
      CHARACTER(LEN=NSMX) :: CREAC(NWRK),CSPEC(NWRK),CREAC_F(NWRK)
      INTEGER I,I1,I2,NOSPEC(NWRK)
      PARAMETER(I1=1,I2=2)
      LOGICAL :: IS_LE_TO,IS_SPEC(NWRK,NWRK)

      RDIR=TRIM(ADJUSTL(DIR))//'rates/'

      CREAC(1:NWRK)=''
      CSPEC(1:NWRK)=''
      CREAC_F(1:NWRK)=''

      !REMOVE TABS FROM INPUT FILES
      CALL SYSTEM("sed -i 's/\t/\ /g'"//' '//DIR//CHEMFL) 
      CALL SYSTEM("sed -i 's/\t/\ /g'"//' '//DIR//SPECFL) 

      CALL READ_INT_AND_STR(DIR//CHEMFL,REAC,NREAC)   
      CALL READ_INT_AND_STR(DIR//SPECFL,SPEC,NSPEC)
     
      IF(.NOT.IS_LE_TO(NREAC,NWRK)) THEN
       WRITE(*,*) '*ERROR: NREAC>NWRK:',NREAC,NWRK
       STOP
      ENDIF
      IF(.NOT.IS_LE_TO(NSPEC,NWRK)) THEN
       WRITE(*,*) '*ERROR: NSPEC>NWRK:',NSPEC,NWRK
       STOP
      ENDIF

      WRITE(*,*) '***READ_CHEM***'
      WRITE(*,*) 'SPECIES:',NSPEC 
      DO I=1,NSPEC
       WRITE(*,'(I5,X,A)') I,TRIM(ADJUSTL(CSPEC(I)))
      ENDDO 
      WRITE(*,*) 'REACTIONS:',NREAC 
      DO I=1,NREAC
       WRITE(*,'(I5,X,A)') I,TRIM(ADJUSTL(CREAC(I)))
      ENDDO 
     
      CALL GET_FORMATTED_REACTIONS()
      CALL GET_REACTION_SPECIES()
     
      WRITE(*,*) '***READ_CHEM***'

      RETURN
      END
      !------------------------------------------------------------------
      SUBROUTINE READ_INT_AND_STR(FL,STRL,NL)
      USE GLOBAL, ONLY : NSMX,NWRK
      CHARACTER(LEN=*) :: FL 
      INTEGER :: N,ID,STAT,I,NL,NA
      CHARACTER(LEN=NSMX) :: STRL(NWRK),CL(2),C
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
      SUBROUTINE GET_REAC_SPEC(NSPEC,NREAC,CREAC,CSPEC,CREAC_SPEC,NC)
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
       CALL GET_FORMATTED_REACTIONS()
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
      SUBROUTINE GET_REACTION_SPECIES()
      USE GLOBAL
      IMPLICIT NONE
      INTEGER :: I,J
      LOGICAL :: IS_SPEC_IN_REACTION,IS_BOLSIG_REACTION, &
                 IS_ANY_NEUTRAL_REACTION,IS_CHARGED_SPECIES

      IS_SPEC_IN_REAC(1:NWRK,1:NWRK)=.FALSE.
      DO J=1,NREAC
       DO I=1,NSPEC
        IF(IS_SPEC_IN_REACTION(REACF(J),SPEC(I))) THEN 
         IS_SPEC_IN_REAC(I,J)=.TRUE.
        ENDIF
       ENDDO
      ENDDO

      !TODO :: IMPLEMENT NEXT TWO AS SEPARATE SUBROUTINES 
      ! TO KEEP ABOVE CODE GENERAL! 
      !ADD MANUALLY SPEC 'E' NOT PRESENT FOR BOLSIG REACTIONS
      DO I=1,NREAC
       IF(IS_BOLSIG_REACTION(REACF(I))) THEN
        IS_SPEC_IN_REAC(1,I)=.TRUE. !TODO: THIS IS HARD-CODED HERE! 
       ENDIF
      ENDDO

      !CHECK FOR 'ANY_NEUTRAL' REACTIONS
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
      IMPLICIT NONE
      CHARACTER(LEN=*) :: REAC,NEUTRAL
      LOGICAL :: IS_ANY_NEUTRAL_REACTION
      PARAMETER(NEUTRAL='ANY_NEUTRAL')

      IS_ANY_NEUTRAL_REACTION=INDEX(REAC,NEUTRAL).GT.0

      END FUNCTION
      !-----------------------------------------------------------------
      FUNCTION IS_BOLSIG_REACTION(REAC)
      IMPLICIT NONE
      LOGICAL :: IS_BOLSIG_REACTION
      CHARACTER(LEN=*) :: REAC
      CHARACTER(LEN=*) :: BOLSIG
      PARAMETER(BOLSIG='bolsig')
      
      IF(INDEX(REAC,BOLSIG).GT.0) THEN
       IS_BOLSIG_REACTION=.TRUE.
      ELSE
       IS_BOLSIG_REACTION=.FALSE.
      ENDIF

      END FUNCTION
      !-----------------------------------------------------------------
      FUNCTION IS_FOREIGN_SPEC_IN_REACTION(CREAC,NSPEC,CSPEC)
      IMPLICIT NONE
      INTEGER :: NSPEC,I,REAC_SPEC(1,NSPEC)
      CHARACTER(LEN=*) :: CREAC,CSPEC(NSPEC)
      LOGICAL :: IS_FOREIGN_SPEC_IN_REACTION

      IS_FOREIGN_SPEC_IN_REACTION=.FALSE.
      
      !CALL GET_REACTION_SPECIES(NSPEC,1,CREAC,CSPEC,REAC_SPEC, &
      !                          CREAC_SPEC,NOSPEC)
      

      !ENDDO

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
      IMPLICIT NONE
      CHARACTER(LEN=*) :: SPEC
      CHARACTER(LEN=LEN(TRIM(ADJUSTL(SPEC)))) :: C
      LOGICAL :: IS_CHARGED_SPECIES

      C=TRIM(ADJUSTL(SPEC))
      IS_CHARGED_SPECIES=(INDEX(C,'^+').GT.0).OR. &
                         (INDEX(C,'^-').GT.0).OR.(C.EQ.'E')
             
      END FUNCTION
      !-----------------------------------------------------------------
      SUBROUTINE REMOVE_REACTIONS_WITH_FOREIGN_SPECIES(NREAC,CREAC, &
                      NSPEC,CSPEC,CREAC_F,NF)
      IMPLICIT NONE
      INTEGER :: NREAC,NSPEC,NF,I,J
      CHARACTER(LEN=*) :: CREAC(NREAC),CSPEC(NSPEC),CREAC_F(NREAC)

      J=0
      DO I=1,NREAC
       !IF(.NOT.IS_FOREIGN_SPECIES_PRESENT(CREAC(I),NSPEC,CSPEC)) THEN
       ! J=J+1
       ! CREAC_F(J)=CREAC(I)
       !ENDIF
      ENDDO
      NF=I

      RETURN
      END
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
      SUBROUTINE GET_NEUTRAL_SPECIES(NSPEC,CSPEC,CNEUT,NC)
      IMPLICIT NONE
      LOGICAL :: IS_CHARGED_SPECIES
      INTEGER :: NSPEC,NC,I,J
      CHARACTER(LEN=*) :: CSPEC(NSPEC),CNEUT(NSPEC)

      J=0
      DO I=1,NSPEC
       IF(.NOT.IS_CHARGED_SPECIES(CSPEC(I))) THEN
        J=J+1
        CNEUT(J)=CSPEC(I)
       ENDIF
      ENDDO
      NC=J

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
      SUBROUTINE GET_FORMATTED_REACTIONS()
      !
      ! AUTHOR: Z. NIKOLAOU  
      !  
      USE GLOBAL
      IMPLICIT NONE
      INTEGER :: I
      CHARACTER(LEN=NSMX) :: C,CWRK

      DO I=1,NREAC
       C=REAC(I)
       CALL RPLTXTE(C,'^+','^POS',CWRK,NSMX)
       C=CWRK
       CALL RPLTXTE(C,':',' ',CWRK,NSMX)
       C=CWRK
       CALL RPLTXTE(C,'^-','^NEG',CWRK,NSMX)
       C=CWRK
       CALL RPLTXTE(C,'+',' ',CWRK,NSMX)
       C=CWRK
       CALL RPLTXTE(C,'=>',' ',CWRK,NSMX)
       C=CWRK
       CALL RPLTXTE(C,'->',' ',CWRK,NSMX)
       C=CWRK
       CALL RPLTXTE(C,'(E-V)',' ',CWRK,NSMX)
       C=CWRK
       CALL RPLTXTE(C,'ANY_NEUTRAL','  ANY_NEUTRAL ',CWRK,NSMX)
       C=CWRK
       CALL RPLTXTE(C,'POS','+',CWRK,NSMX)
       C=CWRK
       CALL RPLTXTE(C,'NEG','-',CWRK,NSMX)
       C='* '//TRIM(ADJUSTL(CWRK))//' *'
       REACF(I)=TRIM(ADJUSTL(C))
      ENDDO
      
      RETURN
      END
      !-----------------------------------------------------------------

