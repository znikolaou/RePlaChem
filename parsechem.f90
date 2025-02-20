      SUBROUTINE READ_CHEM(INDIR,RATEDIR,CHEMFL,SPECFL,CREAC,CSPEC, &
                           CREAC_F,NREAC,NSPEC,REAC_SPEC_MATRIX)
      ! AUTHOR: Z. NIKOLAOU
      USE GLOBAL, ONLY : NWRKC,NSTR_SPMX,NSTR_REMX
      USE PRECIS, ONLY: DBL_P
      IMPLICIT NONE
      CHARACTER(LEN=*) :: INDIR,RATEDIR,CHEMFL,SPECFL
      CHARACTER(LEN=NSTR_REMX) :: CREAC(NWRKC),CREAC_F(NWRKC)
      CHARACTER(LEN=NSTR_SPMX) :: CSPEC(NWRKC)
      INTEGER I,I1,I2,NREAC,NSPEC,REAC_SPEC_MATRIX(NWRKC,NWRKC)
      PARAMETER(I1=1,I2=2)

      !REMOVE TABS FROM INPUT FILES
      CALL SYSTEM("sed -i 's/\t/\ /g'"//' '//INDIR//CHEMFL) 
      CALL SYSTEM("sed -i 's/\t/\ /g'"//' '//INDIR//SPECFL) 

      OPEN(UNIT=I1,FILE=INDIR//CHEMFL,STATUS='OLD',FORM='FORMATTED')
      OPEN(UNIT=I2,FILE=INDIR//SPECFL,STATUS='OLD',FORM='FORMATTED')
      CALL PARSE_INDEX_AND_STRING(I1,NSTR_REMX,NWRKC,NREAC,CREAC)
      CALL PARSE_INDEX_AND_STRING(I2,NSTR_SPMX,NWRKC,NSPEC,CSPEC)
      CLOSE(1)
      CLOSE(2)

      WRITE(*,*) 'PARSE_CHEM: SPECIES:' 
      WRITE(*,*) '--------------------'
      DO I=1,NSPEC
       WRITE(*,'(I5,X,A)') I,TRIM(ADJUSTL(CSPEC(I)))
      ENDDO 
      WRITE(*,*) 'PARSE_CHEM: REACTIONS:' 
      WRITE(*,*) '----------------------'
      DO I=1,NREAC
       WRITE(*,'(I5,X,A)') I,TRIM(ADJUSTL(CREAC(I)))
      ENDDO 
        
      STOP
      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE PARSE_INDEX_AND_STRING(FILEID,NMX,NWRKC,NLIST,STRLIST)
      IMPLICIT NONE
      LOGICAL ISEMPTY,ISCOMMENT
      INTEGER :: FILEID,NMX,NWRKC,NLIST
      CHARACTER(LEN=NMX) :: STRLIST(NWRKC),C
      INTEGER :: STAT,I,IDUM,INDXR,IS,INDX_SPACE
      
      I=0
      STAT=0
      DO WHILE(STAT.EQ.0)
       READ(FILEID,'(A)',IOSTAT=STAT) C 
       IF(.NOT.ISEMPTY(C).AND.(.NOT.ISCOMMENT(C))) THEN
        C=TRIM(ADJUSTL(C))
        IS=INDX_SPACE(C)
        I=I+1
        STRLIST(I)=C(IS+1:)
       ENDIF
      ENDDO      
      NLIST=I

      END SUBROUTINE 
      !-----------------------------------------------------------------
      SUBROUTINE SEPARATE_REACTIONS(N,CREAC,REAC,PROD)
      USE GLOBAL, ONLY : NSTR_REMX
      IMPLICIT NONE
      INTEGER :: N,I
      CHARACTER(LEN=*) :: CREAC(N),REAC(N),PROD(N)
      CHARACTER(LEN=NSTR_REMX) :: CLIST(2)

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
      SUBROUTINE GET_REACTION_SPECIES(NSPEC,NREAC,CSPEC,CREAC, &
                                      REAC_SPEC,NOSPEC)
      IMPLICIT NONE
      INTEGER :: NSPEC,NREAC,NOSPEC(NREAC),I,J,K,REAC_SPEC(NREAC,NSPEC)
      CHARACTER(LEN=*) :: CREAC(NREAC),CSPEC(NSPEC)
      LOGICAL :: IS_SPEC_IN_REACTION,IS_BOLSIG_REACTION, &
                 IS_ANY_NEUTRAL_REACTION,IS_CHARGED_SPECIES

      REAC_SPEC(1:NREAC,1:NSPEC)=0
      DO J=1,NREAC
       K=0
       DO I=1,NSPEC
        IF(IS_SPEC_IN_REACTION(CREAC(J),CSPEC(I))) THEN
         K=K+1
         REAC_SPEC(J,I)=1
         NOSPEC(J)=K
        ENDIF
       ENDDO
      ENDDO

      !ADD MANUALLY SPEC 'E' NOT PRESENT FOR BOLSIG REACTIONS
      DO I=1,NREAC
       IF(IS_BOLSIG_REACTION(CREAC(I))) THEN
        REAC_SPEC(I,1)=1 !TODO: THIS IS HARD-CODED HERE! 
        NOSPEC(I)=NOSPEC(I)+1
       ENDIF
      ENDDO

      !CHECK FOR 'ANY_NEUTRAL' REACTIONS
      DO I=1,NREAC
       IF(IS_ANY_NEUTRAL_REACTION(CREAC(I))) THEN
        DO J=1,NSPEC
         IF(.NOT.IS_CHARGED_SPECIES(CSPEC(J))) THEN
          REAC_SPEC(I,J)=1
         ENDIF
        ENDDO
       ENDIF
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
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
      SUBROUTINE REMOVE_DUPLICATE_REACTIONS(NREAC,CREAC,CREAC_F,NF)
      IMPLICIT NONE
      INTEGER :: NREAC,NF
      CHARACTER(LEN=*) :: CREAC(NREAC),CREAC_F(NREAC)

      CALL REMOVE_DUPLICATE_STRINGS(NREAC,CREAC,CREAC_F,NF)
      WRITE(*,*) 'REMOVED',NREAC-NF,'DUPLICATE REACTIONS'

      RETURN
      END  
      !-----------------------------------------------------------------
      SUBROUTINE GET_FORMATTED_REACTIONS(NREAC,CREAC,CREAC_F)
      !
      ! REFORMATS REACTIONS AS BELOW-INTENDED FOR PLASMA CHEMISTRY
      ! 
      ! AUTHOR: Z. NIKOLAOU (2025)     
      USE GLOBAL, ONLY : NSTR_REMX
      IMPLICIT NONE
      INTEGER :: NREAC,I
      CHARACTER(LEN=*) :: CREAC(NREAC)
      CHARACTER(LEN=*) :: CREAC_F(NREAC)
      CHARACTER(LEN=NSTR_REMX) :: C,CWRK

      DO I=1,NREAC
       C=CREAC(I)
       CALL RPLTXTE(C,'^+','^POS',CWRK,NSTR_REMX)
       C=CWRK
       CALL RPLTXTE(C,'^-','^NEG',CWRK,NSTR_REMX)
       C=CWRK
       CALL RPLTXTE(C,'+',' + ',CWRK,NSTR_REMX)
       C=CWRK
       CALL RPLTXTE(C,'=>',' => ',CWRK,NSTR_REMX)
       C=CWRK
       CALL RPLTXTE(C,'->',' -> ',CWRK,NSTR_REMX)
       C=CWRK
       CALL RPLTXTE(C,':',': ',CWRK,NSTR_REMX)
       C=CWRK
       CALL RPLTXTE(C,'(E-V)',' *(E-V)',CWRK,NSTR_REMX)
       C=CWRK
       CALL RPLTXTE(C,'ANY_NEUTRAL','  ANY_NEUTRAL ',CWRK,NSTR_REMX)
       C=CWRK
       CALL RPLTXTE(C,'POS','+',CWRK,NSTR_REMX)
       C=CWRK
       CALL RPLTXTE(C,'NEG','-',CWRK,NSTR_REMX)
       CREAC_F(I)='$ '//TRIM(ADJUSTL(CWRK))//' $'
      ENDDO
      
      RETURN
      END
      !-----------------------------------------------------------------

