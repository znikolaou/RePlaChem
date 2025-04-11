      !-----------------------------------------------------------------
      !
      ! AUTHOR: Z. NIKOLAOU
      !
      !-----------------------------------------------------------------
      MODULE CHEM_PARSE
      USE GLOBAL
      IMPLICIT NONE
      INTEGER, PRIVATE :: NLINES,NREAC_RAW
      CHARACTER(LEN=NSMX), PRIVATE :: LINES(NLINEMX),REAC_F(NREMX), &
                                      CHEM_LINES(NREMX)
      CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
       KEY_ELEM='ELEMENTS', &
       KEY_SPEC='SPECIES', &
       KEY_REAC='REACTIONS', &
       KEY_BOLS='BOLSIG', &
       KEY_END='END', &
       KEY_SET='SET', &
       KEY_NEUTRAL='ANY_NEUTRAL', &
       KEY_ANY_ION_POS='ANY_ION_POSITIVE', & 
       KEY_ANY_ION_NEG='ANY_ION_NEGATIVE', &
       KEY_ANY_SPECIES='ANY_SPECIES', &
       KEY_DUPLICATE='DUPLICATE', &
       KEY_THIRD_BODY_LIND='+M', &
       KEY_THIRD_BODY_TROE='(+M)', &
       KEY_EQ_SEP_F='->', &
       KEY_EQ_SEP_DOUBLE='<=>', &
       KEY_EQ_SEP_SINGLE='=>', &
       KEY_SLASH='/', &
       KEY_AT='@', &
       KEY_EXCL='!', &
       KEY_EQ='=', &
       KEY_DOLLAR='$'

      CONTAINS

      !-----------------------------------------------------------------
      SUBROUTINE CM_INIT(FL)
      IMPLICIT NONE
      INTEGER :: I,J
      CHARACTER(LEN=*) :: FL

      WRITE(*,*)
      WRITE(*,'(A)') '***MODULE CHEM_PARSE***'
      WRITE(*,*) 
      WRITE(*,'(AXA)') 'FILE:',TRIM(ADJUSTL(FL))
      WRITE(*,*) 
      WRITE(*,'(A)') 'PARSING CHEM. MECHANISM ...' 
      WRITE(*,*) 

      CALL REMOVE_TABS_FROM_FILE(FL)

      NELEM=0
      NSPEC=0
      NSPEC_BOLSIG=0
      NREAC=0
      NBOLS_SET=0
      NREAC_DOLLAR=0
      ELEM(1:NSPMX)=' '
      SPEC(1:NSPMX)=' '
      SPEC_BOLSIG(1:NSPMX)=' '
      SPEC_CHARGE(1:NSPMX)=ZERO
      REAC(1:NREMX)=' '
      REAC_CONST(1:NREMX)=' '
      RSPEC(1:NREMX,1:NSPMX)=0
      IS_SPEC_CHARGED(1:NSPMX)=.FALSE.
      IS_ANY_NEUTRAL_REAC(1:NREMX)=.FALSE.
      IS_ANY_SPEC_REAC(1:NREMX)=.FALSE.
      IS_ANY_ION_POS_REAC(1:NREMX)=.FALSE.
      IS_ANY_ION_NEG_REAC(1:NREMX)=.FALSE.
      IS_THIRD_BODY_REAC(1:NREMX)=.FALSE.
      BOLSIG_SEC_SET_LIST=' '
      REAC_SEC_DOLLAR_LIST=' '
      NUR(1:NREMX,1:NSPMX)=ZERO
      NUP(1:NREMX,1:NSPMX)=ZERO
      DELTANU(1:NREMX,1:NSPMX)=ZERO
      THIRD_BODY_SPEC(1:NREMX,1:NSPMX)=' '

      CALL READ_LINES(FL,LINES,NLINES)
      CALL CM_READ_ELEMENTS()
      CALL CM_READ_SPECIES()
      CALL CM_READ_BOLSIG_SPECIES()
      CALL CM_READ_BOLSIG_SEC_SET_LIST()
      CALL CM_READ_AND_FILTER_REACTIONS() 
      CALL CM_READ_REAC_SEC_DOLLAR_LIST()

      WRITE(*,'(A)') 
      WRITE(*,'(A)') 'FINISHED READING FILE'
      WRITE(*,*) 
      WRITE(*,'(A)') 'POST-PROCESSING ...'
      WRITE(*,*)  
      
      CALL CM_SET_IS_SPEC_CHARGED()
      CALL CM_SET_REAC_INFO()
      CALL CM_EXTRACT_REAC_AND_CONST()
      CALL CM_SET_STOICH_COEFFS() !SET VARS:NUR,NUP,DELTANU
      CALL CM_SET_RSPEC()  !SET VARS: RSPEC
      CALL CM_CHECK_CHARGE()  
      !TODO:CALL CHECK_STOICHIOMETRY()
       
      WRITE(*,'(A)') 'EXTRACTED REACTIONS AND RATE CONSTANTS:'
      WRITE(*,'(A)') '---------------------------------------'
      DO I=1,NREAC
       WRITE(*,'(I5XAXAXA)') I,TRIM(ADJUSTL(REAC(I))), &
                                 'RATE CONST:', &
                                 TRIM(ADJUSTL(REAC_CONST(I)))
      ENDDO

      WRITE(*,*) 
      WRITE(*,'(A)') 'STOICHIOMETRIC COEFFS:'
      WRITE(*,'(A)') '----------------------'
      DO I=1,NREAC
       WRITE(*,'(I5XAXAXAXA)') I,TRIM(ADJUSTL(REAC(I))), &
                                 '( '//TRIM(ADJUSTL(REAC_F(I)))//' )'
       WRITE(*,'(XA)') 'REACTANTS:'
       DO J=1,NSPEC
        IF(NUR(I,J).GT.0) THEN
         WRITE(*,*) TRIM(ADJUSTL(SPEC(J))),NUR(I,J)
        ENDIF
       ENDDO
       WRITE(*,'(XA)') 'PRODUCTS:'
       DO J=1,NSPEC
        IF(NUP(I,J).GT.0) THEN
         WRITE(*,*) TRIM(ADJUSTL(SPEC(J))),NUP(I,J)
        ENDIF
       ENDDO
      ENDDO
      WRITE(*,*) 
      WRITE(*,'(A)') 'PARSING COMPLETE!' 
      WRITE(*,*) 
      WRITE(*,'(A)') '***MODULE CHEM_PARSE***'

      RETURN
      END 
      !-----------------------------------------------------------------        
      SUBROUTINE CM_READ_ELEMENTS() 
      IMPLICIT NONE
      INTEGER :: GET_KEY_INDEX,I,IS,IE

      IS=GET_KEY_INDEX(KEY_ELEM,NLINES,LINES(1:NLINES),1)
      IF(IS.GT.0) THEN
       IE=GET_KEY_INDEX(KEY_END,NLINES,LINES(1:NLINES),IS)
       WRITE(*,*)
       WRITE(*,'(AXI5XI5)') 'ELEMENTS SECTION, LINES:',IS,IE
       CALL CM_READ_SECTION(IS,IE,NSPMX,NSMX,ELEM,NELEM,.TRUE.)
       WRITE(*,'(AXI5)') 'ELEMENTS: NELEM=',NELEM
       DO I=1,NELEM
        WRITE(*,'(I4XA)') I,TRIM(ADJUSTL(ELEM(I)))
       ENDDO
      ELSE
       WRITE(*,'(A)') '*ERROR: NO ELEMENTS SECTION FOUND!'// &
                      ' TERMINATING ...'
       STOP
      ENDIF

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE CM_READ_SPECIES() 
      IMPLICIT NONE
      INTEGER :: GET_KEY_INDEX,I,IS,IE 

      IS=GET_KEY_INDEX(KEY_SPEC,NLINES,LINES(1:NLINES),1)
      IF(IS.GT.0) THEN
       IE=GET_KEY_INDEX(KEY_END,NLINES,LINES(1:NLINES),IS)
       WRITE(*,*)
       WRITE(*,'(AXI5XI5)') 'SPECIES SECTION, LINES:',IS,IE
       CALL CM_READ_SECTION(IS,IE,NSPMX,NSMX,SPEC,NSPEC,.TRUE.)
       WRITE(*,'(AXI5)') 'SPECIES: NSPEC=',NSPEC
       DO I=1,NSPEC
        WRITE(*,'(I4XA)') I,TRIM(ADJUSTL(SPEC(I)))
       ENDDO
      ELSE
       WRITE(*,*) '*ERROR: NO SPECIES SECTION FOUND! TERMINATING ...'
       STOP
      ENDIF

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE CM_READ_BOLSIG_SPECIES() 
      IMPLICIT NONE
      INTEGER :: GET_KEY_INDEX,I,IS,IE 

      IS=GET_KEY_INDEX(KEY_BOLS,NLINES,LINES(1:NLINES),1)
      IF(IS.GT.0) THEN
       IE=GET_KEY_INDEX(KEY_END,NLINES,LINES(1:NLINES),IS)
       WRITE(*,*)
       WRITE(*,'(AXI5XI5)') 'BOLSIG SPECIES SECTION, LINES:',IS,IE
       CALL CM_READ_SECTION(IS,IE,NSPMX,NSMX,SPEC_BOLSIG,NSPEC_BOLSIG, &
              .TRUE.)
       WRITE(*,'(AXI5)') 'BOLSIG SPECIES: NSPEC_BOSLIG=',NSPEC_BOLSIG
       DO I=1,NSPEC_BOLSIG
        WRITE(*,'(I4XA)') I,TRIM(ADJUSTL(SPEC_BOLSIG(I)))
       ENDDO
      ELSE
      WRITE(*,'(A)') '*WARNING: NO BOLSIG SPECIES SECTION FOUND!' 
      ENDIF

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE CM_READ_BOLSIG_SEC_SET_LIST() 
      IMPLICIT NONE
      INTEGER :: GET_KEY_INDEX,I,IS,IE 

      IS=GET_KEY_INDEX(KEY_BOLS,NLINES,LINES(1:NLINES),1)
      IF(IS.GT.0) THEN
       IE=GET_KEY_INDEX(KEY_END,NLINES,LINES(1:NLINES),IS) 
       WRITE(*,*)
       WRITE(*,'(AXI5XI5)') 'BOLSIG SECTION, LINES:',IS,IE
       CALL CM_READ_SECTION_SET_LIST(IS,IE,NLINEMX, &
                                     BOLSIG_SEC_SET_LIST,NBOLS_SET)
       WRITE(*,'(A)') 'BOLSIG SECTION SET LINES:'
       IF(NBOLS_SET.GT.0) THEN
        DO I=1,NBOLS_SET
         WRITE(*,'(A)') TRIM(ADJUSTL(BOLSIG_SEC_SET_LIST(I)))
        ENDDO 
       ELSE
        WRITE(*,'(XA)') 'NO SET KEYWORDS FOUND IN BOLSIG SECTION!'
       ENDIF
      ELSE
       WRITE(*,'(A)') '*WARNING: NO BOLSIG SPECIES SECTION FOUND!'
      ENDIF     

      RETURN
      END
      !----------------------------------------------------------------- 
      SUBROUTINE CM_READ_AND_FILTER_REACTIONS() 
      IMPLICIT NONE
      INTEGER :: GET_KEY_INDEX,I,IS,IE 

      IS=GET_KEY_INDEX(KEY_REAC,NLINES,LINES(1:NLINES),1)
      IF(IS.GT.0) THEN
       IE=GET_KEY_INDEX(KEY_END,NLINES,LINES(1:NLINES),IS)
       WRITE(*,*)
       WRITE(*,'(AXI5XI5)') 'REACTIONS SECTION, LINES:',IS,IE
       WRITE(*,'(A)') 'REACTIONS:'
       CALL CM_READ_SECTION(IS,IE,NREMX,NSMX,CHEM_LINES,NREAC_RAW, &
                            .FALSE.)
       CALL CM_FILTER_REACTIONS()
      ELSE
       WRITE(*,'(A)') '*ERROR: NO REACTIONS SECTION FOUND!,'// &
                      'TERMINATING ...'
       STOP
      ENDIF
     
      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE CM_READ_REAC_SEC_DOLLAR_LIST() 
      IMPLICIT NONE
      INTEGER :: GET_KEY_INDEX,I,IS,IE 

      IS=GET_KEY_INDEX(KEY_REAC,NLINES,LINES(1:NLINES),1)
      IE=GET_KEY_INDEX(KEY_END,NLINES,LINES(1:NLINES),IS)
      WRITE(*,*)
      WRITE(*,'(A)') 'REACTIONS SET DOLLAR LIST:'
      CALL CM_READ_SECTION_DOLLAR_LIST(IS,IE,NLINEMX, &
              REAC_SEC_DOLLAR_LIST,NREAC_DOLLAR)
      IF(NREAC_DOLLAR.EQ.0) THEN
       WRITE(*,'(XA)') 'NO $ EXPRESSIONS FOUND IN REACTIONS SECTION!'
      ELSE
       DO I=1,NREAC_DOLLAR
        WRITE(*,'(A)') TRIM(ADJUSTL(REAC_SEC_DOLLAR_LIST(I)))
       ENDDO
      ENDIF
     
      RETURN
      END
      !----------------------------------------------------------------- 
      SUBROUTINE CM_FILTER_REACTIONS()
      IMPLICIT NONE
      INTEGER :: I,NA,NG,ILINE,IREAC,ITHIRD,NTHIRD
      CHARACTER(LEN=NSMX) :: LINE,ATLIST(NSPMX),GLIST(NSPMX,NSPMX)

      IREAC=0
      ILINE=1
      DO WHILE(ILINE.LE.NREAC_RAW)  
       LINE=CHEM_LINES(ILINE)
       IF(CM_IS_GROUP_SPEC_REAC(LINE)) THEN
        
        WRITE(*,'(4XAXA)') '*GROUP REAC FOUND: '//TRIM(ADJUSTL(LINE))
        CALL CM_EXTRACT_REAC_GROUP_SPEC(ILINE,NSPMX,ATLIST,GLIST,NA,NG) 
        !GENERATE NEW REAC FOR EACH GROUP 
        DO I=1,NG
         IREAC=IREAC+1
         CALL CM_GENERATE_REAC_FOR_GROUP(LINE,NA,ATLIST(1:NA), &
                 GLIST(1:NA,I), REAC(IREAC))
         WRITE(*,'(I4XA)') IREAC,TRIM(ADJUSTL(REAC(IREAC)))
        ENDDO
        ILINE=ILINE+NA+1
       
       ELSEIF(CM_IS_THIRD_BODY_LIND_REACTION(LINE)) THEN
        
        IREAC=IREAC+1
        REAC(IREAC)=LINE
        IS_THIRD_BODY_REAC(IREAC)=.TRUE.
        
        WRITE(*,'(A)') 'THIRD BODY (LINDEMAN FORM) REAC FOUND:'
        WRITE(*,'(XI4XA)') IREAC,TRIM(ADJUSTL(REAC(IREAC)))
       
        CALL CM_READ_THIRD_BODY_SPEC(ILINE+1,NTHIRD, &
                                      THIRD_BODY_SPEC(IREAC,1:NSPEC))
      
        ILINE=ILINE+2
      
        ELSEIF(CM_IS_THIRD_BODY_TROE_REACTION(LINE)) THEN
       
        IREAC=IREAC+1
        REAC(IREAC)=LINE
        IS_THIRD_BODY_REAC(IREAC)=.TRUE.
        
        WRITE(*,'(A)') 'THIRD BODY (TROE FORM) REAC FOUND:'
        WRITE(*,'(XI4XA)') IREAC,TRIM(ADJUSTL(REAC(IREAC)))
        
        IF(CM_IS_THREE_TROE_LINES(ILINE)) THEN
         ITHIRD=ILINE+3
         ILINE=ILINE+4
        ELSE
         ITHIRD=ILINE+2
         ILINE=ILINE+3
        ENDIF
        
        CALL CM_READ_THIRD_BODY_SPEC(ITHIRD,NTHIRD, &
                                      THIRD_BODY_SPEC(IREAC,1:NSPEC))
       
       ELSE 

        IREAC=IREAC+1  
        REAC(IREAC)=LINE   
        ILINE=ILINE+1
        WRITE(*,'(I4XA)') IREAC,TRIM(ADJUSTL(REAC(IREAC)))
       
       ENDIF
      
      ENDDO
      NREAC=IREAC

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE CM_GENERATE_REAC_FOR_GROUP(R,NA,ATSPEC,GLIST,RNEW)
      IMPLICIT NONE
      INTEGER :: NA,I
      CHARACTER(LEN=*) :: R,ATSPEC(NA),GLIST(NA),RNEW
      CHARACTER(LEN=LEN(R)) :: ROLD
     
      ROLD=R
      DO I=1,NA
       CALL REPLACE_TEXT(ROLD,TRIM(ADJUSTL(ATSPEC(I))), &
                 TRIM(ADJUSTL(GLIST(I))),RNEW,NSMX)
       ROLD=RNEW
      ENDDO
      
      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE CM_EXTRACT_REAC_GROUP_SPEC(IREAC,NL,ATSPEC,GLIST,NA, &
                                             NG)
      IMPLICIT NONE
      INTEGER :: IREAC,NL,NG,I,K,IK,NA,ILINE
      CHARACTER(LEN=*) :: ATSPEC(NL),GLIST(NL,NL)
      CHARACTER(LEN=NSMX) :: LINE

      !GET @ SPECIES
      CALL GET_SUBTEXTP1_DISTINCT_LIST(CHEM_LINES(IREAC),KEY_AT,NL, &
              ATSPEC,NA) 
      !FOR EACH @ SPECIES FIND LINE WITH GROUP SPEC
      DO I=1,NA
       IK=0
       ILINE=IREAC
       WRITE(*,'(5XAXA)') 'GROUP SPEC:',TRIM(ADJUSTL(ATSPEC(I)))
       DO WHILE(IK.EQ.0)
        ILINE=ILINE+1 
        IK=INDEX(CHEM_LINES(ILINE),TRIM(ADJUSTL(ATSPEC(I))))
       ENDDO
       IF(IK.EQ.0) THEN
        WRITE(*,*) '*EXTRACT_GROUP_SPEC ERROR: NO GROUP SPEC FOUND!'
        STOP
       ENDIF
       LINE=TRIM(ADJUSTL(CHEM_LINES(ILINE)))
       LINE=TRIM(ADJUSTL(LINE(INDEX(LINE,KEY_EQ)+1:)))
       CALL SPLIT_STRING_WITH_SPACES(LINE,NL,GLIST(I,:),NG)
       DO K=1,NG
        WRITE(*,'(6XA)') TRIM(ADJUSTL(GLIST(I,K)))
       ENDDO
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE CM_EXTRACT_REAC_AND_CONST()
      IMPLICIT NONE
      INTEGER :: I,IP,GET_INDEX_FIRST_SPACE
  
      DO I=1,NREAC
       IP=INDEX(REAC(I),KEY_EXCL)
       IF(IP.LE.0) IP=GET_INDEX_FIRST_SPACE(REAC(I))
        IF(IP.LE.0) THEN
         WRITE(*,*) '*ERROR: NO REAC CONST FOUND, TERMINATING ...'
         STOP
        ENDIF
        REAC_CONST(I)=TRIM(ADJUSTL(REAC(I)(IP+1:)))
        REAC(I)(IP:)=' '
      ENDDO

      RETURN
      END     
      !-----------------------------------------------------------------
      SUBROUTINE CM_SET_REAC_INFO()
      IMPLICIT NONE

      CALL CM_SET_NEUTRAL_IF_ANY()
      CALL CM_SET_ION_POS_IF_ANY()
      CALL CM_SET_ION_NEG_IF_ANY()
      CALL CM_SET_ANY_SPEC_IF_ANY()

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE CM_SET_SPECIES_FROM_LIST()
      IMPLICIT NONE 
      INTEGER :: I,J
      
      DO J=1,NREAC
       DO I=1,NSPEC
        IF(CM_IS_SPEC_IN_REACTION(REAC_F(J),SPEC(I))) THEN 
         RSPEC(J,I)=1
        ENDIF
       ENDDO
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE CM_SET_NEUTRAL_IF_ANY()
      INTEGER :: I,J
      
      DO J=1,NREAC
       IF(CM_IS_ANY_NEUTRAL_REACTION(REAC(J))) THEN
        IS_ANY_NEUTRAL_REAC(J)=.TRUE.
       ENDIF
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE CM_SET_ION_POS_IF_ANY()
      INTEGER :: I,J
      
      DO J=1,NREAC
       IF(CM_IS_ANY_ION_POS_REACTION(REAC(J))) THEN
        IS_ANY_ION_POS_REAC(J)=.TRUE.
       ENDIF
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE CM_SET_ION_NEG_IF_ANY()
      INTEGER :: I,J
      
      DO J=1,NREAC
       IF(CM_IS_ANY_ION_NEG_REACTION(REAC(J))) THEN
        IS_ANY_ION_NEG_REAC(J)=.TRUE.
       ENDIF
      ENDDO

      RETURN
      END
      !----------------------------------------------------------------- 
      SUBROUTINE CM_SET_ANY_SPEC_IF_ANY()
      INTEGER :: I,J
      
      DO J=1,NREAC
       IF(CM_IS_ANY_SPEC_REACTION(REAC(J))) THEN
        IS_ANY_SPEC_REAC(J)=.TRUE.
       ENDIF
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE CM_SET_REAC_F() 
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
       CALL REPLACE_TEXT(C,'<=>',' ',CWRK,NSMX)
       C=CWRK
       CALL REPLACE_TEXT(C,'=>',' ',CWRK,NSMX)
       C=CWRK
       CALL REPLACE_TEXT(C,'->',' ',CWRK,NSMX)
       C=CWRK
       CALL REPLACE_TEXT(C,'(E-V)',' ',CWRK,NSMX)
       C=CWRK
       CALL REPLACE_TEXT(C,KEY_NEUTRAL,' ANY_NEUTRAL ',CWRK,NSMX)
       C=CWRK
       CALL REPLACE_TEXT(C,KEY_ANY_ION_POS,' ANY_ION_POSITIVE ',CWRK, &
               NSMX)
       C=CWRK
       CALL REPLACE_TEXT(C,KEY_ANY_ION_NEG,' ANY_ION_NEGATIVE ',CWRK, &
               NSMX)
       C=CWRK
       CALL REPLACE_TEXT(C,KEY_ANY_SPECIES,' ANY_SPECIES ',CWRK,NSMX)
       C=CWRK
       CALL REPLACE_TEXT(C,'POS','+',CWRK,NSMX)
       C=CWRK
       CALL REPLACE_TEXT(C,'NEG','-',CWRK,NSMX)
       C='* '//TRIM(ADJUSTL(CWRK))//' *'
       REAC_F(I)=TRIM(ADJUSTL(C))
      ENDDO
      
      RETURN
      END
      !-----------------------------------------------------------------
      FUNCTION CM_IS_SPEC_IN_REACTION(REAC,SPEC)
      IMPLICIT NONE
      LOGICAL :: CASEA,CASEB,CASEC,CM_IS_SPEC_IN_REACTION
      CHARACTER(LEN=*) :: REAC,SPEC
      CHARACTER(LEN=LEN(TRIM(ADJUSTL(SPEC)))+2) :: CA
      
      CA=' '//TRIM(ADJUSTL(SPEC))//' ' 
      IF(INDEX(REAC,CA).GT.0) THEN
       CM_IS_SPEC_IN_REACTION=.TRUE.
      ELSE
       CM_IS_SPEC_IN_REACTION=.FALSE.
      ENDIF

      END FUNCTION
      !-----------------------------------------------------------------
      SUBROUTINE CM_SET_IS_SPEC_CHARGED()
      IMPLICIT NONE
      INTEGER :: I
      LOGICAL :: IELE,IPOS,INEG

      DO I=1,NSPEC
       IPOS=INDEX(SPEC(I),'^+').GT.0
       INEG=INDEX(SPEC(I),'^-').GT.0
       IELE=(TRIM(ADJUSTL(SPEC(I))).EQ.'E').OR. &
            (TRIM(ADJUSTL(SPEC(I))).EQ.'e')
       IS_SPEC_CHARGED(I)=IPOS.OR.INEG.OR.IELE
       IF(IPOS) SPEC_CHARGE(I)=ONE
       IF(INEG.OR.IELE) SPEC_CHARGE(I)=-ONE
      ENDDO
             
      END
      !-----------------------------------------------------------------
      FUNCTION CM_IS_THREE_TROE_LINES(ILINE)
      IMPLICIT NONE
      INTEGER :: ILINE
      LOGICAL :: CM_IS_THREE_TROE_LINES

      CM_IS_THREE_TROE_LINES= &
              INDEX(CHEM_LINES(ILINE+1),'LOW').GT.0.AND. &
              INDEX(CHEM_LINES(ILINE+2),'TROE').GT.0

      END FUNCTION
      !-----------------------------------------------------------------
      FUNCTION CM_IS_ANY_NEUTRAL_REACTION(REAC)
      IMPLICIT NONE
      CHARACTER(LEN=*) :: REAC
      LOGICAL :: CM_IS_ANY_NEUTRAL_REACTION

      CM_IS_ANY_NEUTRAL_REACTION=INDEX(REAC,KEY_NEUTRAL).GT.0

      END FUNCTION
      !-----------------------------------------------------------------
      FUNCTION CM_IS_ANY_ION_POS_REACTION(REAC)
      IMPLICIT NONE
      CHARACTER(LEN=*) :: REAC
      LOGICAL :: CM_IS_ANY_ION_POS_REACTION

      CM_IS_ANY_ION_POS_REACTION=INDEX(REAC,KEY_ANY_ION_POS).GT.0

      END FUNCTION
      !-----------------------------------------------------------------
      FUNCTION CM_IS_ANY_ION_NEG_REACTION(REAC)
      IMPLICIT NONE
      CHARACTER(LEN=*) :: REAC
      LOGICAL :: CM_IS_ANY_ION_NEG_REACTION

      CM_IS_ANY_ION_NEG_REACTION=INDEX(REAC,KEY_ANY_ION_NEG).GT.0

      END FUNCTION
      !-----------------------------------------------------------------
      FUNCTION CM_IS_ANY_SPEC_REACTION(REAC)
      IMPLICIT NONE
      CHARACTER(LEN=*) :: REAC
      LOGICAL :: CM_IS_ANY_SPEC_REACTION

      CM_IS_ANY_SPEC_REACTION=INDEX(REAC,KEY_ANY_SPECIES).GT.0

      END FUNCTION
      !-----------------------------------------------------------------
      FUNCTION CM_IS_THIRD_BODY_LIND_REACTION(R)
      IMPLICIT NONE
      LOGICAL :: CM_IS_THIRD_BODY_LIND_REACTION
      CHARACTER(LEN=*) :: R

      CM_IS_THIRD_BODY_LIND_REACTION= &
       (INDEX(R,KEY_THIRD_BODY_LIND).NE.0).AND. &
       (.NOT.CM_IS_THIRD_BODY_TROE_REACTION(R))

      END FUNCTION
      !-----------------------------------------------------------------
      FUNCTION CM_IS_THIRD_BODY_TROE_REACTION(R)
      IMPLICIT NONE
      LOGICAL :: CM_IS_THIRD_BODY_TROE_REACTION
      CHARACTER(LEN=*) :: R

      CM_IS_THIRD_BODY_TROE_REACTION=INDEX(R,KEY_THIRD_BODY_TROE).NE.0

      END FUNCTION
      !-----------------------------------------------------------------
      FUNCTION CM_IS_DOLLAR(LINE)
      IMPLICIT NONE
      CHARACTER(LEN=*) :: LINE
      LOGICAL :: CM_IS_DOLLAR

      CM_IS_DOLLAR=INDEX(LINE,KEY_DOLLAR).GT.0

      END FUNCTION
      !-----------------------------------------------------------------
      FUNCTION CM_IS_BOLSIG_REACTION(REAC)
      IMPLICIT NONE
      LOGICAL :: CM_IS_BOLSIG_REACTION
      CHARACTER(LEN=*) :: REAC
      
      IF(INDEX(REAC,KEY_BOLS).GT.0) THEN
       CM_IS_BOLSIG_REACTION=.TRUE.
      ELSE
       CM_IS_BOLSIG_REACTION=.FALSE.
      ENDIF

      END FUNCTION
      !-----------------------------------------------------------------
      FUNCTION CM_IS_GROUP_SPEC_REAC(R)
      IMPLICIT NONE
      LOGICAL :: CM_IS_GROUP_SPEC_REAC
      CHARACTER(LEN=*) :: R

      CM_IS_GROUP_SPEC_REAC=INDEX(R,KEY_AT).NE.0

      END FUNCTION
      !-----------------------------------------------------------------
      FUNCTION CM_IS_DUPLICATE(R)
      IMPLICIT NONE
      LOGICAL :: CM_IS_DUPLICATE
      CHARACTER(LEN=*) :: R

      CM_IS_DUPLICATE=INDEX(R,KEY_DUPLICATE).NE.0

      END FUNCTION
      !-----------------------------------------------------------------
      FUNCTION CM_ISKEYWORD(LINE)
      IMPLICIT NONE
      CHARACTER(LEN=*) :: LINE
      LOGICAL :: CM_ISKEYWORD

      CM_ISKEYWORD=INDEX(LINE,KEY_SET).GT.0

      END FUNCTION
      !-----------------------------------------------------------------
      SUBROUTINE CM_READ_SECTION(LSTART,LEND,NAR_MX,NSTR_MX,CARR,NAR, &
                      IS_SPLIT) 
      IMPLICIT NONE
      INTEGER :: LSTART,LEND,NAR_MX,NSTR_MX,NAR,GET_KEY_INDEX,IS,IE, &
              I,J,JOFS,NA 
      CHARACTER(LEN=NSTR_MX) :: CARR(NAR_MX)
      CHARACTER(LEN=NSMX) :: COLMS(NAR_MX)
      LOGICAL :: ISCOMMENT,ISEMPTY,IS_SPLIT

      !WRITE(*,'(AXI6XAI6)') 'LINES:',LSTART,'-',LEND
      JOFS=0
      NAR=0
      DO I=LSTART+1,LEND-1
       IF(.NOT.ISCOMMENT(LINES(I))) THEN
        IF(.NOT.ISEMPTY(LINES(I))) THEN
         IF(.NOT.CM_ISKEYWORD(LINES(I))) THEN
          IF(.NOT.CM_IS_DUPLICATE(LINES(I))) THEN
           IF(IS_SPLIT) THEN       
            CALL SPLIT_STRING_WITH_SPACES(LINES(I),NSPMX,COLMS,NA)  
            DO J=1,NA
             CARR(J+JOFS)=TRIM(ADJUSTL(COLMS(J)))
            ENDDO
            JOFS=JOFS+NA
            NAR=NAR+NA
            ELSE
            NAR=NAR+1
            CARR(NAR)=TRIM(ADJUSTL(LINES(I)))
           ENDIF
          ENDIF
         ENDIF 
        ENDIF
       ENDIF
      ENDDO

      RETURN
      END 
      !-----------------------------------------------------------------
      SUBROUTINE CM_READ_SECTION_SET_LIST(LSTART,LEND,NL,LIST,NAL) 
      IMPLICIT NONE
      INTEGER :: LSTART,LEND,NL,NAL,I,J
      CHARACTER(LEN=NSMX) :: LIST(NL)
      LOGICAL :: ISCOMMENT,ISEMPTY

      NAL=0
      J=0
      DO I=LSTART+1,LEND-1
       IF(.NOT.ISCOMMENT(LINES(I))) THEN
        IF(.NOT.ISEMPTY(LINES(I))) THEN
         IF(CM_ISKEYWORD(LINES(I))) THEN
          J=J+1
          LIST(J)=TRIM(ADJUSTL(LINES(I)))
         ENDIF
        ENDIF
       ENDIF
      ENDDO
      NAL=J
      IF(NAL.GT.NL) THEN
       WRITE(*,'(A)') '*ERROR: CM_READ_SECTION_SET_LIST'
       WRITE(*,'(A)') 'NAL>NL, TERMINATING ...'
       STOP
      ENDIF

      RETURN
      END 
      !-----------------------------------------------------------------
      SUBROUTINE CM_READ_THIRD_BODY_SPEC(ILINE,NTHIRD,THIRD_BODY_SPEC)
      IMPLICIT NONE
      INTEGER :: IREAC,ILINE,NTHIRD,JSLASH,I,J
      CHARACTER(LEN=NSMX) :: COLMS(2*NSPEC),THIRD_BODY_SPEC(NSPEC)
           
      CALL SPLIT_STRING(CHEM_LINES(ILINE),KEY_SLASH,2*NSPEC,COLMS, &
                                    NTHIRD)
      
      J=0
      I=1
      DO WHILE(TRIM(ADJUSTL(COLMS(I))).NE.' ')
       J=J+1
       THIRD_BODY_SPEC(J)=TRIM(ADJUSTL(COLMS(I)))
       I=I+2
      ENDDO
      NTHIRD=J
      WRITE(*,'(XA)') 'THIRD BODY SPECIES:'
      DO I=1,NTHIRD
       WRITE(*,'(I3XA)') I,TRIM(ADJUSTL(THIRD_BODY_SPEC(I)))
      ENDDO 

 
      END
      !-----------------------------------------------------------------
      SUBROUTINE CM_READ_SECTION_DOLLAR_LIST(LSTART,LEND,NL,LIST,NAL) 
      IMPLICIT NONE
      INTEGER :: LSTART,LEND,NL,NAL,I,J
      CHARACTER(LEN=NSMX) :: LIST(NL)
      LOGICAL :: ISCOMMENT,ISEMPTY

      NAL=0
      J=0
      DO I=LSTART+1,LEND-1
       IF(CM_IS_DOLLAR(LINES(I))) THEN
         J=J+1
         LIST(J)=TRIM(ADJUSTL(LINES(I)))
       ENDIF
      ENDDO
      NAL=J
      IF(NAL.GT.NL) THEN
       WRITE(*,'(A)') '*ERROR: CM_READ_SECTION_DOLLAR_LIST'
       WRITE(*,'(A)') 'NAL>NL, TERMINATING ...'
       STOP
      ENDIF

      RETURN
      END 
      !-----------------------------------------------------------------
      FUNCTION CM_GET_LINE_NO(KEY,ISTART)
      IMPLICIT NONE
      CHARACTER(LEN=*) KEY
      INTEGER :: I,ISTART,CM_GET_LINE_NO

      DO I=ISTART,NLINES
       IF(TRIM(ADJUSTL(LINES(I))).EQ.KEY) EXIT
      ENDDO
      CM_GET_LINE_NO=I

      END FUNCTION
      !----------------------------------------------------------------- 
      SUBROUTINE CM_SET_STOICH_COEFFS()
      IMPLICIT NONE
      INTEGER :: I

      DO I=1,NREAC
       CALL CM_SET_REAC_STOICH_COEFFS(I)
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE CM_SET_RSPEC()
      IMPLICIT NONE
      INTEGER :: I,J

      DO I=1,NREAC
       DO J=1,NSPEC
        IF(NUR(I,J).NE.ZERO.OR.NUP(I,J).NE.ZERO) THEN 
         RSPEC(I,J)=1
        ELSE
         RSPEC(I,J)=0
        ENDIF
       ENDDO
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE CM_SET_REAC_STOICH_COEFFS(IR)
      IMPLICIT NONE
      INTEGER :: IR,NA,NR,NP,I,J
      CHARACTER(LEN=NSMX) :: RANDP(2),RLIST(NSPMX),PLIST(NSPMX),CH
      DOUBLE PRECISION :: NUM

      CALL CM_FORMAT_REACTION(IR)
      CALL SPLIT_STRING(REAC_F(IR),KEY_EQ_SEP_F,2,RANDP,NA)
      CALL SPLIT_STRING_WITH_SPACES(RANDP(1),NSPMX,RLIST,NR)
      CALL SPLIT_STRING_WITH_SPACES(RANDP(2),NSPMX,PLIST,NP)
     
      CALL CM_SET_NU_FROM_LIST(IR,NR,RLIST(1:NR),NUR(1:NREAC,1:NSPEC))
      CALL CM_SET_NU_FROM_LIST(IR,NP,PLIST(1:NP),NUP(1:NREAC,1:NSPEC))
      CALL CM_SET_NU_FOR_ANY_NEUTRAL_REAC()
      CALL CM_SET_NU_FOR_ANY_ION_POS_REAC()
      CALL CM_SET_NU_FOR_ANY_ION_NEG_REAC()
      CALL CM_SET_NU_FOR_ANY_SPEC_REAC()
      CALL CM_SET_NU_FOR_THIRD_BODY_REAC()

      DELTANU(1:NREAC,1:NSPEC)=NUP(1:NREAC,1:NSPEC)-NUR(1:NREAC,1:NSPEC)

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE CM_SET_NU_FOR_ANY_NEUTRAL_REAC()
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
      SUBROUTINE CM_SET_NU_FOR_ANY_ION_POS_REAC()
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
      SUBROUTINE CM_SET_NU_FOR_ANY_ION_NEG_REAC()
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
      SUBROUTINE CM_SET_NU_FOR_ANY_SPEC_REAC()
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
      SUBROUTINE CM_SET_NU_FOR_THIRD_BODY_REAC()
      IMPLICIT NONE
      INTEGER :: I,J,INDX

      DO I=1,NREAC
       IF(IS_THIRD_BODY_REAC(I)) THEN
        J=1       
        WRITE(*,*) TRIM(ADJUSTL(REAC(I)))
        DO WHILE(THIRD_BODY_SPEC(I,J).NE.' ')
         INDX=CM_GET_SPECIES_INDEX(THIRD_BODY_SPEC(I,J))
         NUR(I,INDX)=ONE
         NUP(I,INDX)=ONE
         J=J+1
        ENDDO
       ENDIF        
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE CM_SET_NU_FROM_LIST(IR,NL,NO_AND_CHAR_LIST,NU)
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
      SUBROUTINE CM_FORMAT_REACTION(I)
      IMPLICIT NONE
      INTEGER :: I
      CHARACTER(LEN=NSMX) :: CWRK,C
    
      CALL REPLACE_TEXT(REAC(I),'^+','^POS',CWRK,NSMX)
      C=CWRK
      CALL REPLACE_TEXT(C,'^-','^NEG',CWRK,NSMX)
      C=CWRK
      CALL REPLACE_TEXT(C,KEY_EQ_SEP_DOUBLE,' -> ',CWRK,NSMX)
      C=CWRK
      CALL REPLACE_TEXT(C,KEY_EQ_SEP_SINGLE,' -> ',CWRK,NSMX)
      C=CWRK
      CALL REPLACE_TEXT(C,'(+M)','+(M)',CWRK,NSMX)
      C=CWRK
      CALL REPLACE_TEXT(C,'+',' ',CWRK,NSMX)
      C=CWRK
      CALL REPLACE_TEXT(C,'NEG','-',CWRK,NSMX)
      C=CWRK
      CALL REPLACE_TEXT(C,'POS','+',CWRK,NSMX)
      REAC_F(I)=' '//TRIM(ADJUSTL(CWRK))//' '
      
      RETURN
      END
      !-----------------------------------------------------------------
      FUNCTION CM_GET_SPECIES_INDEX(SP)
      IMPLICIT NONE
      CHARACTER(LEN=*) :: SP
      INTEGER :: I,CM_GET_SPECIES_INDEX

      CM_GET_SPECIES_INDEX=0
      DO I=1,NSPEC
       IF(TRIM(ADJUSTL(SP)).EQ.TRIM(ADJUSTL(SPEC(I)))) THEN
         CM_GET_SPECIES_INDEX=I
         EXIT
       ENDIF
      ENDDO
      
      END FUNCTION
      !-----------------------------------------------------------------
      SUBROUTINE CM_CHECK_CHARGE()
      IMPLICIT NONE
      DOUBLE PRECISION :: CHARGE
      INTEGER :: I

      DO I=1,NREAC
       CHARGE=CM_GET_REACTION_CHARGE(I)
       IF(CHARGE.NE.ZERO) THEN
        WRITE(*,*) '*CM_CHECK_CHARGE: ERROR, CHARGE NOT ZERO FOR REAC:'
        WRITE(*,*) TRIM(ADJUSTL(REAC(I)))
        WRITE(*,*) 'CHARGE=',CHARGE
        WRITE(*,*) 'CHECK CHEM. MECH. TERMINATING ...'
        STOP
       ENDIF
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
      FUNCTION CM_GET_REACTION_CHARGE(I)
      IMPLICIT NONE
      INTEGER :: I,J
      DOUBLE PRECISION :: CM_GET_REACTION_CHARGE

      CM_GET_REACTION_CHARGE=ZERO
      DO J=1,NSPEC
       CM_GET_REACTION_CHARGE=(NUR(I,J)-NUP(I,J))*SPEC_CHARGE(J)+ &
                           CM_GET_REACTION_CHARGE
      ENDDO

      END FUNCTION
      !-----------------------------------------------------------------

      END MODULE
      !-----------------------------------------------------------------
