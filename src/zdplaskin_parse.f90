      !-----------------------------------------------------------------
      MODULE ZDPLASKIN_PARSE
      USE GLOBAL, ONLY : NELEM,NSPEC,NSPEC_BOLSIG,NREAC,ELEM,SPEC, &
                         SPEC_BOLSIG,REAC,REAC_CONST,IS_SPEC_CHARGED, &
                         REAC_SPEC,RSPEC, &
                         NSFLMX,NSMX,NLINEMX,NSPMX,NREMX,NSMX
      IMPLICIT NONE
      CHARACTER(LEN=*), PARAMETER, PRIVATE :: KEY_ELEM='ELEMENTS', &
       KEY_SPEC='SPECIES', KEY_REAC='REACTIONS',KEY_BOLS='BOLSIG', &
       KEY_END='END', KEY_SET='SET', KEY_AT='@',KEY_EXCL='!', &
       KEY_EQ='=',KEY_NEUTRAL='ANY_NEUTRAL'
      CHARACTER(LEN=NSMX), PRIVATE :: LINES(NLINEMX)
      CHARACTER(LEN=NSMX), PRIVATE :: CHEM_LINES(NREMX)
      CHARACTER(LEN=NSMX) :: REAC_F(NREMX)
      INTEGER, PRIVATE :: NLINES,NREAC_RAW
      
      CONTAINS

      !-----------------------------------------------------------------
      SUBROUTINE ZDP_INIT(FL)
      IMPLICIT NONE
      INTEGER :: I,J
      CHARACTER(LEN=*) :: FL
      
      WRITE(*,*) '***MODULE ZPDLASKIN_PARSE***'

      CALL REMOVE_TABS_FROM_FILE(FL)

      NSPEC=0
      NREAC=0
      NELEM=0
      ELEM(1:NSPMX)=' '
      SPEC(1:NSPMX)=' '
      SPEC_BOLSIG(1:NSPMX)=' '
      REAC(1:NREMX)=' '
      REAC_CONST(1:NREMX)=' '
      REAC_SPEC(1:NREMX,1:NSPMX)=' '
      IS_SPEC_CHARGED(1:NSPMX)=.FALSE.
      
      !READ 
      CALL READ_LINES(FL,LINES,NLINES)
      CALL ZDP_READ_ELEMENTS()
      CALL ZDP_READ_SPECIES()
      CALL ZDP_READ_BOLSIG_SPECIES()
      CALL ZDP_READ_AND_FILTER_REACTIONS() 

      !POST-PROCESS
      CALL ZDP_SET_IS_SPEC_CHARGED()
      CALL ZDP_EXTRACT_REAC_AND_CONST()
      CALL ZDP_SET_REAC_SPEC()

      WRITE(*,*) 'IS SPEC CHARGED:'
      DO I=1,NSPEC
       WRITE(*,*) I,TRIM(ADJUSTL(SPEC(I))), &
                         IS_SPEC_CHARGED(I)
      ENDDO
      WRITE(*,*) 'EXTRACTED REACTIONS AND RATE CONSTANTS:'
      DO I=1,NREAC
       WRITE(*,'(I5XAXAXA)') I,TRIM(ADJUSTL(REAC(I))),'RATE CONST:', &
                           TRIM(ADJUSTL(REAC_CONST(I)))
      ENDDO
      WRITE(*,*) 'REACTION SPECIES'
      DO I=1,NREAC
       WRITE(*,'(I5XAXAXA)') I,TRIM(ADJUSTL(REAC(I)))
       DO J=1,NSPEC
        IF(RSPEC(I,J).EQ.1) WRITE(*,*) TRIM(ADJUSTL(SPEC(J)))
       ENDDO
      ENDDO

      WRITE(*,*) 
      WRITE(*,*) '***MODULE ZPDLASKIN_PARSE***'

      RETURN
      END 
      !-----------------------------------------------------------------        
      SUBROUTINE ZDP_READ_ELEMENTS() 
      IMPLICIT NONE
      INTEGER :: GET_KEY_INDEX,I,IS,IE

      IS=GET_KEY_INDEX(KEY_ELEM,NLINES,LINES(1:NLINES),1)
      IE=GET_KEY_INDEX(KEY_END,NLINES,LINES(1:NLINES),IS)
      CALL ZDP_READ_SECTION(IS,IE,NSPMX,NSMX,ELEM,NELEM,.TRUE.)
      WRITE(*,'(AXI4)') 'ELEMENTS: NELEM=',NELEM
      DO I=1,NELEM
       WRITE(*,*) I,TRIM(ADJUSTL(ELEM(I)))
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE ZDP_READ_SPECIES() 
      IMPLICIT NONE
      INTEGER :: GET_KEY_INDEX,I,IS,IE 

      IS=GET_KEY_INDEX(KEY_SPEC,NLINES,LINES(1:NLINES),1)
      IE=GET_KEY_INDEX(KEY_END,NLINES,LINES(1:NLINES),IS)
      CALL ZDP_READ_SECTION(IS,IE,NSPMX,NSMX,SPEC,NSPEC,.TRUE.)
      WRITE(*,'(AXI4)') 'SPECIES: NSPEC=',NSPEC
      DO I=1,NSPEC
       WRITE(*,*) I,TRIM(ADJUSTL(SPEC(I)))
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE ZDP_READ_BOLSIG_SPECIES() 
      IMPLICIT NONE
      INTEGER :: GET_KEY_INDEX,I,IS,IE 

      IS=GET_KEY_INDEX(KEY_BOLS,NLINES,LINES(1:NLINES),1)
      IE=GET_KEY_INDEX(KEY_END,NLINES,LINES(1:NLINES),IS)
      CALL ZDP_READ_SECTION(IS,IE,NSPMX,NSMX,SPEC_BOLSIG,NSPEC_BOLSIG, &
              .TRUE.)
      WRITE(*,'(AXI4)') 'BOLSIG SPECIES: NSPEC_BOSLIG=',NSPEC_BOLSIG
      DO I=1,NSPEC_BOLSIG
       WRITE(*,*) I,TRIM(ADJUSTL(SPEC_BOLSIG(I)))
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE ZDP_READ_AND_FILTER_REACTIONS() 
      IMPLICIT NONE
      INTEGER :: GET_KEY_INDEX,I,IS,IE 

      IS=GET_KEY_INDEX(KEY_REAC,NLINES,LINES(1:NLINES),1)
      IE=GET_KEY_INDEX(KEY_END,NLINES,LINES(1:NLINES),IS)
      WRITE(*,'(A)') 'REACTIONS:'
      CALL ZDP_READ_SECTION(IS,IE,NREMX,NSMX,CHEM_LINES,NREAC_RAW,.FALSE.)
      CALL ZDP_FILTER_REACTIONS()
     
      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE ZDP_FILTER_REACTIONS()
      IMPLICIT NONE
      INTEGER :: I,NA,NG,ILINE,IREAC
      CHARACTER(LEN=NSMX) :: R,ATLIST(NSPMX),GLIST(NSPMX,NSPMX)

      IREAC=0
      ILINE=1
      DO WHILE(ILINE.LE.NREAC_RAW)  
       R=CHEM_LINES(ILINE)
       IF(ZDP_IS_GROUP_SPEC_REAC(R)) THEN
        WRITE(*,'(4XAXA)') '*GROUP REAC FOUND: '//TRIM(ADJUSTL(R))
        CALL ZDP_EXTRACT_REAC_GROUP_SPEC(ILINE,NSPMX,ATLIST,GLIST,NA,NG) 
        !GENERATE NEW REAC FOR EACH GROUP 
        DO I=1,NG
         IREAC=IREAC+1
         CALL ZDP_GENERATE_REAC_FOR_GROUP(R,NA,ATLIST(1:NA), &
                 GLIST(1:NA,I), REAC(IREAC))
         WRITE(*,'(I4XA)') IREAC,TRIM(ADJUSTL(REAC(IREAC)))
        ENDDO
        ILINE=ILINE+NA+1
       ELSE   
        IREAC=IREAC+1  
        REAC(IREAC)=R   
        ILINE=ILINE+1
        WRITE(*,'(I4XA)') IREAC,TRIM(ADJUSTL(REAC(IREAC)))
       ENDIF
      ENDDO
      NREAC=IREAC

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE ZDP_GENERATE_REAC_FOR_GROUP(R,NA,ATSPEC,GLIST,RNEW)
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
      SUBROUTINE ZDP_EXTRACT_REAC_GROUP_SPEC(IREAC,NL,ATSPEC,GLIST,NA, &
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
      SUBROUTINE ZDP_EXTRACT_REAC_AND_CONST()
      IMPLICIT NONE
      INTEGER :: I,IP
  
      DO I=1,NREAC
       IP=INDEX(REAC(I),KEY_EXCL)
       IF(IP.GT.0) THEN
        REAC_CONST(I)=TRIM(ADJUSTL(REAC(I)(IP+1:)))
        REAC(I)(IP:)=' '
       ENDIF
      ENDDO

      RETURN
      END     
      !-----------------------------------------------------------------
      SUBROUTINE ZDP_SET_REAC_SPEC()
      IMPLICIT NONE

      CALL ZDP_SET_REAC_F()
      CALL ZDP_SET_SPECIES_FROM_LIST()
      CALL ZDP_SET_NEUTRAL_IF_ANY()

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE ZDP_SET_SPECIES_FROM_LIST()
      IMPLICIT NONE 
      INTEGER :: I,J
      
      DO J=1,NREAC
       DO I=1,NSPEC
        IF(ZDP_IS_SPEC_IN_REACTION(REAC_F(J),SPEC(I))) THEN 
         RSPEC(J,I)=1
        ENDIF
       ENDDO
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE ZDP_SET_NEUTRAL_IF_ANY()
      INTEGER :: I,J
      
      DO J=1,NREAC
       IF(ZDP_IS_ANY_NEUTRAL_REACTION(REAC_F(J))) THEN
        DO I=1,NSPEC
         IF(.NOT.IS_SPEC_CHARGED(I)) THEN
          RSPEC(J,I)=1
         ENDIF
        ENDDO
       ENDIF
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE ZDP_SET_REAC_F()
      !
      ! AUTHOR: Z. NIKOLAOU  
      !  
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
       REAC_F(I)=TRIM(ADJUSTL(C))
      ENDDO
      
      RETURN
      END
      !-----------------------------------------------------------------
      FUNCTION ZDP_IS_SPEC_IN_REACTION(REAC,SPEC)
      IMPLICIT NONE
      LOGICAL :: CASEA,CASEB,CASEC,ZDP_IS_SPEC_IN_REACTION
      CHARACTER(LEN=*) :: REAC,SPEC
      CHARACTER(LEN=LEN(TRIM(ADJUSTL(SPEC)))+2) :: CA
      
      CA=' '//TRIM(ADJUSTL(SPEC))//' ' 
      IF(INDEX(REAC,CA).GT.0) THEN
       ZDP_IS_SPEC_IN_REACTION=.TRUE.
      ELSE
       ZDP_IS_SPEC_IN_REACTION=.FALSE.
      ENDIF

      END FUNCTION
      !-----------------------------------------------------------------
      SUBROUTINE ZDP_SET_IS_SPEC_CHARGED()
      IMPLICIT NONE
      INTEGER :: I

      DO I=1,NSPEC
       IS_SPEC_CHARGED(I)=(INDEX(SPEC(I),'^+').GT.0).OR. &
                          (INDEX(SPEC(I),'^-').GT.0).OR. &
                          (TRIM(ADJUSTL(SPEC(I))).EQ.'E')
      ENDDO
             
      END
      !-----------------------------------------------------------------
      FUNCTION ZDP_IS_ANY_NEUTRAL_REACTION(REAC)
      IMPLICIT NONE
      CHARACTER(LEN=*) :: REAC
      LOGICAL :: ZDP_IS_ANY_NEUTRAL_REACTION

      ZDP_IS_ANY_NEUTRAL_REACTION=INDEX(REAC,KEY_NEUTRAL).GT.0

      END FUNCTION
      !-----------------------------------------------------------------
      FUNCTION ZDP_IS_BOLSIG_REACTION(REAC)
      IMPLICIT NONE
      LOGICAL :: ZDP_IS_BOLSIG_REACTION
      CHARACTER(LEN=*) :: REAC
      
      IF(INDEX(REAC,KEY_BOLS).GT.0) THEN
       ZDP_IS_BOLSIG_REACTION=.TRUE.
      ELSE
       ZDP_IS_BOLSIG_REACTION=.FALSE.
      ENDIF

      END FUNCTION
      !-----------------------------------------------------------------
      FUNCTION ZDP_IS_GROUP_SPEC_REAC(R)
      IMPLICIT NONE
      LOGICAL :: ZDP_IS_GROUP_SPEC_REAC
      CHARACTER(LEN=*) :: R

      ZDP_IS_GROUP_SPEC_REAC=INDEX(R,KEY_AT).NE.0

      END FUNCTION
      !-----------------------------------------------------------------
      FUNCTION ZDP_ISKEYWORD(LINE)
      IMPLICIT NONE
      CHARACTER(LEN=*) :: LINE
      LOGICAL :: ZDP_ISKEYWORD

      ZDP_ISKEYWORD=INDEX(LINE,KEY_SET).GT.0

      END FUNCTION
      !-----------------------------------------------------------------
      SUBROUTINE ZDP_READ_SECTION(LSTART,LEND,NAR_MX,NSTR_MX,CARR,NAR, &
                      IS_SPLIT) 
      IMPLICIT NONE
      INTEGER :: LSTART,LEND,NAR_MX,NSTR_MX,NAR,GET_KEY_INDEX,IS,IE, &
              I,J,JOFS,NA 
      CHARACTER(LEN=NSTR_MX) :: CARR(NAR_MX)
      CHARACTER(LEN=NSMX) :: COLMS(NAR_MX)
      LOGICAL :: ISCOMMENT,ISEMPTY,IS_SPLIT

      WRITE(*,*) '***ZDP_READ_SECTION***'
      WRITE(*,'(AXI6XAI6)') 'LINES:',LSTART,'-',LEND
      JOFS=0
      NAR=0
      DO I=LSTART+1,LEND-1
       IF(.NOT.ISCOMMENT(LINES(I))) THEN
        IF(.NOT.ISEMPTY(LINES(I))) THEN
         IF(.NOT.ZDP_ISKEYWORD(LINES(I))) THEN
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
      ENDDO
      WRITE(*,*) '***ZDP_READ_SECTION***'

      RETURN
      END 
      !-----------------------------------------------------------------
      FUNCTION ZDP_GET_LINE_NO(KEY,ISTART)
      IMPLICIT NONE
      CHARACTER(LEN=*) KEY
      INTEGER :: I,ISTART,ZDP_GET_LINE_NO

      DO I=ISTART,NLINES
       IF(TRIM(ADJUSTL(LINES(I))).EQ.KEY) EXIT
      ENDDO
      ZDP_GET_LINE_NO=I

      END FUNCTION
      !-----------------------------------------------------------------     
      END MODULE
      !-----------------------------------------------------------------
