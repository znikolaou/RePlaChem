      !-----------------------------------------------------------------
      MODULE ZDPLASKIN_PARSE
      USE GLOBAL, ONLY : ELEM,SPEC,SPEC_BOLSIG,REAC,NELEM,NSPEC, &
                         NSPEC_BOLSIG,NREAC,NSFLMX,NSMX,NLINEMX, & 
                         NSPMX,NREMX,NSMX
      IMPLICIT NONE
      CHARACTER(LEN=*), PARAMETER, PRIVATE :: KEY_ELEM='ELEMENTS', &
       KEY_SPEC='SPECIES', KEY_REAC='REACTIONS',KEY_BOLS='BOLSIG', &
       KEY_END='END', KEY_SET='SET', KEY_AT='@',KEY_EXCL='!', KEY_EQ='='
      CHARACTER(LEN=NSMX), PRIVATE :: LINES(NLINEMX)
      CHARACTER(LEN=NSMX), PRIVATE :: CHEM_LINES(NREMX)
      INTEGER, PRIVATE :: NLINES,NREAC_RAW
      
      CONTAINS

      !-----------------------------------------------------------------
      SUBROUTINE ZDP_INIT(FL)
      IMPLICIT NONE
      CHARACTER(LEN=*) FL
      
      WRITE(*,*) '***MODULE ZPDLASKIN_PARSE***'

      CALL REMOVE_TABS_FROM_FILE(FL)

      NSPEC=0
      NREAC=0
      NELEM=0
      ELEM(1:NSPMX)=' '
      SPEC(1:NSPMX)=' '
      SPEC_BOLSIG(1:NSPMX)=' '
      REAC(1:NREMX)=' '
      
      CALL READ_LINES(FL,LINES,NLINES)
      CALL ZDP_READ_ELEMENTS()
      CALL ZDP_READ_SPECIES()
      CALL ZDP_READ_BOLSIG_SPECIES()
      CALL ZDP_READ_AND_FILTER_REACTIONS() 

      WRITE(*,*) 
      WRITE(*,*) 'VARS SET: ELEM, SPEC, SPEC_BOLSIG, REAC, NSPEC,' &
                 //' NREAC, NSPEC_BOLSIG'
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
      CALL FILTER_REACTIONS()
     
      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE FILTER_REACTIONS()
      IMPLICIT NONE
      INTEGER :: I,NA,NG,ILINE,IREAC
      CHARACTER(LEN=NSMX) :: R,ATLIST(NSPMX),GLIST(NSPMX,NSPMX)

      IREAC=0
      ILINE=1
      DO WHILE(ILINE.LE.NREAC_RAW)  
      !CALL EXTRACT_REAC_AND_CONS(CHEM_LINES(I),R,KR)
       R=CHEM_LINES(ILINE)
       IF(IS_GROUP_SPEC_REAC(R)) THEN
        WRITE(*,'(4XAXA)') '*GROUP REAC FOUND: '//TRIM(ADJUSTL(R))
        CALL EXTRACT_REAC_GROUP_SPEC(ILINE,NSPMX,ATLIST,GLIST,NA,NG) 
        !GENERATE NEW REAC FOR EACH GROUP 
        DO I=1,NG
         IREAC=IREAC+1
         CALL GENERATE_REAC_FOR_GROUP(R,NA,ATLIST(1:NA),GLIST(1:NA,I), &
                 REAC(IREAC))
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

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE GENERATE_REAC_FOR_GROUP(R,NA,ATSPEC,GLIST,RNEW)
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
      SUBROUTINE EXTRACT_REAC_GROUP_SPEC(IREAC,NL,ATSPEC,GLIST,NA,NG)
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
      SUBROUTINE EXTRACT_REAC_AND_CONS(RRAW,R,KR)
      IMPLICIT NONE
      INTEGER :: I
      CHARACTER(LEN=*) :: RRAW,R,KR
      
      IF(IS_CONST_GIVEN(RRAW)) THEN
       I=INDEX(RRAW,KEY_EXCL)
       KR=TRIM(ADJUSTL(RRAW(I+1:)))
       R=TRIM(ADJUSTL(RRAW(1:I-1)))
      ELSE
       R=TRIM(ADJUSTL(RRAW))
       KR=' '
      ENDIF

      RETURN
      END     
      !-----------------------------------------------------------------
      FUNCTION IS_CONST_GIVEN(R)
      IMPLICIT NONE
      LOGICAL IS_CONST_GIVEN
      CHARACTER(LEN=*) :: R

      IS_CONST_GIVEN=INDEX(R,KEY_EXCL).NE.0

      END FUNCTION
      !-----------------------------------------------------------------
      FUNCTION IS_GROUP_SPEC_REAC(R)
      IMPLICIT NONE
      LOGICAL :: IS_GROUP_SPEC_REAC
      CHARACTER(LEN=*) :: R

      IS_GROUP_SPEC_REAC=INDEX(R,KEY_AT).NE.0

      END FUNCTION
      !-----------------------------------------------------------------
      FUNCTION ISKEYWORD(LINE)
      IMPLICIT NONE
      CHARACTER(LEN=*) :: LINE
      LOGICAL :: ISKEYWORD

      ISKEYWORD=INDEX(LINE,KEY_SET).GT.0

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
         IF(.NOT.ISKEYWORD(LINES(I))) THEN
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
      FUNCTION GET_LINE_NO(KEY,ISTART)
      IMPLICIT NONE
      CHARACTER(LEN=*) KEY
      INTEGER :: I,ISTART,GET_LINE_NO

      DO I=ISTART,NLINES
       IF(TRIM(ADJUSTL(LINES(I))).EQ.KEY) EXIT
      ENDDO
      GET_LINE_NO=I

      END FUNCTION
      !-----------------------------------------------------------------     
      END MODULE
      !-----------------------------------------------------------------
