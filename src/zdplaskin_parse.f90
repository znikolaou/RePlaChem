      !TODO: FIX 'SET' KEYWORD WHEN READING BOLSIG SECTION
      !TODO: ACCOUNT FOR REACTIONS WITH @
      !TODO: SPLIT REACTION LINE INTO REACTION PART AND RATE CONSTANT PART
      !-----------------------------------------------------------------
      MODULE ZDPLASKIN_PARSE
      USE GLOBAL, ONLY : ELEM,SPEC,SPEC_BOLSIG,REAC,NELEM,NSPEC, &
                         NSPEC_BOLSIG,NREAC,NSFLMX,NSMX_LINE,NLINEMX, & 
                         NSPMX,NREMX,NSMX,REAC_CONS
      IMPLICIT NONE
      CHARACTER(LEN=*), PARAMETER :: KEY_ELEM='ELEMENTS'
      CHARACTER(LEN=*), PARAMETER :: KEY_SPEC='SPECIES'
      CHARACTER(LEN=*), PARAMETER :: KEY_REAC='REACTIONS'
      CHARACTER(LEN=*), PARAMETER :: KEY_BOLS='BOLSIG'
      CHARACTER(LEN=*), PARAMETER :: KEY_END='END'
      CHARACTER(LEN=*), PARAMETER :: KEY_AT='@'
      CHARACTER(LEN=*), PARAMETER :: KEY_EXCL='!'
      CHARACTER(LEN=*), PARAMETER :: KEY_EQ='='


      CHARACTER(LEN=NSFLMX) :: FLCHEM
      CHARACTER(LEN=NSMX_LINE) :: LINES(NLINEMX)
      CHARACTER(LEN=NSMX) :: REAC_RAW(NREMX)
      INTEGER :: NLINES,NREAC_RAW
      
      CONTAINS
      !-----------------------------------------------------------------
      SUBROUTINE ZDP_INIT(FL)
      IMPLICIT NONE
      CHARACTER(LEN=*) FL
      
      WRITE(*,*) '***MODULE ZPDLASKING_PARSE***'

      NSPEC=0
      NREAC=0
      NELEM=0
      ELEM(1:NSPMX)=' '
      SPEC(1:NSPMX)=' '
      REAC(1:NREMX)=' '
      REAC_CONS(1:NREMX)=' '
      FLCHEM=FL

      CALL READ_LINES(FL,LINES,NLINES)
      CALL ZDP_READ_ELEMENTS()
      CALL ZDP_READ_SPECIES()
      CALL ZDP_READ_AND_FILTER_REACTIONS()
      
      WRITE(*,*) '***MODULE ZPDLASKING_PARSE***'

      RETURN
      END 
      !-----------------------------------------------------------------        
      SUBROUTINE ZDP_READ_ELEMENTS() 
      IMPLICIT NONE
      INTEGER :: GET_KEY_INDEX,I,IS,IE

      IS=GET_KEY_INDEX(KEY_ELEM,NLINES,LINES(1:NLINES),1)
      IE=GET_KEY_INDEX(KEY_END,NLINES,LINES(1:NLINES),IS)
      CALL ZDP_READ_SECTION(IS,IE,NSPMX,NSMX,ELEM,NELEM,.TRUE.)
      WRITE(*,'(AXI4)') 'NO OF ELEMENTS: ',NELEM
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
      WRITE(*,'(AXI4)') 'NO OF SPECIES: ',NSPEC
      DO I=1,NSPEC
       WRITE(*,*) I,TRIM(ADJUSTL(SPEC(I)))
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE ZDP_READ_BOLSIG() 
      IMPLICIT NONE
      INTEGER :: GET_KEY_INDEX,I,IS,IE 

      IS=GET_KEY_INDEX(KEY_BOLS,NLINES,LINES(1:NLINES),1)
      IE=GET_KEY_INDEX(KEY_END,NLINES,LINES(1:NLINES),IS)
      CALL ZDP_READ_SECTION(IS,IE,NSPMX,NSMX,SPEC_BOLSIG,NSPEC_BOLSIG, &
              .TRUE.)
      WRITE(*,'(AXI4)') 'NO OF BOLSIG SPECIES: ',NSPEC_BOLSIG
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
      CALL ZDP_READ_SECTION(IS,IE,NREMX,NSMX,REAC_RAW,NREAC_RAW,.FALSE.)
      CALL FILTER_REACTIONS()
      STOP
      WRITE(*,'(AXI4)') 'NO OF REACTIONS: ',NREAC
      DO I=1,NREAC
       WRITE(*,*) I,TRIM(ADJUSTL(REAC(I)))
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE FILTER_REACTIONS()
      IMPLICIT NONE
      INTEGER :: GET_KEY_INDEX
      INTEGER :: I,J,K,L,IREAC,ICONS,IEQ,NA,IARR(NSPMX,2),NG,IKEY,NC
      CHARACTER(LEN=NSMX) :: R,KR,CGLINE,THIRD_SPEC,ATLIST(NSPMX), &
                             GLIST(NSPMX),ROLD,RNEW,C
      CHARACTER(LEN=10) :: OLDTXT,NEWTXT
      CHARACTER(LEN=2) :: GSYMB
      CHARACTER(LEN=NSMX) :: COLMS(50)

      IREAC=0
      I=1
      L=0
      DO WHILE(I.LE.NREAC_RAW)
       !CALL EXTRACT_REAC_AND_CONS(REAC_RAW(I),R,KR)
       R=REAC_RAW(I)
       !WRITE(*,'(I4XAXAXA)') I,TRIM(ADJUSTL(R))  
       IF(IS_GROUP_SPEC_REAC(R)) THEN
        WRITE(*,*) '   ***'//TRIM(ADJUSTL(R))
        CALL GET_SUBTEXTP1_DISTINCT_LIST(R,KEY_AT,NSPMX,ATLIST,NA)
        
        !READ NEXT NA LINES FOR GROUP SPECIES
        GLIST(1:NSPMX)=' '
        DO J=1,NA
         I=I+1 
         CALL GET_GSPEC_LIST(REAC_RAW(I),NSPMX,C,NG)
         !!!Build single string for output
         GLIST(I)=TRIM(ADJUSTL(C))
        ENDDO

         !GET INDEX OF @ IN GLIST
         DO J=1,NA
          IKEY=GET_KEY_INDEX(TRIM(ADJUSTL(ATLIST(J))),NSPMX,GLIST,1)
          CALL SPLIT_STRING_WITH_SPACES(GLIST(IKEY),50,COLMS,NC)
          C=' '
          DO K=1,NG
           C=TRIM(ADJUSTL(C))//' '//TRIM(ADJUSTL(GLIST(IKEY,K)))
          ENDDO
          WRITE(*,*) TRIM(ADJUSTL(ATLIST(J))),'   :', &
                  TRIM(ADJUSTL(C))

         !!!
         
        ENDDO
        DO K=1,NG
         L=L+1
         ROLD=R
         DO J=1,NA
          IKEY=GET_KEY_INDEX(TRIM(ADJUSTL(ATLIST(J))),NSPMX,GLIST,1)
          OLDTXT=ATLIST(J)
          NEWTXT=GLIST(IKEY,K)
          CALL REPLACE_TEXT(ROLD,TRIM(ADJUSTL(OLDTXT)), &
                  TRIM(ADJUSTL(NEWTXT)),RNEW,NSMX)
          ROLD=RNEW
         ENDDO
         REAC(L)=RNEW
         WRITE(*,'(I4XA)') L,TRIM(ADJUSTL(REAC(L)))
        ENDDO
 
       ELSE
        L=L+1       
        REAC(L)=R    
        WRITE(*,'(I4XA)') L,TRIM(ADJUSTL(REAC(L)))
       ENDIF
       I=I+1
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
      SUBROUTINE GET_GSPEC_LIST(GLINE,N,GLIST,NA)
      IMPLICIT NONE
      INTEGER :: I,N,NA,IE
      CHARACTER(LEN=*) :: GLINE,GLIST(N)

      !IE=INDEX(GLINE,KEY_EQ)
      CALL SPLIT_STRING_WITH_SPACES(GLINE,N,GLIST,NA)

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
      SUBROUTINE ZDP_READ_SECTION(LSTART,LEND,NAR_MX,NSTR_MX,CARR,NAR, &
                      IS_SPLIT) 
      IMPLICIT NONE
      INTEGER :: LSTART,LEND,NAR_MX,NSTR_MX,NAR,GET_KEY_INDEX,IS,IE, &
              I,J,JOFS,NA 
      CHARACTER(LEN=NSTR_MX) :: CARR(NAR_MX)
      CHARACTER(LEN=NSMX) :: COLMS(NAR_MX)
      LOGICAL :: ISCOMMENT,ISEMPTY,IS_SPLIT

      WRITE(*,*) 'READING SECTION: LINES ',LSTART,'-',LEND
      JOFS=0
      NAR=0
      DO I=LSTART+1,LEND-1
       IF(.NOT.ISCOMMENT(LINES(I))) THEN
        IF(.NOT.ISEMPTY(LINES(I))) THEN
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
      ENDDO

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
      SUBROUTINE ZDP_READ_SECTION_SPECIES() 
      RETURN
      END
      !-----------------------------------------------------------------        
      SUBROUTINE ZDP_READ_SECTION_REACTIONS()
      RETURN
      END
      !-----------------------------------------------------------------        
      SUBROUTINE ZDP_READ_SECTION_BOLSIG()   
      RETURN
      END  
      !-----------------------------------------------------------------        
      END MODULE
      !-----------------------------------------------------------------
