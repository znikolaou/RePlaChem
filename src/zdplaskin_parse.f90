      !TODO: FIX 'SET' KEYWORD WHEN READING BOLSIG SECTION
      !TODO: ACCOUNT FOR REACTIONS WITH @
      !TODO: SPLIT REACTION LINE INTO REACTION PART AND RATE CONSTANT PART
      !-----------------------------------------------------------------
      MODULE ZDPLASKIN_PARSE
      USE GLOBAL, ONLY : ELEM,SPEC,SPEC_BOLSIG,REAC,NELEM,NSPEC, &
                         NSPEC_BOLSIG,NREAC,NSFLMX,NSMX_LINE,NLINEMX, & 
                         NSPMX,NREMX,NSMX
      IMPLICIT NONE
      CHARACTER(LEN=*), PARAMETER :: KEY_ELEM='ELEMENTS'
      CHARACTER(LEN=*), PARAMETER :: KEY_SPEC='SPECIES'
      CHARACTER(LEN=*), PARAMETER :: KEY_REAC='REACTIONS'
      CHARACTER(LEN=*), PARAMETER :: KEY_BOLS='BOLSIG'
      CHARACTER(LEN=*), PARAMETER :: KEY_END='END'
      CHARACTER(LEN=NSFLMX) :: FLCHEM
      CHARACTER(LEN=NSMX_LINE) :: LINES(NLINEMX)
      INTEGER :: NLINES
      
      CONTAINS
      !-----------------------------------------------------------------
      SUBROUTINE ZDP_INIT(FL)
      IMPLICIT NONE
      CHARACTER(LEN=*) FL
      
      WRITE(*,*) '***MODULE ZPDLASKING_PARSE***'

      FLCHEM=FL
      CALL READ_LINES(FL,LINES,NLINES)
      CALL ZDP_READ_ELEMENTS()
      CALL ZDP_READ_SPECIES()
      CALL ZDP_READ_BOLSIG()
      CALL ZDP_READ_REACTIONS()
      
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
      WRITE(*,'(AXI4') 'NO OF ELEMENTS: ',NELEM
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
       SUBROUTINE ZDP_READ_REACTIONS() 
      IMPLICIT NONE
      INTEGER :: GET_KEY_INDEX,I,IS,IE 

      IS=GET_KEY_INDEX(KEY_REAC,NLINES,LINES(1:NLINES),1)
      IE=GET_KEY_INDEX(KEY_END,NLINES,LINES(1:NLINES),IS)
      CALL ZDP_READ_SECTION(IS,IE,NREMX,NSMX,REAC,NREAC,.FALSE.)
      WRITE(*,'(AXI4)') 'NO OF REACTIONS: ',NREAC
      DO I=1,NREAC
       WRITE(*,*) I,TRIM(ADJUSTL(REAC(I)))
      ENDDO

      RETURN
      END
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
