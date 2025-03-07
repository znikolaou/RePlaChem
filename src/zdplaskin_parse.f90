      !-----------------------------------------------------------------
      MODULE ZDPLASKIN_PARSE
      USE GLOBAL, ONLY : ELEM,SPEC,REAC,NELEM,NSPEC,NREAC,NSFLMX, &
                         NSMX_LINE,NLINEMX
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
      
      FLCHEM=FL
      CALL READ_LINES(FL,LINES,NLINES)
      CALL ZDP_READ_SECTION_ELEMENTS()
      CALL ZDP_READ_SECTION_SPECIES()
      CALL ZDP_READ_SECTION_REACTIONS()
      CALL ZDP_READ_SECTION_BOLSIG()

      RETURN
      END 
      !-----------------------------------------------------------------        
      SUBROUTINE ZDP_READ_SECTION_ELEMENTS() 
      IMPLICIT NONE
      INTEGER :: GET_KEY_INDEX,IS,IE,I

      IS=GET_KEY_INDEX(KEY_ELEM,NLINES,LINES(1:NLINES),1)
      IE=GET_KEY_INDEX(KEY_END,NLINES,LINES(1:NLINES),IS)
      WRITE(*,*) 'ELEMENTS SECTION, LINES:',IS,'-',IE
      DO I=IS+1,IE-1
       !IF(.NOT.ISCOMMENT(LINES(I))) THEN
       ! IF(.NOT.ISEMPTY(LINES(I))) THEN
         !READ ELEMENTS LINES, SPLIT AND SAVE INTOP ARRAY.
       ! ENDIF
       !ENDIF
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
