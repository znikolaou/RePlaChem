      !-----------------------------------------------------------------
      MODULE ZDPLASKIN_PARSE
      USE GLOBAL, ONLY : REAC,SPEC,NSPEC,NREAC
      IMPLICIT NONE
      CHARACTER(LEN=1000) :: FLCHEM

      CONTAINS
      !-----------------------------------------------------------------
      SUBROUTINE ZDP_INIT(FL)
      IMPLICIT NONE
      CHARACTER(LEN=*) FL
      
      FLCHEM=FL

      RETURN
      END 
      !-----------------------------------------------------------------        
      SUBROUTINE ZDP_READ_SECTION_ELEMENTS() 
      RETURN
      END
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
