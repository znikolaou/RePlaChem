      !-----------------------------------------------------------------
      SUBROUTINE READ_CHEM(DIR,CHEMFL)
      !      
      ! AUTHOR: Z. NIKOLAOU
      !
      USE CHEM_PARSE, ONLY : CM_INIT

      IMPLICIT NONE
      CHARACTER(LEN=*) :: DIR,CHEMFL
    
      WRITE(*,*) '***READ_CHEM***'
   
      CALL CM_INIT(DIR//CHEMFL)
       !SET VARS:NELEM,NSPEC,NREAC,ELEM,SPEC,REAC,REAC_SPEC,RSPEC
      
      WRITE(*,*) '***READ_CHEM***'

      RETURN
      END
      !------------------------------------------------------------------
      
