      !-----------------------------------------------------------------
      SUBROUTINE READ_CHEM(DIR,CHEMFL)
      !      
      ! AUTHOR: Z. NIKOLAOU
      !
      USE ZDPLASKIN_PARSE, ONLY : ZDP_INIT
      !TODO
      !USE CHEMKIN_PARSE, ONLY : CK_INIT

      IMPLICIT NONE
      CHARACTER(LEN=*) :: DIR,CHEMFL
    
      WRITE(*,*) '***READ_CHEM***'
   
      CALL ZDP_INIT(DIR//CHEMFL)
       !SET VARS:NELEM,NSPEC,NREAC,ELEM,SPEC,REAC,REAC_SPEC,RSPEC
      
      WRITE(*,*) '***READ_CHEM***'

      RETURN
      END
      !------------------------------------------------------------------
      
