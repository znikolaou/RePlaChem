      !-----------------------------------------------------------------
      !
      ! AUTHOR: Z. NIKOLAOU
      !
      !-----------------------------------------------------------------
      PROGRAM MAIN_REDCHEM        
      USE GLOBAL
      USE CHEM_PARSE, ONLY : CM_INIT
      IMPLICIT NONE
      INTEGER :: NREAC_BOLSIG_SKEL,I,J,K,N,L
      INTEGER, ALLOCATABLE :: TRGD_SET(:,:),TRGD_USET(:,:),SP_USET(:), &
       RE_USET(:),EL_USET(:),SPBOLS_USET(:)
      DOUBLE PRECISION, ALLOCATABLE :: RR(:),OIC(:,:),STATS(:,:,:)
      CHARACTER(LEN=NSFLMX) :: FLCASE,BUILD_CASE_DIR
      
      CALL READ_CONTROL()           
      
      CALL CM_INIT(INDIR//CHEMFL)
      
      ALLOCATE(TRGD_SET(NTRG,NSPEC))
      ALLOCATE(TRGD_USET(NTRG,NSPEC))
      ALLOCATE(SP_USET(NSPEC))
      ALLOCATE(RE_USET(NREAC))
      ALLOCATE(EL_USET(NELEM))
      ALLOCATE(SPBOLS_USET(NSPEC))
      ALLOCATE(RR(NREAC))
      ALLOCATE(OIC(NTRG,NSPEC))
      ALLOCATE(STATS(NTRG,NSPEC,3))
     
      TRGD_SET(1:NTRG,1:NSPEC)=0
      TRGD_USET(1:NTRG,1:NSPEC)=0
      SP_USET(1:NSPEC)=0
      RE_USET(1:NREAC)=0
      EL_USET(1:NELEM)=1
      SPBOLS_USET(1:NSPEC)=0
      RR(1:NREAC)=ZERO
      OIC(1:NTRG,1:NSPEC)=ZERO
      STATS(:,:,:)=ZERO
      
      WRITE(*,'(A)') 'READING REACTION RATE DATA ...'       
      DO N=1,NCASE
       FLCASE=BUILD_CASE_DIR(N)        
       DO I=1,NDATA
        CALL READ_REACTION_RATES(I,FLCASE,RR)
        CALL DRIVER_DRG(RR,TRGD_SET,OIC)   
        CALL GET_SPECIES_UNION_SET_FOR_TARGET(TRGD_SET,TRGD_USET)
        CALL GET_STATS(OIC,STATS)
       ENDDO
      ENDDO
      
      CALL GET_SPECIES_SET(TRGD_USET,SP_USET) 
      CALL GET_REACTION_SET(SP_USET,RE_USET,NREAC_BOLSIG_SKEL)
      CALL GET_BOLSIG_SET(SP_USET,SPBOLS_USET) 

      CALL WRITE_SUMMARY(TRGD_USET,STATS,SP_USET,RE_USET)  
      CALL WRITE_CHEM_MECH_FMT_ZDP(EL_USET,SP_USET,SPBOLS_USET, &
                                   RE_USET,NREAC_BOLSIG_SKEL)
 
      WRITE(*,'(A)') 'RUN COMPLETED SUCCESSFULLY!'
      
      DEALLOCATE(TRGD_SET)
      DEALLOCATE(TRGD_USET)
      DEALLOCATE(SP_USET)
      DEALLOCATE(RE_USET)
      DEALLOCATE(EL_USET)
      DEALLOCATE(SPBOLS_USET)
      DEALLOCATE(RR)
      DEALLOCATE(OIC)
      DEALLOCATE(STATS)

      STOP
      END
      !-----------------------------------------------------------------
