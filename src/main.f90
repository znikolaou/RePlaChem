      !-----------------------------------------------------------------
      !
      ! AUTHOR: Z. NIKOLAOU
      !
      !-----------------------------------------------------------------
      PROGRAM MAIN_REDCHEM        
      USE GLOBAL
      USE CHEM_PARSE, ONLY : CM_INIT
      IMPLICIT NONE
      LOGICAL :: IS_STRING_PRESENT
      INTEGER :: NDATA,NCASE,NTRG,CHECK,IRDMETH,NSPEC_SKEL,NREAC_SKEL, &
                 INDX_TRG(NSPMX),I,J,K,N,IPROD
      INTEGER, ALLOCATABLE :: SET_TRG(:,:),SPSET_TRG(:,:), &
                              SPSET_TRGUP(:,:),SPSET_UNION(:), &
                              RESET_UNION(:),ELEMSET_UNION(:), &
                              SPBOLSET_UNION(:)
      DOUBLE PRECISION, ALLOCATABLE :: WIJ(:,:),RR(:),JIJW(:,:)
      DOUBLE PRECISION :: ETOL(NSPMX)
      CHARACTER(LEN=NSMX) :: COMMAND
      CHARACTER(LEN=2) :: JCASE
      CHARACTER(LEN=NSFLMX) :: FLCASE,CHEMFL,SPECFL
      
      WRITE(*,*) 'MAIN:'

      CALL READ_CONTROL(NSPMX,NTRG,NCASE,NDATA,IRDMETH,INDX_TRG,ETOL, &
                        CHEMFL,SPECFL)  
      CALL CM_INIT(INDIR//CHEMFL)
      
      ALLOCATE(WIJ(NREAC,NSPEC))
      ALLOCATE(RR(NREAC))
      ALLOCATE(JIJW(NSPEC,NSPEC))
      ALLOCATE(SET_TRG(NTRG,NSPEC))
      ALLOCATE(SPSET_TRGUP(NTRG,NSPEC))
      ALLOCATE(SPSET_UNION(NSPEC))
      ALLOCATE(RESET_UNION(NREAC))
      ALLOCATE(ELEMSET_UNION(NELEM))
      ALLOCATE(SPBOLSET_UNION(NSPEC))
      
      WIJ(1:NREAC,1:NSPEC)=ZERO
      RR(1:NREAC)=ZERO
      JIJW(1:NSPEC,1:NSPEC)=ZERO
      SPSET_TRGUP(1:NTRG,1:NSPEC)=0
      SPSET_UNION(1:NSPEC)=0
      RESET_UNION(1:NREAC)=0
      DO I=1,NELEM
       ELEMSET_UNION(I)=ONE
       SPBOLSET_UNION(I)=ZERO
      ENDDO

      WRITE(*,*) 'GOING THROUGH DATASETS'       
      WRITE(*,*) 
      DO N=1,NCASE
       WRITE(*,*) 'CASE',N
       WRITE(JCASE,'(I2.2)') N
       FLCASE=INDIR//RATE_DIR//'case_'//JCASE//'/'         
       DO I=1,NDATA
        WRITE(*,*) ' -->DATASET',I
        WRITE(*,*) '-------------------'
        !TODO: DO NOT PASS RATEFL
        CALL READ_REACTION_RATES(I,FLCASE,REACTION_RATE_FL,NREAC,RR)
        !1=MY WAY
        !2=GRAPH SEARCH (PATH DEPENDENT)
        !TODO: DO NOT PASS COM VARS
        CALL DRIVER_DRG(IRDMETH,NSPEC,NREAC,NTRG,INDX_TRG(1:NTRG), &
                        ETOL(1:NTRG),DELTANU(1:NREAC,1:NSPEC), &
                        RSPEC(1:NREAC,1:NSPEC),WIJ,RR, &
                        JIJW,NSMX,NSMX,SPEC,REAC,SET_TRG)  
        !SAVE PICs
        !UPDATE_SETS
        DO K=1,NSPEC
         DO J=1,NTRG
          IF(SET_TRG(J,K).EQ.1.AND.SPSET_TRGUP(J,K).EQ.0) THEN
           SPSET_TRGUP(J,K)=1
          ENDIF
         ENDDO  
        ENDDO   
       
       ENDDO!DATA
      ENDDO!CASE
      !-----------------------------------------------------------------     
      
      WRITE(*,*) 'FINAL TARGET SPECIES SETS:'
      DO I=1,NTRG
       WRITE(*,*) I,TRIM(SPEC(INDX_TRG(I))),SUM(SPSET_TRGUP(I,:))
      ENDDO
      WRITE(*,*) '-------------------'

      !TODO: WRITE AS SEPARATE ROUTINE
      !TODO: GET BOSLIG_SPECIES UNION SET FROM REDUCED SPECIES SET
      !GET UNION SET
      SPSET_UNION(1:NSPEC)=0
      DO I=1,NTRG
       SPSET_UNION(INDX_TRG(I))=1
      ENDDO
      DO J=1,NSPEC
       DO I=1,NTRG
        IF(INDX_TRG(I).EQ.J) CYCLE 
        IF(SPSET_TRGUP(I,J).EQ.1.AND.SPSET_UNION(J).EQ.0) THEN
         SPSET_UNION(J)=1
        ENDIF
       ENDDO
      ENDDO      

      !GET PRELIMINARY REACTION SET
      WRITE(*,*) 'EXTRACTING PRELIMINARY SKEL MECH REACTIONS'
      WRITE(*,*) '------------------------------------------'
      DO J=1,NREAC
       IPROD=1
       WRITE(*,*) J,TRIM(REAC(J))
       DO I=1,NSPEC
        IF(RSPEC(J,I).EQ.1) THEN !SPEC IN R
         IPROD=IPROD*SPSET_UNION(I)
         WRITE(*,*) TRIM(ADJUSTL(SPEC(I))),SPSET_UNION(I),IPROD
        ENDIF
       ENDDO
       IF(IPROD.EQ.1) THEN
        RESET_UNION(J)=1
        WRITE(*,*) '------------------------------KEEP'
       ELSE 
        WRITE(*,*) '------------------------------REMOVE'
       ENDIF
      ENDDO

      !GOTO 1 
      !CHECK REACTIONS
      WRITE(*,*)
      WRITE(*,*) 'CHECKING PRELIMINARY SKEL MECH REACTIONS'
      WRITE(*,*)
      N=0
      DO I=1,NSPEC
       IF(SPSET_UNION(I).EQ.1) THEN !UNION SPEC
        
        DO J=1,NREAC
         IF(RESET_UNION(J).EQ.1.AND.RSPEC(J,I).EQ.1) THEN !EXISTS IN SOME REACTION
          !WRITE(*,*) TRIM(ADJUSTL(REAC(J))) 
          CHECK=1
          !WRITE(*,'(A,X,A,X,I1)') 'CHECK ',SPEC(I),CHECK
          EXIT
         ELSE
          CHECK=0
         ENDIF
        ENDDO 
        IF(CHECK.EQ.0) THEN !SPEC IN SKEL SET NOT IN ANY REACTION IN SKEL SET -> REMOVE SPEC
         N=N+1
         !WRITE(*,*) 'ERROR: REDUCED REACTION SET NOT POSSIBLE'
         !WRITE(*,*) TRIM(ADJUSTL(REAC(J))) 
         !STOP
         WRITE(*,*) 'CHECK REACTIONS-REMOVING ', &
                     TRIM(ADJUSTL(SPEC(I))), &
                     ' FROM SKEL MECH SET'         
         SPSET_UNION(I)=0         
        ENDIF
       ENDIF
      ENDDO
      IF(N.NE.0) THEN
       WRITE(*,*) 'WARNING: SOME SKEL MECH SPECIES REMOVED'
      ENDIF
      WRITE(*,*) 

1     CONTINUE

      WRITE(*,*) 'ELIMINATED SPECIES:'
      DO I=1,NSPEC
       IF(SPSET_UNION(I).EQ.0) THEN
        WRITE(*,*) TRIM(SPEC(I)) 
       ENDIF
      ENDDO
      WRITE(*,*) '-------------------'

      WRITE(*,*) 'SPECIES UNION SET:'
      DO I=1,NSPEC
       IF(SPSET_UNION(I).EQ.1) THEN 
        WRITE(*,*) TRIM(SPEC(I))
       ENDIF
      ENDDO
      WRITE(*,*) '-------------------'

      WRITE(*,*) 'REACTIONS UNION SET:'
      DO I=1,NREAC
       IF(RESET_UNION(I).EQ.1) THEN   
        WRITE(*,*) TRIM(REAC(I))
       ENDIF
      ENDDO
      WRITE(*,*) '-------------------'
      
      DO I=1,NSPEC
       IF(SPSET_UNION(I).EQ.1) THEN
        IF(IS_STRING_PRESENT(NSPEC_BOLSIG,SPEC_BOLSIG, &
                TRIM(ADJUSTL(SPEC(I))))) THEN
         SPBOLSET_UNION(I)=1.0E0
        ENDIF
       ENDIF
      ENDDO 

      NSPEC_SKEL=SUM(SPSET_UNION)
      NREAC_SKEL=SUM(RESET_UNION)
      IF(NSPEC_SKEL.LE.0.OR.NREAC_SKEL.LE.0) THEN
       WRITE(*,*) 'ERROR: NO SKEL MECH CREATED. TERMINATING ...'
       STOP
      ENDIF
       
      CALL WRITE_CHEM_MECH_FMT_ZDP(OUTDIR,CHEMRED_FL, &
                                   NELEM,NSPEC,NREAC, &
                                   ELEM,SPEC,REAC,REAC_CONST, &
                                   ELEMSET_UNION,SPSET_UNION, &
                                   SPBOLSET_UNION,RESET_UNION)

      WRITE(*,*) 'SKELETAL MECHANISM SIZE:'
      WRITE(*,*) 'NUNION SPEC=',SUM(SPSET_UNION)
      WRITE(*,*) 'NUNION REAC=',SUM(RESET_UNION)
      WRITE(*,*) 
      WRITE(*,*) 'RUN COMPLETED SUCCESSFULLY!'
      WRITE(*,*) '-------------------'

      DEALLOCATE(WIJ)
      DEALLOCATE(RR)
      DEALLOCATE(JIJW)
      DEALLOCATE(SET_TRG)
      DEALLOCATE(SPSET_TRGUP)
      DEALLOCATE(SPSET_UNION)
      DEALLOCATE(RESET_UNION)
      DEALLOCATE(ELEMSET_UNION)
      DEALLOCATE(SPBOLSET_UNION)

      STOP
      END
