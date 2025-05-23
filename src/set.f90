      !-----------------------------------------------------------------
      !
      ! AUTHOR: Z. NIKOLAOU
      !
      !-----------------------------------------------------------------
      SUBROUTINE GET_SPECIES_UNION_SET_FOR_TARGET(SET_TRG_DEP,SET_UNION)
      USE GLOBAL, ONLY : NTRG,NSPEC
      IMPLICIT NONE
      INTEGER :: I,J,SET_TRG_DEP(NTRG,NSPEC),SET_UNION(NTRG,NSPEC)

      DO J=1,NSPEC
       DO I=1,NTRG
        IF(SET_TRG_DEP(I,J).EQ.1.AND.SET_UNION(I,J).EQ.0) THEN
         SET_UNION(I,J)=1
        ENDIF
       ENDDO
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE GET_SPECIES_SET(INDX_TRG,SPSET_TRGUP,SPSET_UNION)
      USE GLOBAL, ONLY : NSPEC,NTRG
      IMPLICIT NONE
      INTEGER :: I,J,SPSET_UNION(NSPEC),INDX_TRG(NTRG), &
                 SPSET_TRGUP(NTRG,NSPEC)
     
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

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE GET_REACTION_SET(SPSET_UNION,RESET_UNION)
      USE GLOBAL, ONLY : NSPEC,NREAC,SPEC,REAC,RSPEC
      IMPLICIT NONE
      INTEGER :: I,J,N,IPROD,CHECK,SPSET_UNION(NSPEC), &
                 RESET_UNION(NREAC)

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

      END
      !-----------------------------------------------------------------
      SUBROUTINE GET_BOLSIG_SET(SPSET_UNION,SPEC_BOLSIG,SPBOLSET_UNION)
      USE GLOBAL, ONLY : NSPEC,NSPEC_BOLSIG,SPEC
      IMPLICIT NONE
      CHARACTER(LEN=*) :: SPEC_BOLSIG(NSPEC)
      INTEGER :: I,SPSET_UNION(NSPEC),SPBOLSET_UNION(NSPEC)
      LOGICAL :: IS_STRING_PRESENT

      DO I=1,NSPEC
       IF(SPSET_UNION(I).EQ.1) THEN
        IF(IS_STRING_PRESENT(NSPEC_BOLSIG,SPEC_BOLSIG(1:NSPEC_BOLSIG), &
                TRIM(ADJUSTL(SPEC(I))))) THEN
         SPBOLSET_UNION(I)=1
        ENDIF
       ENDIF
      ENDDO 

      END
      !-----------------------------------------------------------------
