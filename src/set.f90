      !-----------------------------------------------------------------
      !
      ! AUTHOR: Z. NIKOLAOU
      !
      !-----------------------------------------------------------------
      SUBROUTINE GET_SPECIES_UNION_SET_FOR_TARGET(SET_TRG_DEP,SET_UNION)
      USE GLOBAL, ONLY : NTRG,NSPEC
      IMPLICIT NONE
      INTEGER :: I,J,SET_TRG_DEP(NTRG,NSPEC),SET_UNION(NTRG,NSPEC)

      WRITE(*,'(A)') 'CALCULATING SPECIES UNION SET FOR TARGETS ...'
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
      SUBROUTINE GET_SPECIES_SET(SPSET_TRGUP,SPSET_UNION)
      USE GLOBAL, ONLY : INDX_TRG,NSPEC,NTRG
      IMPLICIT NONE
      INTEGER :: I,J,SPSET_UNION(NSPEC),SPSET_TRGUP(NTRG,NSPEC)
    
      WRITE(*,'(A)') 'CALCULATING FINAL SPECIES UNION SET ...'  
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
      SUBROUTINE GET_REACTION_SET(SPSET_UNION,RESET_UNION, &
                                  NREAC_BOLSIG_SKEL)
      USE GLOBAL, ONLY : NSPEC,NREAC,SPEC,RSPEC,IS_BOLSIG_REAC
      IMPLICIT NONE
      INTEGER :: I,J,NREAC_SKEL,NREAC_BOLSIG_SKEL, &
                 NREAC_NONBOLSIG_SKEL, IPROD,CHECK, &
                 SPSET_UNION(NSPEC),RESET_UNION(NREAC)

      WRITE(*,'(A)') 'CALCULATING FINAL REACTION UNION SET ...'  
      RESET_UNION(1:NREAC)=0
      DO J=1,NREAC
       IPROD=1
       DO I=1,NSPEC
        IF(RSPEC(J,I).EQ.1) IPROD=IPROD*SPSET_UNION(I)
       ENDDO
       IF(IPROD.EQ.1) RESET_UNION(J)=1
      ENDDO

      DO I=1,NSPEC
       IF(SPSET_UNION(I).EQ.1) THEN
        DO J=1,NREAC
         IF(RESET_UNION(J).EQ.1.AND.RSPEC(J,I).EQ.1) THEN 
          CHECK=1
          EXIT
         ELSE
          CHECK=0
         ENDIF
        ENDDO 
        IF(CHECK.EQ.0) THEN
         WRITE(*,'(AXAXA)') '*WARNING-GET_REACTION_SET: SKEL. SPEC', &
                            TRIM(ADJUSTL(SPEC(I))), &
                            'NOT IN ANY SKEL. REACTION, REMOVING ...'
         SPSET_UNION(I)=0         
        ENDIF
       ENDIF
      ENDDO
     
      NREAC_SKEL=SUM(RESET_UNION)
      NREAC_BOLSIG_SKEL=0
      NREAC_NONBOLSIG_SKEL=0
      DO I=1,NREAC
       IF(RESET_UNION(I).EQ.1) THEN
        IF(IS_BOLSIG_REAC(I)) THEN 
         NREAC_BOLSIG_SKEL=NREAC_BOLSIG_SKEL+1
        ELSE 
         NREAC_NONBOLSIG_SKEL=NREAC_NONBOLSIG_SKEL+1
        ENDIF
       ENDIF
      ENDDO

      IF(NREAC_BOLSIG_SKEL+NREAC_NONBOLSIG_SKEL.NE.NREAC_SKEL) THEN
       WRITE(*,'(A)') '*ERROR (GET_REACTION_SET)'
       WRITE(*,'(XAXI4XAXI4XAXI4)') &
        'NREAC-BOLSIG+NREAC-NON_BOLSIG.NE.NREAC_SKEL', &
        NREAC_BOLSIG_SKEL,'+',NREAC_NONBOLSIG_SKEL,'NE',NREAC_SKEL
       WRITE(*,'(XA)') 'TERMINATING ...'
       STOP
      ENDIF

      END
      !-----------------------------------------------------------------
      SUBROUTINE GET_BOLSIG_SET(SPSET_UNION,SPBOLSET_UNION)
      USE GLOBAL, ONLY : NSPEC,NSPEC_BOLSIG,SPEC,SPEC_BOLSIG
      IMPLICIT NONE
      INTEGER :: I,SPSET_UNION(NSPEC),SPBOLSET_UNION(NSPEC)
      LOGICAL :: IS_STRING_PRESENT

      WRITE(*,'(A)') 'CALCULATING FINAL BOLSIG SPECIES UNION SET ...'  
      DO I=1,NSPEC
       IF(SPSET_UNION(I).EQ.1) THEN
        IF(IS_STRING_PRESENT(NSPEC_BOLSIG,SPEC_BOLSIG(1:NSPEC_BOLSIG), &
                TRIM(ADJUSTL(SPEC(I))))) THEN
         SPBOLSET_UNION(I)=1
        ENDIF
       ENDIF
      ENDDO

      DO I=1,NSPEC
       IF(SPBOLSET_UNION(I).EQ.1) THEN
        IF(.NOT.IS_STRING_PRESENT(NSPEC,SPEC(1:NSPEC), &
         TRIM(ADJUSTL(SPEC(I))))) THEN     
         WRITE(*,'(A)') '*ERROR (GET_BOLSIG_SET)'
         WRITE(*,'(3(XA))') 'BOLSIG SPECIES',TRIM(ADJUSTL(SPEC(I))), &
                    'NOT IN SPECIES LIST' 
         WRITE(*,'(XA)') 'TERMINATING ...'
         STOP
        ENDIF  
       ENDIF
      ENDDO 

      END
      !-----------------------------------------------------------------
