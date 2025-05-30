      !-----------------------------------------------------------------
      !
      ! AUTHOR: Z. NIKOLAOU
      !
      !-----------------------------------------------------------------
      SUBROUTINE DRIVER_DRG(RR,SET_TRG,OIC) 
      USE GLOBAL, ONLY : NSMX,NSPEC,NREAC,NTRG,SPEC,DELTANU,RSPEC, &
                         INDX_TRG,ETOL
      USE GRAPH_SEARCH
      IMPLICIT NONE
      INTEGER :: SET_TRG(NTRG,NSPEC),I,J,NO_NEIGH(NSPEC), & 
                 NEIGH_INDEX(NSPEC,NSPEC)  
      DOUBLE PRECISION :: RR(NREAC),DIC(NSPEC,NSPEC),OIC(NTRG,NSPEC)
  
      SET_TRG(1:NTRG,1:NSPEC)=0
      CALL GET_NEIGH(NEIGH_INDEX,NO_NEIGH)
      CALL GET_DIC(RR,DIC) 
      DO I=1,NTRG
       SET_TRG(I,INDX_TRG(I))=1
       CALL GRPH_DIJKSTRA(NSPEC,DIC,NEIGH_INDEX,NO_NEIGH,INDX_TRG(I), &
                          OIC(I,:))      
       DO J=1,NSPEC
        IF(OIC(I,J).GT.ETOL(I)) SET_TRG(I,J)=1
       ENDDO         
      ENDDO

      END
      !-----------------------------------------------------------------
      SUBROUTINE GET_DIC(RR,DIC)
      USE GLOBAL, ONLY : NSPEC,NREAC,DELTANU,RSPEC,SPEC,ZERO,ONE
      IMPLICIT NONE
      INTEGER :: I,J,K,N
      DOUBLE PRECISION :: RR(NREAC),DIC(NSPEC,NSPEC),FT,DTRM,PTRM,MXVL, &
                          WKI(NREAC,NSPEC)

      DIC(1:NSPEC,1:NSPEC)=ZERO

      DO I=1,NSPEC !TRG
       CALL GET_WKI(I,RR,PTRM,DTRM,WKI)
       MXVL=MAX(DTRM,PTRM)
       DO J=1,NSPEC !DEP
        CALL GET_WKIDIJ(I,J,WKI,FT)
        IF(MXVL.NE.ZERO) DIC(I,J)=FT/MXVL         
         IF(DIC(I,J).GT.ONE) THEN
          WRITE(*,*) '*ERROR:GET_DIC: DIC>1.0. TERMINATING ...'
          WRITE(*,*) TRIM(ADJUSTL(SPEC(I))),'->', &
                     TRIM(ADJUSTL(SPEC(J))), &
                     'DTRM=',DTRM,'PTRM=',PTRM,'MXVL=',MXVL, &
                     'FT=',FT,'DIC=',DIC(I,J)
          STOP
         ENDIF
       ENDDO
      ENDDO 
      
      END
      !-----------------------------------------------------------------
      SUBROUTINE GET_WKI(I,RR,PTRM,DTRM,WKI)
      USE GLOBAL, ONLY : NSPEC,NREAC,ZERO,DELTANU
      IMPLICIT NONE
      INTEGER :: I,K
      DOUBLE PRECISION :: RR(NREAC),PTRM,DTRM,WKI(NREAC,NSPEC)
       
      DTRM=ZERO
      PTRM=ZERO
      DO K=1,NREAC    
       WKI(K,I)=DELTANU(K,I)*RR(K)          
       DTRM=DTRM+MAX(-WKI(K,I),ZERO)
       PTRM=PTRM+MAX(WKI(K,I),ZERO)    
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE GET_WKIDIJ(I,J,WKI,S)
      USE GLOBAL, ONLY : NSPEC,NREAC,RSPEC,ZERO
      IMPLICIT NONE
      INTEGER :: I,J,K
      DOUBLE PRECISION :: S,WKI(NREAC,NSPEC)
       
      S=ZERO
      DO K=1,NREAC
       IF(RSPEC(K,J).EQ.1.AND.J.NE.I) S=WKI(K,I)+S                    
      ENDDO
      S=ABS(S)

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE GET_NEIGH(NEIGH_INDEX,NO_NEIGH)
      USE GLOBAL, ONLY : NSPEC,NREAC,RSPEC
      IMPLICIT NONE
      INTEGER :: I,J,K,NEIGH_INDEX(NSPEC,NSPEC),NO_NEIGH(NSPEC)
      LOGICAL :: FLAG(NSPEC,NSPEC)

      FLAG(1:NSPEC,1:NSPEC)=.FALSE.
      NO_NEIGH(1:NSPEC)=0
      DO I=1,NSPEC
       DO J=1,NSPEC
        DO K=1,NREAC
         IF(RSPEC(K,J).EQ.1.AND.J.NE.I) THEN                                       
          IF(.NOT.(FLAG(I,J))) THEN
           FLAG(I,J)=.TRUE.
           NO_NEIGH(I)=NO_NEIGH(I)+1
           NEIGH_INDEX(I,NO_NEIGH(I))=J       
          ENDIF
         ENDIF 
        ENDDO
       ENDDO
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
