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
      CALL GET_DIC(RR,DIC,NEIGH_INDEX,NO_NEIGH) 

      DO I=1,NTRG
       SET_TRG(I,INDX_TRG(I))=1
       CALL GRPH_DIJKSTRA(NSPEC,DIC,NEIGH_INDEX,NO_NEIGH,INDX_TRG(I), &
                          OIC(I,:))      
       DO J=1,NSPEC
        IF(OIC(I,J).GT.ETOL(I)) THEN
         SET_TRG(I,J)=1
        ENDIF
       ENDDO         
      ENDDO

      END
      !-----------------------------------------------------------------
      SUBROUTINE GET_DIC(RR,DIC,NEIGH_INDEX,NO_NEIGH)
      USE GLOBAL, ONLY : NSPEC,NREAC,DELTANU,RSPEC,SPEC,ZERO,ONE
      IMPLICIT NONE
      INTEGER :: NEIGH_INDEX(NSPEC,NSPEC),NO_NEIGH(NSPEC),I,J,K,N
      DOUBLE PRECISION :: RR(NREAC),DIC(NSPEC,NSPEC),FT,DTRM,PTRM,MXVL, &
                          WKI(NREAC,NSPEC)
      LOGICAL :: RIJ_FLAG(NSPEC,NSPEC)

      RIJ_FLAG(1:NSPEC,1:NSPEC)=.FALSE.
      NO_NEIGH(1:NSPEC)=0
      NEIGH_INDEX(1:NSPEC,1:NSPEC)=0
      DIC(1:NSPEC,1:NSPEC)=ZERO

      DO I=1,NSPEC !TRG
       CALL GET_SPEC_RATE(I,RR,PTRM,DTRM,WKI)
       MXVL=MAX(DTRM,PTRM)
       DO J=1,NSPEC !DEP
        FT=ZERO
        DO K=1,NREAC
         IF(RSPEC(K,J).EQ.1.AND.J.NE.I) THEN                   
          FT=WKI(K,I)+FT                    
          IF(.NOT.(RIJ_FLAG(I,J))) THEN
           RIJ_FLAG(I,J)=.TRUE.
           NO_NEIGH(I)=NO_NEIGH(I)+1
           NEIGH_INDEX(I,NO_NEIGH(I))=J       
          ENDIF
         ENDIF 
        ENDDO
        FT=ABS(FT)
        IF(MXVL.NE.ZERO) THEN 
         DIC(I,J)=FT/MXVL         
         IF(DIC(I,J).GT.ONE) THEN
          WRITE(*,*) '*ERROR:GET_DIC: DIC > 1.0, TERMINATING ...!'
          WRITE(*,*) TRIM(ADJUSTL(SPEC(I))),' ', &
                      ' ',TRIM(ADJUSTL(SPEC(J))), &
                     DTRM,PTRM,MXVL,DTRM+PTRM-FT,DIC(I,J)

          STOP
         ENDIF
        ENDIF
       ENDDO !DEP
      ENDDO !TRG 
      
      END
      !-----------------------------------------------------------------
      SUBROUTINE GET_SPEC_RATE(I,RR,PTRM,DTRM,WKI)
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
