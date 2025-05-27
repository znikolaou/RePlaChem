      !-----------------------------------------------------------------
      !
      ! AUTHOR: Z. NIKOLAOU
      !
      !-----------------------------------------------------------------
      SUBROUTINE DRIVER_DRG(RR,SET_TRG,DIC_PATH) 
      USE GLOBAL, ONLY : NSMX,NSPEC,NREAC,NTRG,SPEC,DELTANU,RSPEC, &
                         INDX_TRG,ETOL
      USE GRAPH_SEARCH
      IMPLICIT NONE
      INTEGER :: SETA(NSPEC),SET_TRG(NTRG,NSPEC), &      
                 I,J,K,NEIGHB(NSPEC,NSPEC),N_NEIGHB(NSPEC)   
      CHARACTER(LEN=NSMX) :: SRT_SPECNM(NSPEC)
      DOUBLE PRECISION :: RR(NREAC),MXVL,DIC(NSPEC,NSPEC), &
                          DIC_PATH(NTRG,NSPEC),SRT_DIC(NSPEC)
  
      SET_TRG(1:NTRG,1:NSPEC)=0
      CALL GET_DIRECT_INTER_COEFF(NSPEC,NREAC,RR, &
       DELTANU(1:NREAC,1:NSPEC),RSPEC(1:NREAC,1:NSPEC),SPEC(1:NSPEC), &
       DIC,NEIGHB,N_NEIGHB)       
      
      DO I=1,NTRG
       SET_TRG(I,INDX_TRG(I))=1
       SET_TRG(I,1:NSPEC)=0
       CALL GRPH_DIJKSTRA(NSPEC,DIC,NEIGHB,N_NEIGHB,INDX_TRG(I), &
                         DIC_PATH(I,:))      
       DO J=1,NSPEC
        IF(DIC_PATH(I,J).GT.ETOL(I)) THEN
         SET_TRG(I,J)=1
        ENDIF
       ENDDO         

      ENDDO !FOR TARGET SPECIES 
            
      END
      !-----------------------------------------------------------------
      SUBROUTINE GET_DIRECT_INTER_COEFF(NSPEC,NREAC,RR,DELTANU,IDB, &
                                        CSPECNM,DIC,NEIGHB,N_NEIGHB)
      IMPLICIT NONE
      INTEGER :: NSPEC,NREAC,NEIGHB(NSPEC,NSPEC),N_NEIGHB(NSPEC), &
                 IDB(NREAC,NSPEC),I,J,K,N
      CHARACTER(LEN=*) :: CSPECNM(NSPEC)
      DOUBLE PRECISION :: DELTANU(NREAC,NSPEC),RR(NREAC), &
                          DIC(NSPEC,NSPEC),FT,DTRM,PTRM,MXVL, &
                          WKI(NREAC,NSPEC),FTT
      DOUBLE PRECISION, PARAMETER :: ZERO=0.0E0,ONE=1.0E0
      LOGICAL :: RIJ_FLAG(NSPEC,NSPEC)

      RIJ_FLAG(1:NSPEC,1:NSPEC)=.FALSE.
      N_NEIGHB(1:NSPEC)=0
      NEIGHB(1:NSPEC,1:NSPEC)=0
      DIC(1:NSPEC,1:NSPEC)=ZERO

      !TRG SPEC
      DO I=1,NSPEC
       DTRM=ZERO
       PTRM=ZERO
       DO K=1,NREAC    
        WKI(K,I)=DELTANU(K,I)*RR(K)          
        DTRM=DTRM+MAX(-WKI(K,I),ZERO)
        PTRM=PTRM+MAX(WKI(K,I),ZERO)                 
       ENDDO
       MXVL=MAX(DTRM,PTRM)
       !LOOP THROUGH DEP SPEC
       DO J=1,NSPEC
        FT=ZERO
        DO K=1,NREAC
         IF(IDB(K,J).EQ.1.AND.J.NE.I) THEN                   
          FT=WKI(K,I)+FT                    
          IF(.NOT.(RIJ_FLAG(I,J))) THEN
           RIJ_FLAG(I,J)=.TRUE.
           N_NEIGHB(I)=N_NEIGHB(I)+1
           NEIGHB(I,N_NEIGHB(I))=J       
          ENDIF
         ENDIF 
        ENDDO
        FT=ABS(FT)
        IF(MXVL.NE.ZERO) THEN 
         DIC(I,J)=FT/MXVL         
         IF(DIC(I,J).GT.ONE) THEN
          WRITE(*,*) 'GET_DIRECT_INTER_COEFF:ERROR, DIC > 1.0!'
          WRITE(*,*) TRIM(ADJUSTL(CSPECNM(I))),' ', &
                      ' ',TRIM(ADJUSTL(CSPECNM(J))), &
                       DTRM,PTRM,MXVL,DTRM+PTRM-FT,DIC(I,J)

          STOP
         ENDIF
        ENDIF
       ENDDO !FOR DEP

      ENDDO !FOR TRG 
      
      END
      !-----------------------------------------------------------------
