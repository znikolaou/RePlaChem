      !-----------------------------------------------------------------
      !
      ! AUTHOR: Z. NIKOLAOU
      !
      !-----------------------------------------------------------------
      SUBROUTINE DRIVER_DRG(NSPEC,NREAC,NTRG,INDX_TRG,ETOL, &
                            DNU,IDB,RR, &
                            LEN_CSP,LEN_CRE,CSPECNM,CREACNM,SET_TRG)        
      USE GRAPH_SEARCH
      IMPLICIT NONE
      INTEGER :: NSPEC,NREAC,NTRG,INDX_TRG(NTRG), &
                 LEN_CSP,LEN_CRE,IDB(NREAC,NSPEC), &
                 SETA(NSPEC),SET_TRG(NTRG,NSPEC)      
      CHARACTER(LEN=LEN_CSP) :: CSPECNM(NSPEC),SRT_SPECNM(NSPEC)
      CHARACTER(LEN=LEN_CRE) :: CREACNM(NREAC)
      DOUBLE PRECISION :: ETOL(NTRG),RR(NREAC), &
                          DNU(NREAC,NSPEC),MXVL
      INTEGER :: I,J,K,NEIGHB(NSPEC,NSPEC),N_NEIGHB(NSPEC)   
      DOUBLE PRECISION :: DIC(NSPEC,NSPEC),DIC_PATH(NTRG,NSPEC), &
                          SRT_DIC(NSPEC)
  
      WRITE(*,*)
      WRITE(*,*) '***DRIVER_DRG***'
      WRITE(*,*) 

      SET_TRG(1:NTRG,1:NSPEC)=0
      CALL GET_DIRECT_INTER_COEFF(NSPEC,NREAC,RR,DNU,IDB,CSPECNM, &
                                 CREACNM, &
                                 DIC,NEIGHB,N_NEIGHB)       
      WRITE(*,*) 'DIRECT INTER. COEFFS:'
      DO J=1,NSPEC
       WRITE(*,*) TRIM(ADJUSTL(CSPECNM(J))),N_NEIGHB(J)
        DO I=1,N_NEIGHB(J)
         WRITE(*,*) TRIM(ADJUSTL(CSPECNM(J))),'->', &
                    TRIM(ADJUSTL(CSPECNM(NEIGHB(J,I)))), &
                    DIC(J,NEIGHB(J,I))
        ENDDO
      ENDDO

      WRITE(*,*) 'SEARCHING FOR STRONGEST PATH ...'

      DO I=1,NTRG
       CALL GRPH_DIJKSTRA(NSPEC,DIC,NEIGHB,N_NEIGHB,INDX_TRG(I), &
                         DIC_PATH(I,:)) !STRONGEST PATH FOUND (MAXIMUM PRODUCT OF DICs)
          
        !TODO: WRITE AS SEPARATE ROUTINE         
        DO J=1,NSPEC
         IF(DIC_PATH(I,J).GT.ETOL(I)) THEN
          SET_TRG(I,J)=1
         ELSE
          SET_TRG(I,J)=0
         ENDIF
         SET_TRG(I,INDX_TRG(I))=1
        ENDDO         

        WRITE(*,*) 'STRONGEST PATH FOUND, TARGET: ', &
                   TRIM(ADJUSTL(CSPECNM(INDX_TRG(I)))),',', &
                   'NO OF CONNECTIONS=',SUM(SET_TRG(I,:))         
        DO J=1,NSPEC
         IF(SET_TRG(I,J).EQ.1) THEN
          WRITE(*,*) ' ->',TRIM(ADJUSTL(CSPECNM(J))), &
                           DIC_PATH(I,J)
         ENDIF
        ENDDO 
       
        !SORT OVERALL DICs FOR EACH TARGET     
        CALL SORT(NSPEC,LEN_CSP,DIC_PATH(I,:),CSPECNM,SRT_DIC, &
                  SRT_SPECNM)

        WRITE(*,*) 'SORTED OICs'
        DO J=1,NSPEC
         WRITE(*,*) TRIM(ADJUSTL(SRT_SPECNM(J))),SRT_DIC(J)
        ENDDO

      ENDDO !FOR TARGET SPECIES 
            
            !SET_TRG(NTRG,NSPEC)

      WRITE(*,*) 
      WRITE(*,*) '***DRIVER DRG***'
      WRITE(*,*)       

      END
      !-----------------------------------------------------------------
      SUBROUTINE GET_RIJ(NSPEC,NREAC,WKI,IDB,RIJ)
      IMPLICIT NONE
      INTEGER :: NSPEC,NREAC,IDB(NREAC,NSPEC),I,J,K
      DOUBLE PRECISION :: WKI(NREAC,NSPEC),RIJ(NSPEC,NSPEC),FT,FB,ETHR
      PARAMETER(ETHR=1.0E-20)

      DO J=1,NSPEC
       DO I=1,NSPEC
        RIJ(I,J)=0.0E0
        IF(I.EQ.J) RIJ(I,J)=1.0E0
       ENDDO
      ENDDO

      !CALCULATE DIRECT INTERACTION COEFF (DIC)
      DO J=1,NSPEC
       DO I=1,NSPEC
         FT=0.0E0
         FB=0.0E0 
         DO K=1,NREAC
          FT=IDB(K,J)*ABS(WKI(K,I))+FT
          FB=ABS(WKI(K,I))+FB
          FT=MAX(FT,0.0E0)
          FB=MAX(FB,0.0E0)
         ENDDO
         IF(FB.GT.ETHR) RIJ(I,J)=MAX(FT/FB,0.0E0) 
         !WRITE(*,*) J,RIJ(J,1:NSPEC)
       ENDDO
      ENDDO
      
      END
      !-----------------------------------------------------------------
      SUBROUTINE GET_DIRECT_INTER_COEFF(NSPEC,NREAC,RR,DELTANU,IDB, &
                                        CSPECNM,CREACNM, &
                                        DIC,NEIGHB,N_NEIGHB)
      IMPLICIT NONE
      INTEGER :: NSPEC,NREAC,NEIGHB(NSPEC,NSPEC),N_NEIGHB(NSPEC), &
                 IDB(NREAC,NSPEC),I,J,K,N
      CHARACTER(LEN=*) :: CSPECNM(NSPEC)
      CHARACTER(LEN=*) :: CREACNM(NREAC)
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
         !IF(ABS(DELTANU(K,J)).NE.ZERO.AND.J.NE.I) THEN !SPEC J IN REAC K
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
          WRITE(*,*) 'GET_DIRECT_INTER_COEFF:ERROR, DIC > 1.0'
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
      SUBROUTINE GET_SETA(ITRG,NSPEC,ETOL,RAB,LEN_CSP,CSPECNM,SETA)
      IMPLICIT NONE
      INTEGER :: ITRG,NSPEC,LEN_CSP,I,J,K,DIFF,SETA(NSPEC), &
                 SET_TRG(NSPEC),SET_TRGO(NSPEC)
      CHARACTER(LEN=LEN_CSP), DIMENSION(NSPEC) :: CSPECNM
      DOUBLE PRECISION :: ETOL,RAB(NSPEC,NSPEC)

      SETA(1:NSPEC)=0
      SET_TRG(1:NSPEC)=0
      SET_TRG(ITRG)=1 
      SET_TRGO(1:NSPEC)=0
      DIFF=1
      K=0

      DO WHILE(DIFF.NE.0) 

       K=K+1

       !WRITE(*,*) 'LEVEL',K,'TARGET SET SIZE=',SUM(SET_TRG)
       !WRITE(*,*) 'TARGET SET:'
       ! DO I=1,NSPEC
       ! IF(SET_TRG(I).EQ.1) WRITE(*,*) CSPECNM(I)
       !ENDDO      

       DO I=1,NSPEC !FOR EACH NEW SPEC IN TARGET SET        
        
        IF( (SET_TRG(I).EQ.1).AND.(SET_TRGO(I).EQ.0) ) THEN           
         DO J=1,NSPEC
          IF( RAB(I,J).GT.ETOL ) THEN
           SETA(J)=1
           !WRITE(*,'(3A,E12.5)') TRIM(CSPECNM(I)),'->', &
           !                      TRIM(CSPECNM(J)),RAB(I,J)
          ENDIF
         ENDDO
        ENDIF  

       ENDDO !SETA UPDATED
       SET_TRGO(1:NSPEC)=SET_TRG(1:NSPEC)
       !WRITE(*,*) 'NEW SET:'
       !DO I=1,NSPEC
       ! IF(SETA(I).EQ.1) WRITE(*,*) CSPECNM(I)
       !ENDDO              
       DIFF=ABS(SUM(SETA(:)-SET_TRG(:)))
       SET_TRG(1:NSPEC)=SETA(1:NSPEC)
       !DO I=1,NSPEC 
       ! SETAO(I)=SETA(I)-SETAO(I)
       !ENDDO
       ! IF( (SETA(I).EQ.1).AND.(SETAO(I).EQ.0) ) THEN
       !  SETAO(I)=1
       ! ELSE
       !  SETAO(I)=0
       ! ENDIF
       !ENDDO
      ENDDO 
      !WRITE(*,*) 'FINAL SET SIZE=',SUM(SETA)
      !DO I=1,NSPEC
      ! IF(SETA(I).EQ.1) WRITE(*,*) CSPECNM(I)
      !ENDDO              
      
      END
      !-----------------------------------------------------------------
      !TODO: UNUSED (REMOVE?)
      SUBROUTINE GET_SETS(NSPEC,NREAC,IDB,ISP_SET,ISP_SET_UPD, &
                          IRE_SET_UPD)      
      IMPLICIT NONE 
      INTEGER :: NSPEC,NREAC,IDB(NREAC,NSPEC),I,J,ISP_SET(NSPEC), &
                 IRE_SET(NREAC),ISP_SET_UPD(NSPEC), &
                 IRE_SET_UPD(NREAC),IWRK,NRSPEC
      
      !*SPECIES AND REACTION ARRAYS ZEROED IN MAIN FUNCTION 

      !SPECIES SET
      DO I=1,NSPEC 
       IF( (ISP_SET(I).EQ.1).AND.(ISP_SET_UPD(I).EQ.0) ) THEN
        ISP_SET_UPD(I)=1
       ENDIF
      ENDDO !SPEC IN SET

      !REACTIONS SET
      DO J=1,NREAC
       IWRK=0
       NRSPEC=SUM(IDB(J,:))
       DO I=1,NSPEC
         IF( (ISP_SET_UPD(I).EQ.1).AND.(IDB(J,I).EQ.1) ) THEN
          IWRK=IWRK+1
         ENDIF
       ENDDO
       IF(IWRK.EQ.NRSPEC) IRE_SET_UPD(J)=1
      ENDDO
         
      END
      !-----------------------------------------------------------------
