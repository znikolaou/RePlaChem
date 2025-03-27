      SUBROUTINE DRIVER_DRG(ISEARCH,NSPEC,NREAC,NTRG,INDX_TRG,ETOL, &
                            DNU,IDB,WKJ,RR,JIJW, &
                            LEN_CSP,LEN_CRE,CSPECNM,CREACNM,SET_TRG)        
!
!     AUTHOR: ZACHARIAS M. NIKOLAOU
! 
!     DESCRIPTION: DRIVER SUBROUTINE FOR DRG. 
!
!     INPUT:
!     ISEARCH          ~SEARCH METHOD
!     NSPEC            ~NUMBER OF SPECIES
!     NREAC            ~NUMBER OF REACTIONS
!     NTRG             ~NUMBER OF TARGET SPECIES
!     INDX_TRG(NSPEC)  ~TARGET SPECIES INDICES
!     ETOL(NTRG)       ~THRESHOLD DRG ERROR
!     DNU(NSPEC,NREAC) ~DIFFERENCE IN STOICH. COEFFS. (NUP-NUR)
!     IDB(NREAC,NSPEC) ~1 IF SPEC IN REACTION, O OTHERWISE (NOT USED)
!     WKJ              ~SPECIES J RATE FROM REACTION K
!     RR(NREAC)        ~REACTION RATE
!     JIJW(NSPEC,NSPEC)~JACOBIAN DEFINED DIC
!     LEN_CSP          ~LENGTH OF SPECIES STRING NAMES
!     LEN_CRE          ~LENGTH OF REACTION STRING NAMES  
!     CSPECNM(NSPEC)   ~SPECIES NAMES
!     CREACNM(NREAC)   ~REACTION NAMES
!
!     OUTPUT:
!     SET_TRG(NSPEC)   ~TARGET SPECIES SET
! 
!-----------------------------------------------------------------------
!
!     This file is part of <REDCHEM_v0.0>     
!     Copyright (C) <2018>  <Zacharias M. Nikolaou>
!
!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program.  If not, see <https://www.gnu.org/licenses/>.
!
!     Contact details: ZachariasMNic@gmail.com
!
!-----------------------------------------------------------------------
      USE PRECIS, ONLY: DBL_P
      USE GRAPH_SEARCH
      IMPLICIT NONE
!
      INTEGER :: ISEARCH,NSPEC,NREAC,NTRG,INDX_TRG(NTRG), &
                 LEN_CSP,LEN_CRE,IDB(NREAC,NSPEC), &
                 SETA(NSPEC),SET_TRG(NTRG,NSPEC)      
      CHARACTER(LEN=LEN_CSP) :: CSPECNM(NSPEC),SRT_SPECNM(NSPEC)
      CHARACTER(LEN=LEN_CRE) :: CREACNM(NREAC)
      REAL(KIND=DBL_P) :: ETOL(NTRG),WKJ(NREAC,NSPEC),RR(NREAC), &
                          DNU(NSPEC,NREAC),JIJW(NSPEC,NSPEC),MXVL
      !
      INTEGER :: I,J,K,NEIGHB(NSPEC,NSPEC),N_NEIGHB(NSPEC)   
      REAL(KIND=DBL_P) :: DIC(NSPEC,NSPEC),DIC_PATH(NTRG,NSPEC), &
                          SRT_DIC(NSPEC)
      !
!
      WRITE(*,*) 'DRIVER_DRG'
      WRITE(*,*) '----------'
      WRITE(*,*) 'TARGET SPECIES INDEX, NAME:'
      DO I=1,NTRG
       WRITE(*,*) I,INDX_TRG(I),TRIM(ADJUSTL(CSPECNM(INDX_TRG(I))))
      ENDDO
      
      SET_TRG(:,:)=0


      IF(ISEARCH.EQ.1) THEN

       !DIC-DRG
       CALL GET_RIJ(NSPEC,NREAC,WKJ,IDB,DIC)

      ELSEIF(ISEARCH.EQ.2) THEN

       !DIC-DRGEP 
       WRITE(*,*) 'TEST FOR ISEARCH 2'
       CALL GET_DICEP(NSPEC,NREAC,RR,DNU,WKJ,CSPECNM,CREACNM, &
                      DIC,NEIGHB,N_NEIGHB)   
              
       WRITE(*,*) 'DICs:'
       DO J=1,NSPEC
        WRITE(*,*) TRIM(ADJUSTL(CSPECNM(J))),N_NEIGHB(J)
        DO I=1,N_NEIGHB(J)
         WRITE(*,*) TRIM(ADJUSTL(CSPECNM(J))),'->', &
                    TRIM(ADJUSTL(CSPECNM(NEIGHB(J,I)))), &
                    DIC(J,NEIGHB(J,I))
        ENDDO
       ENDDO

      ELSEIF(ISEARCH.EQ.3) THEN !JACOBIAN-BASED DICs
        
       CALL GET_DICEP(NSPEC,NREAC,RR,DNU,WKJ,CSPECNM,CREACNM, &
                      DIC,NEIGHB,N_NEIGHB) !ONLY NEIGH,N_NEIGH USED      

       !OVER-WRITE DICs
       DIC(:,:)=0.0D0
       DO I=1,NSPEC
        MXVL=MAXVAL(ABS(JIJW(I,:)))
        IF(MXVL.NE.0.0) THEN
         DIC(I,:)=ABS(JIJW(I,:))/MXVL
        ELSE
         DIC(I,:)=0.0D0
        ENDIF         
       ENDDO

       WRITE(*,*) 'DICs:'
       DO J=1,NSPEC
        WRITE(*,*) TRIM(ADJUSTL(CSPECNM(J))),N_NEIGHB(J)
        DO I=1,N_NEIGHB(J)
         WRITE(*,*) TRIM(ADJUSTL(CSPECNM(J))),'->', &
                    TRIM(ADJUSTL(CSPECNM(NEIGHB(J,I)))), &
                    DIC(J,NEIGHB(J,I))
        ENDDO
       ENDDO

      ENDIF

                                                    !DIC(NSPEC,NSPEC)
                                                    !NEIGHB(NSPEC,NSPEC)
                                                    !N_NEIGHB(NSPEC)
      !-----------------------------------------------------------------
 
      !GET TARGET SPEC SET
      IF(ISEARCH.EQ.1) THEN !MY WAY

       WRITE(*,*) 'SEARCHING...MY WAY'
       DO I=1,NTRG
        CALL GET_SETA(INDX_TRG(I),NSPEC,ETOL(I),DIC,LEN_CSP, &
                      CSPECNM,SETA) 
        SET_TRG(I,:)=SETA(:)
        SET_TRG(I,INDX_TRG(I))=1
       ENDDO

      ELSEIF(ISEARCH.EQ.2.OR.ISEARCH.EQ.3) THEN !GRAPH SEARCH

       WRITE(*,*) 'SEARCHING...GRAPH SEARCH'            
       DO I=1,NTRG
        CALL GRPH_DIJKSTRA(NSPEC,DIC,NEIGHB,N_NEIGHB,INDX_TRG(I), &
                         DIC_PATH(I,:)) !STRONGEST PATH FOUND (MAXIMUM PRODUCT OF DICs)
          
        DO J=1,NSPEC
         IF(DIC_PATH(I,J).GT.ETOL(I)) THEN
          SET_TRG(I,J)=1
         ELSE
          SET_TRG(I,J)=0
         ENDIF
         SET_TRG(I,INDX_TRG(I))=1
        ENDDO         
        WRITE(*,*) 'TARGET=',TRIM(ADJUSTL(CSPECNM(INDX_TRG(I)))),',', &
                   'NO OF CONNECTIONS=',SUM(SET_TRG(I,:))         
        DO J=1,NSPEC
         IF(SET_TRG(I,J).EQ.1) THEN
          WRITE(*,*) '->',TRIM(ADJUSTL(CSPECNM(J))), &
                           DIC_PATH(I,J)
         ENDIF
        ENDDO 
        WRITE(*,*)
        ! 
        ! 
        !SORT OVERALL DICs FOR EACH TARGET
        
        CALL SORT(NSPEC,LEN_CSP,DIC_PATH(I,:),CSPECNM,SRT_DIC, &
                  SRT_SPECNM)
        WRITE(*,*) 'SORTED OICs'
        DO J=1,NSPEC
         WRITE(*,'(A,G12.5)') TRIM(ADJUSTL(SRT_SPECNM(J))), &
                                SRT_DIC(J)
        ENDDO

        !DO J=1,N_NEIGHB(INDX_TRG(I))
        ! WRITE(*,*) J,'->', &
        !              TRIM(ADJUSTL(CSPECNM(NEIGHB(INDX_TRG(I),J)))), &
        !              DIC_PATH(I,NEIGHB(INDX_TRG(I),J)), &
        !              DIC(INDX_TRG(I),NEIGHB(INDX_TRG(I),J)), &
        !              SET_TRG(I,NEIGHB(INDX_TRG(I),J))
        !ENDDO
        
       ENDDO 
       
      ENDIF
                                                    !SET_TRG(NTRG,NSPEC)
      !-----------------------------------------------------------------

      !OUTPUT SETS
      !DO I=1,NTRG             
      ! WRITE(*,'(A,X,A,X,I5)') 'SET SIZE FOR TARGET ', &
      !             TRIM(ADJUSTL(CSPECNM(INDX_TRG(I)))),SUM(SET_TRG(I,:))
      ! DO J=1,NSPEC
      !  IF(SET_TRG(I,J).EQ.1) THEN
      !   WRITE(*,*) TRIM(ADJUSTL(CSPECNM(J))),'->',DIC_PATH(I,J)
      !  ENDIF
      ! ENDDO              
      ! WRITE(*,*) '--------------------------'
      !ENDDO
 
      WRITE(*,*) 'EXITING DRIVER DRG'
      WRITE(*,*) '------------------'
            
      END SUBROUTINE
!-----------------------------------------------------------------------
!
      SUBROUTINE GET_RIJ(NSPEC,NREAC,WKI,IDB,RIJ)
!
!     AUTHOR: ZACHARIAS M. NIKOLAOU
!
!     DESCRIPTION:  
!     CALCULATES ALL SPECIES INTER-DEPENDENCIES RIJ
!     RIJ=SUM_K |DELTA_{JK}*W_{KI}| / SUM_K |W_KI|
!
!     INPUT: 
!     NSPEC           ~TOTAL NUMBER OF SPECIES
!     NREAC           ~TOTAL NUMBER OF REACTIONS
!     WKI(NREAC,NSPEC) ~SPECIES I RATE IN REACTION K
!     IDB(NREAC,NSPEC)    ~INDEX:
!                         ~1=SPECIES IN REACTION K
!                         ~0=SPECIES NOT IN REACTION K
!
!     OUTPUT:            
!     RIJ(NSPEC,NSPEC)
! 
      USE PRECIS, ONLY: DBL_P
!
      IMPLICIT NONE
      INTEGER :: NSPEC,NREAC,IDB(NREAC,NSPEC)
      REAL(KIND=DBL_P) :: WKI(NREAC,NSPEC)
      !
      INTEGER :: I,J,K
      REAL(KIND=DBL_P) :: RIJ(NSPEC,NSPEC),FT,FB,ETHR
      PARAMETER(ETHR=1.0E-20)
!
!
      !INITIALISE
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
!      
      END SUBROUTINE
!-----------------------------------------------------------------------
!
      SUBROUTINE GET_DICEP(NSPEC,NREAC,RR,DELTANU,WIKMATRIX,CSPECNM, &
                           CREACNM,DIC,NEIGHB,N_NEIGHB)
!
!     AUTHOR: ZACHARIAS M. NIKOLAOU
!
!     DESCRIPTION:  
!     CALCULATES ALL SPECIES INTER-DEPENDENCIES RIJ
!     DIJ=SUM_K |DELTA_{JK}*W_{KI}| / MAX(|PI|,|DI|)
!
!     INPUT: 
!     NSPEC                ~TOTAL NUMBER OF SPECIES
!     NREAC                ~TOTAL NUMBER OF REACTIONS
!     RR(NREAC)            ~REACTION RATE
!     DELTANU(NREAC,NSPEC) 
!
!     OUTPUT:            
!     DIC(NSPEC,NSPEC)    ~DIC
!     NEIGHB(NSPEC,NSPEC) ~NEIGHBOURS LIST
!     N_NEIGHB(NSPEC)     ~NO OF NEIGHBOURS 
! 
      USE PRECIS, ONLY: DBL_P
!
      IMPLICIT NONE
      INTEGER :: NSPEC,NREAC,NEIGHB(NSPEC,NSPEC),N_NEIGHB(NSPEC)
      CHARACTER(LEN=*) :: CSPECNM(NSPEC)
      CHARACTER(LEN=*) :: CREACNM(NREAC)
      REAL(KIND=DBL_P) :: DELTANU(NSPEC,NREAC),RR(NREAC), &
                          WIKMATRIX(NREAC,NSPEC),DIC(NSPEC,NSPEC)

      !
      LOGICAL :: RIJ_FLAG(NSPEC,NSPEC)
      INTEGER :: I,J,K,N
      REAL(KIND=DBL_P) :: FT,DTRM,PTRM,MXVL,WIK,FTT

      !INITIALISE
      RIJ_FLAG(:,:)=.FALSE.
      DO J=1,NSPEC
       N_NEIGHB(J)=0
       DO I=1,NSPEC
        NEIGHB(I,J)=0
        DIC(I,J)=0.0E0
       ENDDO
      ENDDO

      !CALCULATE DIRECT INTERACTION COEFF (DIC_IJ)
      DO J=1,NSPEC
       DO I=1,NSPEC
         DTRM=0.0E0
         PTRM=0.0E0
         FT=0.0E0 
         FTT=0.0
         DO K=1,NREAC                                 
          IF( ABS(DELTANU(I,K)).NE.0.0 ) THEN !SPEC I IN REAC K

           !WIK=WIKMATRIX(K,I)
           WIK=DELTANU(I,K)*RR(K)
           
           DTRM=DTRM+MAX(-WIK,0.0E0)
           PTRM=PTRM+MAX(WIK,0.0E0)

            !FTT=ABS(WIK)+FTT                    

           IF( ABS(DELTANU(J,K)).NE.0.0.AND.J.NE.I ) THEN !SPEC J IN REAC K 
 
            IF( .NOT.(RIJ_FLAG(I,J)) ) THEN
             RIJ_FLAG(I,J)=.TRUE.
             
             N_NEIGHB(I)=N_NEIGHB(I)+1
             NEIGHB(I,N_NEIGHB(I))=J       
            ENDIF

            FT=WIK+FT                    

           ENDIF !SPEC J IN REAC K 

          ENDIF !SPEC I IN REAC K
                  
         ENDDO!REAC LOOP

         FT=ABS(FT)

         !CHEKSUM
         !IF(PTRM+DTRM-FTT.GT.1.0D-6) THEN
         ! WRITE(*,*) 'GET_DICEP:ERROR, CHECKSUM'
         ! WRITE(*,*) TRIM(ADJUSTL(CSPECNM(I))),' ', &
         !             ' ',TRIM(ADJUSTL(CSPECNM(J))), &
         !             ' ','CHECKSUM ',(PTRM+DTRM-FTT),FTT,PTRM+DTRM
         ! STOP
         !ENDIF

         !CALCULATE DIC
         MXVL=MAX(DTRM,PTRM)

         IF(MXVL.LT.0.0) THEN
          WRITE(*,*) 'GET_DICEP:ERROR, MXVL.LT.0.0'
          STOP
         ENDIF

         IF(MXVL.NE.0.0) THEN 
          DIC(I,J)=FT/MXVL         
          IF(DIC(I,J).GT.1.0D0) THEN
           WRITE(*,*) 'GET_DICEP:ERROR, DIC > 1.0'
           WRITE(*,*) TRIM(ADJUSTL(CSPECNM(I))),' ', &
                      ' ',TRIM(ADJUSTL(CSPECNM(J))), &
                       DTRM,PTRM,MXVL,DTRM+PTRM-FT,DIC(I,J)

           STOP
          ENDIF
         ENDIF

       ENDDO
      ENDDO
!      
      END SUBROUTINE
!-----------------------------------------------------------------------
!
      SUBROUTINE GET_SETA(ITRG,NSPEC,ETOL,RAB,LEN_CSP,CSPECNM,SETA)
!
!     AUTHOR: ZACHARIAS M. NIKOLAOU
!       
!     DESCRIPTION:
!     SEARCHES THROUGH DEPENENCY MATRIX TO FORM THE DRG
!
!     INPUT:
!     ITRG  ~TARGET SPECIES INDEX
!     NSPEC ~NO OF SPECIES
!     ETOL  ~DRG ERROR
!     RAB   ~DEPENDENCE INDICATOR A -> B
!     LEN_CSP ~LENGTH OF SPECIES STRING
!     CSPECNM ~SPECIES NAMES
!     OUTPUT:
!
      USE PRECIS, ONLY: DBL_P
!
      IMPLICIT NONE
!
      INTEGER :: ITRG,NSPEC,LEN_CSP
      CHARACTER(LEN=LEN_CSP), DIMENSION(NSPEC) :: CSPECNM
      REAL(KIND=DBL_P) :: ETOL,RAB(NSPEC,NSPEC)
! 
      INTEGER :: I,J,K,DIFF,SETA(NSPEC),SET_TRG(NSPEC),SET_TRGO(NSPEC)
!
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
      
      END SUBROUTINE
!-----------------------------------------------------------------------
!
      SUBROUTINE GET_SETS(NSPEC,NREAC,IDB,ISP_SET,ISP_SET_UPD, &
                          IRE_SET_UPD)      
!
!     AUTHOR: ZACHARIAS M. NIKOLAOU
!
!     DESCRIPTION: 
!     UPDATES SPECIES AND REACTION SETS GIVEN THE SETS OF IMPORTANT
!     SPECIES AND REACTIONS FOR EVERY SAMPLE POINT.
!
!     INPUT: 
!     NSPEC ~NO OF SPECIES
!     NREAC ~NO OF REACTIONS
!     IDB(NSPEC,NSPEC) ~ 1 IF SPEC IN REAC, 0 OTHERWISE
!     ISP_SET(NSPEC) ~DEPEDENT SPEC SET
!
      IMPLICIT NONE
! 
      INTEGER :: NSPEC,NREAC,IDB(NREAC,NSPEC)
      INTEGER :: I,J,ISP_SET(NSPEC),IRE_SET(NREAC)
      INTEGER :: ISP_SET_UPD(NSPEC),IRE_SET_UPD(NREAC),IWRK,NRSPEC
      
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
         
!
      END SUBROUTINE
