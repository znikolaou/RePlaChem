      !-----------------------------------------------------------------
      !
      ! AUTHOR: Z. NIKOLAOU
      ! 
      !-----------------------------------------------------------------
      MODULE GRAPH_SEARCH
      IMPLICIT NONE
     
      CONTAINS
      !-----------------------------------------------------------------
      SUBROUTINE GRPH_DIJKSTRA(NODES,WEIGHTS,NEIGHBS,N_NEIGHBS,TARG,PIC)
      IMPLICIT NONE
      ! PERFORMS GRAPH SEARCH USING DIJKSTRA'S ALGORITHM
      !             ->MODIFIED VERSION FOR CHEMISTRY APPLICATIONS.
      !
      ! REFERENCES: 
      !  E. Dijkstra. Numer. Math. 1 (1959) 261–271.
      !  T. Lu, C. Law. Proc. Combust. Inst. 30 (2005) 1333-1341.
      !  K. Niemeyer, C. Sung., M. Raju. Combust. Flame 157 (2010) 1760-1770.
      !  P. Desjardins, H. Pitsch. Combust. Flame 154 (2008) 67-81. 
      !
      INTEGER :: I,J,K,N,NODES,II,KK
      INTEGER :: NEIGHBS(NODES,NODES),N_NEIGHBS(NODES),TARG,IARR(NODES)
      DOUBLE PRECISION :: WEIGHTS(NODES,NODES),PIC(NODES),MXVL
      !
      IF(NODES.LE.0) THEN
       WRITE(*,*) ' '
       WRITE(*,*) 'GRAPH_SEARCH: NODES<=0 !'
       WRITE(*,*) ' '
       STOP
      ENDIF

      !INITIALISE
      PIC(1:NODES)=0.0D0
      PIC(TARG)=1.0D0 !THIS IS TARGET SPECIES
      DO I=1,NODES
       IARR(I)=I
      ENDDO
      N=NODES
      
      DO WHILE(IARR(1).NE.0)

       !FIND INDEX OF MAX VAL->MOST IMPORTANT SPECIES/NODE.            
       MXVL=0.0D0
       DO K=1,N
        IF(IARR(K).EQ.0) GOTO 1
        !WRITE(*,*) K,IARR(K),PIC(IARR(K)),MXVL
         IF(PIC(IARR(K)).GE.MXVL) THEN 
          II=IARR(K)
          MXVL=PIC(II)
          KK=K
         ENDIF
1       CONTINUE
       ENDDO!! II SET       
       !WRITE(*,*) II,MXVL
                            !--------------
       !REMOVE II FROM LIST
       IARR(KK)=0
       IF(KK.LE.(N-1)) THEN
        DO K=KK,N-1
         IARR(K)=IARR(K+1)
        ENDDO
       ENDIF 
       IARR(N)=0
                            !---------------
       N=N-1 !REDUCE SEARCH SIZE BY 1 EVERY TIME. 
                            !---------------

       !RUN THROUGH NEIGHBS OF II->GET PIC
       DO I=1,N_NEIGHBS(II)
        J=NEIGHBS(II,I)
        IF(WEIGHTS(II,J).GT.(1.0D-30)) THEN
         PIC(J)=MAX(PIC(J),PIC(II)*WEIGHTS(II,J))
        ENDIF
       ENDDO      

      ENDDO

      RETURN

      END 
      !-----------------------------------------------------------------
      END MODULE GRAPH_SEARCH
      !-----------------------------------------------------------------
