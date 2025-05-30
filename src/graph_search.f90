      !-----------------------------------------------------------------
      !
      ! AUTHOR: Z. NIKOLAOU
      ! 
      !-----------------------------------------------------------------
      MODULE GRAPH_SEARCH
      USE GLOBAL, ONLY : ZERO,ONE
      IMPLICIT NONE
     
      CONTAINS
      !-----------------------------------------------------------------
      SUBROUTINE DIJKSTRA_GSEARCH(NODES,WEIGHTS,NEIGH,N_NEIGH,ITRG,OIC)
      IMPLICIT NONE
      ! REFERENCES: E. Dijkstra. Numer. Math. 1 (1959) 261–271.
      INTEGER :: I,J,K,N,NODES,INODE,INDX_FOR_INODE
      INTEGER :: NEIGH(NODES,NODES),N_NEIGH(NODES),ITRG,IARR(NODES)
      DOUBLE PRECISION :: WEIGHTS(NODES,NODES),OIC(NODES),MXVL
      !
      IF(NODES.LE.0) THEN
       WRITE(*,*) ' '
       WRITE(*,*) 'GRAPH_SEARCH: NODES<=0. TERMINATING ...'
       WRITE(*,*) ' '
       STOP
      ENDIF

      OIC(1:NODES)=ZERO
      OIC(ITRG)=ONE
      DO I=1,NODES
       IARR(I)=I
      ENDDO
      N=NODES
      
      DO WHILE(IARR(1).NE.0)       
       CALL GET_MAX_LOC(N,NODES,IARR,OIC,INODE,INDX_FOR_INODE,MXVL)
       CALL UPDATE_LIST(NODES,INDX_FOR_INODE,N,IARR)
       DO I=1,N_NEIGH(INODE)
        J=NEIGH(INODE,I)
        IF(WEIGHTS(INODE,J).GT.(1.0D-30)) THEN
         OIC(J)=MAX(OIC(J),OIC(INODE)*WEIGHTS(INODE,J))
        ENDIF
       ENDDO      

      ENDDO

      RETURN

      END 
      !-----------------------------------------------------------------
      SUBROUTINE GET_MAX_LOC(N,NT,IARR,Y,II,KK,MXVL)
      IMPLICIT NONE
      INTEGER :: N,NT,II,KK,K,IARR(NT)
      DOUBLE PRECISION :: MXVL,Y(NT)

      MXVL=ZERO
      DO K=1,N
       IF(IARR(K).EQ.0) GOTO 1
        IF(Y(IARR(K)).GE.MXVL) THEN 
         II=IARR(K)
         MXVL=Y(II)
         KK=K
        ENDIF
1      CONTINUE
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE UPDATE_LIST(NT,ITORM,N,IARR)
      IMPLICIT NONE
      INTEGER :: N,K,NT,ITORM,IARR(NT)

      IARR(ITORM)=0
      IF(ITORM.LE.(N-1)) THEN
       DO K=ITORM,N-1
        IARR(K)=IARR(K+1)
       ENDDO
      ENDIF 
      IARR(N)=0
      N=N-1

      RETURN
      END
      !-----------------------------------------------------------------
      END MODULE GRAPH_SEARCH
      !-----------------------------------------------------------------
