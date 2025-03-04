      MODULE GRAPH_SEARCH
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
      !
      USE PRECIS, ONLY : DBL_P
      !
      IMPLICIT NONE
      !
      CONTAINS
      !
      SUBROUTINE GRPH_DIJKSTRA(NODES,WEIGHTS,NEIGHBS,N_NEIGHBS,TARG,PIC)
      IMPLICIT NONE
      !
      !AUTHOR: Z.M. NIKOLAOU
      !
      !DESCRIPTION: PERFORMS GRAPH SEARCH USING DIJKSTRA'S ALGORITHM
      !             ->MODIFIED VERSION FOR CHEMISTRY APPLICATIONS.
      !
      !REFERENCES: 
      ! E. Dijkstra. Numer. Math. 1 (1959) 261–271.
      ! T. Lu, C. Law. Proc. Combust. Inst. 30 (2005) 1333-1341.
      ! K. Niemeyer, C. Sung., M. Raju. Combust. Flame 157 (2010) 1760-1770.
      ! P. Desjardins, H. Pitsch. Combust. Flame 154 (2008) 67-81. 
      !
      !INPUT:
      !        NODES                  ~NO OF SPECIES 
      !        WEIGHTS(NODES,NODES)   ~WIJ I->J
      !        NEIGHBS(NODES,NODES)   ~I->J,I->K,I->L ETC.
      !        N_NEIGHBS(NODES)       ~NO OF NEIGHBS FOR EACH NODE
      !        TARG                   ~TARGET NODE INDEX
      !OUTPUT:
      !        PIC(NODES)             ~MAX PATH INTERACTION COEFF. FOR TRG NODE
      ! 
      INTEGER :: I,J,K,N,NODES,II,KK
      INTEGER :: NEIGHBS(NODES,NODES),N_NEIGHBS(NODES),TARG,IARR(NODES)
      REAL(KIND=DBL_P) :: WEIGHTS(NODES,NODES),PIC(NODES),MXVL
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

      END SUBROUTINE
      !-----------------------------------------------------------------
      !
      END MODULE GRAPH_SEARCH
