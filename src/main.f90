      PROGRAM MAIN_REDCHEM
!
!     <REDCHEM_v0.0>
! 
!     AUTHOR: Z.M. NIKOLAOU
!
!     DESCRIPTION: PRODUCES SKELETAL CHEMISTRY FOR KPP-STYLE
!                  CHEMICAL MECHANISMS.
!
!     INPUT: AS BELOW.
!
!-----------------------------------------------------------------------
!     
!     <REDCHEM_v0.0>
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
          
      USE GLOBAL, ONLY: NREMX,NSMX,NSPMX,NREMX 
      USE PRECIS, ONLY: DBL_P
      IMPLICIT NONE
      INTEGER :: NDATA,NCASE,NTRG,CHECK,IRDMETH,NSPEC_SKEL,NREAC_SKEL
      INTEGER, ALLOCATABLE :: INDX_TRG(:)
      REAL(KIND=DBL_P), ALLOCATABLE :: ETOL(:)
      CHARACTER(LEN=*) :: INDIR,OUTDIR,RATEDIR
      CHARACTER(LEN=*) :: TRGFILE,CHEMFILE,SPECFILE,RATEFILE,  & 
                          SPECFILE_KPP,      &
                          SKELETAL_CHEM_KPP, &
                          SKELETAL_SPEC_KPP, &
                          DNUFILE
      !INPUT/OUTPUT DIRS
      PARAMETER(INDIR='../input_test/')
      PARAMETER(OUTDIR='../output/')
      PARAMETER(RATEDIR=INDIR//'rates/')
      !INPUT FILES
      PARAMETER(TRGFILE="target.txt")        !TARGET SPEC
      PARAMETER(DNUFILE="DELTANU.dat")       !NUP-NUR FILE
      PARAMETER(SPECFILE="species.txt")      !SPEC FILE
      PARAMETER(CHEMFILE="reactions.txt")    !REAC FILE 
      PARAMETER(RATEFILE="speciesRR.dat")         !RATES FILE
      PARAMETER(SPECFILE_KPP="SPEC_KPP.txt") !INPUT KPP SPEC FILE
      !OUTPUT FILES
      PARAMETER(SKELETAL_CHEM_KPP=OUTDIR//'/'//'reactions_reduced.txt')
      PARAMETER(SKELETAL_SPEC_KPP=OUTDIR//'/'//'species_reduced.txt')
!
!-----------------------------------------------------------------------
!
      CHARACTER(LEN=NSMX) :: CREACNM_S(NREMX),CREACNM_F(NREMX)
      CHARACTER(LEN=NSMX) :: CSPECNM_S(NSPMX),CSPECNM_F(NSPMX)
      INTEGER :: I,J,K,N,NSPEC,NREAC,LCREAC,LCSPEC,IPROD
!
      CHARACTER(LEN=NSMX), ALLOCATABLE :: CREAC(:)
      CHARACTER(LEN=NSMX), ALLOCATABLE :: CSPEC(:),CSPEC_KPP(:)
      INTEGER, ALLOCATABLE :: RSPEC(:,:),SET_TRG(:,:),SPSET_TRG(:,:), &
       SPSET_TRGUP(:,:),SPSET_UNION(:),RESET_UNION(:)

      REAL(KIND=DBL_P), ALLOCATABLE :: WIJ(:,:),DELTANU(:,:),RR(:), &
                                       JIJW(:,:)
      CHARACTER(LEN=8)  :: DATE
      CHARACTER(LEN=10) :: TIME
      CHARACTER(LEN=5)  :: ZONE
      CHARACTER(LEN=8+10+LEN(OUTDIR)+1) :: DIRCASE
      CHARACTER(LEN=500) :: COMMAND
      CHARACTER(LEN=2) :: JCASE
      CHARACTER(LEN=LEN(RATEDIR)+10) :: FLCASE
!
!-----------------------------------------------------------------------

      WRITE(*,*) 'MAIN:'

      FLCASE=RATEDIR//'output_01/'

      !READ IN ORIGINAL CHEM FILE/SPEC FILE AND KPP PRODUCED SPEC FILE.
      CALL READ_CHEM(INDIR,CHEMFILE,SPECFILE)

      !TODO
      !READ_CONTROL(TRGFILE)
      STOP

      !READ IN DELTANU *SPEC NUMBERING AS IN SPEC_KPP.txt
      !ALLOCATE( DELTANU(NSPEC,NREAC) )
      !OPEN(UNIT=3,FILE=FLCASE//DNUFILE, &
      !     STATUS='OLD',FORM='UNFORMATTED') 
      !READ(3) DELTANU
      !CLOSE(3)
                                !CREACNM_S,CSPECNM_S,CREACNM_F,CSPECNM_F
      !ALLOCATE( CSPEC_KPP(NSPEC) )	
      !CALL PARSE_SPEC_KPP(NSPEC,FLCASE,SPECFILE_KPP,CSPEC_KPP)
      !*THIS IS USED FOR SPECIES ORDERING SINCE KPP OUTPUT SPECIES 
      !ORDERING IS NOT SAME AS KPP INPUT SPECIES ORDERING.


      

      WRITE(*,'(AXI6)') 'NSPEC=',NSPEC
      WRITE(*,'(AXI6)') 'NREAC=',NREAC
      WRITE(*,'(AXI6)') 'NCASE=',NCASE
      WRITE(*,'(AXI6)') 'NSETS=',NDATA
 
      IF(IRDMETH.EQ.1) THEN
       WRITE(*,*) 'REDUCTION USING DRG WITH DFS'
      ELSEIF(IRDMETH.EQ.2) THEN
       WRITE(*,*) 'REDUCTION USING DRGEP WITH DIJIKSTRAS ALGORITHM'
      ELSEIF(IRDMETH.EQ.3) THEN
       WRITE(*,*) 'REDUCTION USING JAC-DRGEP WITH DIJIKSTRAS ALGORITHM'
      ENDIF

      WRITE(*,*) '-------------------'

      ALLOCATE( WIJ(NREAC,NSPEC) )
      ALLOCATE( RR(NREAC) )
      ALLOCATE( JIJW(NSPEC,NSPEC) )
 
      !LOOP THROUGH DATA SETS-PERFORM DRG-UPDATE SPEC AND REAC

      ALLOCATE( SET_TRG(NTRG,NSPEC) )
      ALLOCATE( SPSET_TRGUP(NTRG,NSPEC) )
      ALLOCATE( SPSET_UNION(NSPEC) )
      ALLOCATE( RESET_UNION(NREAC) )
      SPSET_TRGUP(1:NTRG,1:NSPEC)=0
      SPSET_UNION(1:NSPEC)=0
      RESET_UNION(1:NREAC)=0

      WRITE(*,*) 'DRG-LOOPING THROUGH DATASETS'       
      WRITE(*,*) 

      !CASE LOOP
      DO N=1,NCASE
      
       WRITE(*,*) 'CASE',N

       WRITE(JCASE,'(I2.2)') N
       FLCASE=RATEDIR//'output_'//JCASE//'/'
         
       !DATA LOOP
       DO I=1,NDATA

        WRITE(*,*) 'DATASET',I
        WRITE(*,*) '-------------------'

        !WRITE(*,*) 'CHECK JIJW ',JIJW

        CALL READ_RATES(I,FLCASE,RATEFILE,NSPEC,NREAC,WIJ,RR,JIJW)                      
        !1=MY WAY
        !2=GRAPH SEARCH (PATH DEPENDENT)
        CALL DRIVER_DRG(IRDMETH,NSPEC,NREAC,NTRG,INDX_TRG,ETOL, &
                        DELTANU,RSPEC,WIJ,RR,JIJW, &
                        LCSPEC,LCREAC,CSPEC,CREAC,SET_TRG)  
        !SAVE PICs

        !UPDATE_SETS
        DO K=1,NSPEC
         DO J=1,NTRG
          IF(SET_TRG(J,K).EQ.1.AND.SPSET_TRGUP(J,K).EQ.0) THEN
           SPSET_TRGUP(J,K)=1
          ENDIF
         ENDDO  
        ENDDO   
       
        
       ENDDO!DATA

      ENDDO!CASE
      !-----------------------------------------------------------------     
      
      WRITE(*,*) 'FINAL TARGET SPECIES SETS:'
      DO I=1,NTRG
       WRITE(*,*) I,TRIM(CSPEC(INDX_TRG(I))),SUM(SPSET_TRGUP(I,:))
      ENDDO
      WRITE(*,*) '-------------------'

      !GET UNION SET
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

      !GET PRELIMINARY REACTION SET
      WRITE(*,*) 'EXTRACTING PRELIMINARY SKEL MECH REACTIONS'
      WRITE(*,*) '------------------------------------------'
      DO J=1,NREAC
       IPROD=1
       WRITE(*,*) J,TRIM(CREAC(J))
       DO I=1,NSPEC
        IF(RSPEC(J,I).EQ.1) THEN !SPEC IN R
         IPROD=IPROD*SPSET_UNION(I)
         WRITE(*,*) TRIM(ADJUSTL(CSPEC(I))),SPSET_UNION(I),IPROD
        ENDIF
       ENDDO
       IF(IPROD.EQ.1) THEN
        RESET_UNION(J)=1
        WRITE(*,*) '------------------------------KEEP'
       ELSE 
        WRITE(*,*) '------------------------------REMOVE'
       ENDIF
      ENDDO

      !GOTO 1 
      !CHECK REACTIONS
      WRITE(*,*)
      WRITE(*,*) 'CHECKING PRELIMINARY SKEL MECH REACTIONS'
      WRITE(*,*)
      N=0
      DO I=1,NSPEC
       IF(SPSET_UNION(I).EQ.1) THEN !UNION SPEC
        
        DO J=1,NREAC
         IF(RESET_UNION(J).EQ.1.AND.RSPEC(J,I).EQ.1) THEN !EXISTS IN SOME REACTION
          !WRITE(*,*) TRIM(ADJUSTL(CREAC(J))) 
          CHECK=1
          !WRITE(*,'(A,X,A,X,I1)') 'CHECK ',CSPEC(I),CHECK
          EXIT
         ELSE
          CHECK=0
         ENDIF
        ENDDO 
        IF(CHECK.EQ.0) THEN !SPEC IN SKEL SET NOT IN ANY REACTION IN SKEL SET -> REMOVE SPEC
         N=N+1
         !WRITE(*,*) 'ERROR: REDUCED REACTION SET NOT POSSIBLE'
         !WRITE(*,*) TRIM(ADJUSTL(CREAC(J))) 
         !STOP
         WRITE(*,*) 'CHECK REACTIONS-REMOVING ', &
                     TRIM(ADJUSTL(CSPEC(I))), &
                     ' FROM SKEL MECH SET'         
         SPSET_UNION(I)=0         
        ENDIF
       ENDIF
      ENDDO
      IF(N.NE.0) THEN
       WRITE(*,*) 'WARNING: SOME SKEL MECH SPECIES REMOVED'
      ENDIF
      WRITE(*,*) 

1     CONTINUE

      WRITE(*,*) 'ELIMINATED SPECIES:'
      DO I=1,NSPEC
       IF(SPSET_UNION(I).EQ.0) THEN
        WRITE(*,*) TRIM(CSPEC(I)) 
       ENDIF
      ENDDO
      WRITE(*,*) '-------------------'

      WRITE(*,*) 'SPECIES UNION SET:'
      DO I=1,NSPEC
       IF(SPSET_UNION(I).EQ.1) THEN 
        WRITE(*,*) TRIM(CSPEC(I))
       ENDIF
      ENDDO
      WRITE(*,*) '-------------------'

      WRITE(*,*) 'REACTIONS UNION SET:'
      DO I=1,NREAC
       IF(RESET_UNION(I).EQ.1) THEN   
        WRITE(*,*) TRIM(CREAC(I))
       ENDIF
      ENDDO
      WRITE(*,*) '-------------------'
      

      NSPEC_SKEL=SUM(SPSET_UNION)
      NREAC_SKEL=SUM(RESET_UNION)
      IF(NSPEC_SKEL.LE.0.OR.NREAC_SKEL.LE.0) THEN
       WRITE(*,*) 'ERROR: NO SKEL MECH CREATED. TERMINATING ...'
       STOP
      ENDIF
     
      !OUTPUT SKELETAL MECHANISM
      CALL OUTPUT(NSPEC,NREAC,SPSET_UNION,RESET_UNION,CSPEC,CSPECNM_F, &
                  CREACNM_F,SKELETAL_SPEC_KPP,SKELETAL_CHEM_KPP	)


      DEALLOCATE(CREAC)
      DEALLOCATE(CSPEC)	
      DEALLOCATE(RSPEC)
      DEALLOCATE(WIJ)
      DEALLOCATE(INDX_TRG)
      DEALLOCATE(SET_TRG)
      DEALLOCATE(ETOL)

     
      WRITE(*,*) 'SKELETAL MECHANISM SIZE:'
      WRITE(*,*) 'NUNION SPEC=',SUM(SPSET_UNION)
      WRITE(*,*) 'NUNION REAC=',SUM(RESET_UNION)
      WRITE(*,*) 
      WRITE(*,*) 'RUN COMPLETED SUCCESSFULLY!'
      WRITE(*,*) '-------------------'

      CALL SYSTEM('cp log_redchem.txt output/')
      COMMAND='cp'//' '//INDIR//'target.txt'//' '//'output/'
      CALL SYSTEM(COMMAND)
!
      STOP
      END
