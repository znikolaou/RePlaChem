      SUBROUTINE OUTPUT(NSPEC,NREAC,SETSP,SETRE,CSPNM,CSPNMF, &
                        CRENMF,FLSP_KPP,FLRE_KPP)
!
      USE GLOBAL, ONLY :NSMX
      IMPLICIT NONE 
!
      INTEGER :: I,J,K,INDX,NSPEC,NREAC,LS
      INTEGER :: SETSP(NSPEC),SETRE(NREAC)
      CHARACTER(LEN=*) :: CSPNM(NSPEC),CSPNMF(NSPEC),CRENMF(NREAC), &
                          FLSP_KPP,FLRE_KPP
      CHARACTER(LEN=NSMX) :: C,CF
      INTEGER :: NOSTR,NSPEC_SKEL,NREAC_SKEL
      CHARACTER(LEN=8)  :: DATE
      CHARACTER(LEN=10) :: TIME
      CHARACTER(LEN=5)  :: ZONE

!     
      CALL DATETIME(DATE,TIME,ZONE)

      NSPEC_SKEL=SUM(SETSP)  
      NREAC_SKEL=SUM(SETRE)  

      !REACTION FILE
      OPEN(UNIT=1,FILE=FLRE_KPP,STATUS='NEW',FORM='FORMATTED')

      WRITE(1,'(A)') '// AUTHOR: Z.M. NIKOLAOU'
      WRITE(1,'(A)') '// SKELETAL MECHANISM-OUTPUT FROM DRG.'
      WRITE(1,'(4(AX))') '// CREATED: ',TRIM(DATE),TRIM(TIME),TRIM(ZONE)
      WRITE(1,'(A)') '//' 
      WRITE(1,'(A)') '#EQUATIONS' 

      J=0
      DO I=1,NREAC
       IF(SETRE(I).EQ.1) THEN
        J=J+1
        C=CRENMF(I)
        !WRITE(*,*) C
        !C(1:INDEX(C,'}'))=' '
        C=TRIM(ADJUSTL(C))
        WRITE(1,'(A,I4.4,A,X,A)') '{',J,'}',TRIM(ADJUSTL(C))
       ENDIF 
      ENDDO
      !
      GOTO 1

      !OUTPUT DUMMY REACTIONS FOR ELIMINATED SPECIES
      WRITE(1,'(A)') '{DUMMY REACTIONS-ELIMINATED SPEC}'
      DO I=1,NSPEC
       IF(SETSP(I).EQ.0) THEN
        C=CSPNM(I)
        WRITE(1,'(2(A,X,A,X))') TRIM(ADJUSTL(C)),'=',TRIM(ADJUSTL(C)), &
                           ': 0.0D0 ;' 
       ENDIF
      ENDDO
      !
      WRITE(1,'(A)') '//' 
      CLOSE(1)

1     CONTINUE

      !SPECIES FILE
      OPEN(UNIT=1,FILE=FLSP_KPP,STATUS='NEW',FORM='FORMATTED')

      WRITE(1,'(A)') '// AUTHOR: Z.M. NIKOLAOU'
      WRITE(1,'(A)') '// SKELETAL MECHANISM-OUTPUT FROM DRG.'
      WRITE(1,'(4(AX))') '// CREATED: ',TRIM(DATE),TRIM(TIME),TRIM(ZONE)
      WRITE(1,'(A)') '//' 
      WRITE(1,'(A)') '#DEFVAR' 
      DO J=1,NSPEC
       C=CSPNM(J)
       CALL FLSHL(LEN(C),C,LS)
       
       IF(SETSP(J).EQ.1) THEN
        !WRITE(*,'(A,X,I3)') C
        K=0
        I=0         
        DO WHILE(K.EQ.0.AND.I.LE.NSPEC-1)
         I=I+1
         CF=CSPNMF(I)
         CF=CF(1:INDEX(CF,'='))
         K=NOSTR(CF,C)         
         !K=NOSTR(CSPNMF(I),C)
         !write(*,*) K
        ENDDO            
        !WRITE(*,*) CSPNMF(I)
        WRITE(1,'(A)') TRIM(ADJUSTL(CSPNMF(I)))
       ENDIF
      ENDDO
      WRITE(1,'(A)') '//' 
      CLOSE(1)

      END SUBROUTINE
