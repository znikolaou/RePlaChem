      SUBROUTINE READ_CHEM(INDIR,RATEDIR,CHEMFILE,SPECFILE,DNUFILE, &
                           CREAC_S,CSPEC_S,CREAC_F,CSPEC_F,NREAC,NSPEC)
!
!     AUTHOR: ZACHARIAS M. NIKOLAOU
! 
!     DESCRIPTION: PARSES KPP-STYLE CHEM AND SPEC FILES.
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
      USE GLOBAL, ONLY : NWRKC,NSTR_SPMX,NSTR_REMX
      USE PRECIS, ONLY: DBL_P
!
      IMPLICIT NONE
!
      CHARACTER(LEN=*) :: INDIR,RATEDIR,CHEMFILE,SPECFILE,DNUFILE
      CHARACTER(LEN=NSTR_REMX) :: CREAC_S(NWRKC),CREAC_F(NWRKC)
      CHARACTER(LEN=NSTR_SPMX) :: CSPEC_S(NWRKC),CSPEC_F(NWRKC)
      INTEGER I,NREAC,NSPEC!,I_REAC(NWRKC,NWRKC),I_PROD(NWRKC,NWRKC)
!

      !REMOVE TABS FROM INPUT FILES
      CALL SYSTEM("sed -i 's/\t/\ /g'"//' '//INDIR//CHEMFILE) 
      CALL SYSTEM("sed -i 's/\t/\ /g'"//' '//INDIR//SPECFILE) 

      OPEN(UNIT=1,FILE=INDIR//CHEMFILE,STATUS='OLD',FORM='FORMATTED')
      OPEN(UNIT=2,FILE=INDIR//SPECFILE,STATUS='OLD',FORM='FORMATTED')

      !PARSE REACTIONS: KPP FORMAT
      !CALL PARSE_CHEM(1,NSTR_REMX,NWRKC,NREAC,CREAC_S,CREAC_F)

      !PARSE REACTIONS: SIMPLE FORMAT
      CALL PARSE_INDEX_AND_STRING(1,NSTR_REMX,NWRKC,NREAC,CREAC_S)
      WRITE(*,*) 'PARSE_CHEM: REACTIONS:' 
      DO I=1,NREAC
       WRITE(*,'(I5,X,A)') I,TRIM(ADJUSTL(CREAC_S(I)))
      ENDDO 

      !PARSE SPECIES: KPP FORMAT
      !CALL PARSE_SPEC(2,NSTR_SPMX,NWRKC,NSPEC,CSPEC_S,CSPEC_F)
      !PARSE SPECIES: SIMPLE FORMAT
      CALL PARSE_INDEX_AND_STRING(2,NSTR_SPMX,NWRKC,NSPEC,CSPEC_S)
      WRITE(*,*) 'PARSE_CHEM: SPECIES:' 
      DO I=1,NSPEC
       WRITE(*,'(I5,X,A)') I,TRIM(ADJUSTL(CSPEC_S(I)))
      ENDDO 
      
      CLOSE(1)
      CLOSE(2)
       
      !REMOVE SPECIES APPEAR. ONLY AS PRODUCTS
      !CALL GET_REACPROD(NWRKC,NWRKC,NSTR_SPMX,NSTR_REMX, &
      !                  CSPEC_S,CREAC_S,I_REAC,I_PROD)

      !REACTION SPECIES
      !TODO
      !CALL GET_REACTION_SPECIES(NSPEC,NREAC,CSPEC_S,CREAC_S,REAC_SPEC)


      !SPEC_NU_MATRIX
      !TODO
      !CALL GET_SPEC_NU_MATRIX(NSPEC,NREAC,CSPEC_S,CREAC_S,SPEC_NU_MATRIX)
      !SPEC_NU_MATRIX(NSPEC,NREAC,1) !REAC
      !SPECI_NU_MATRIX(NSPEC,NREAC,2) !PROD

      STOP
      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE PARSE_INDEX_AND_STRING(FILEID,NMX,NWRKC,NLIST,STRLIST)
      IMPLICIT NONE
      LOGICAL ISEMPTY,ISCOMMENT
      INTEGER :: FILEID,NMX,NWRKC,NLIST
      CHARACTER(LEN=NMX) :: STRLIST(NWRKC),C
      INTEGER :: STAT,I,IDUM,INDXR,IS,INDX_SPACE
      
      I=0
      STAT=0
      DO WHILE(STAT.EQ.0)
       READ(FILEID,'(A)',IOSTAT=STAT) C  !FORMAT IS: INTEGER STRING
       IF(.NOT.ISEMPTY(C).AND.(.NOT.ISCOMMENT(C))) THEN
        C=TRIM(ADJUSTL(C))
        IS=INDX_SPACE(C)
        I=I+1
        STRLIST(I)=C(IS+1:)
       ENDIF
      ENDDO      
      NLIST=I

      END SUBROUTINE 
      !----------------------------------------------------------------
      SUBROUTINE PARSE_CHEM(ID,NREMX,NWRKC,NREAC,CREACS,CREACF)
!
!     AUTHOR: ZACHARIAS M. NIKOLAOU
!
!     DESCRIPTION:
!     READS CHEM.DAT FILE, AND RETURNS NUMBER OF REACTIONS AND
!     REACTION NAMES
!    
!     INPUT: 
!     ID ~FILE ID
!     NWRKC ~CHARACTER ARRAY SIZE
!     
!     OUTPUT:
!     NREAC ~NUMBER OF REACTIONS
!     CREAC ~REACTION NAMES
! 
      IMPLICIT NONE
      INTEGER :: ID,NREMX,NWRKC,NREAC
      CHARACTER(LEN=NREMX) :: CREACS(NWRKC),CREACF(NWRKC)
      !
      CHARACTER(LEN=NREMX) :: C,CWRK
      INTEGER :: I,STAT,IREAC,INDXR
      
      I=0
      STAT=0
      IREAC=0
      DO WHILE( STAT.EQ.0 )
       I=I+1
       READ(ID,'(A)',IOSTAT=STAT) C
        
       C=TRIM(ADJUSTL(C))
       IF(C(1:1).EQ.'#') THEN !COMMENT IGNORE
        !WRITE(*,*) 'COMMENT-LINE',I         
       ELSE
        IF(C.EQ.' ') THEN !EMPTY SPACE IGNORE
         !WRITE(*,*) 'EMPTY SPACE-LINE',I 
        ELSE
         IF(C(1:1).EQ.'{') THEN !EQUATION
          IREAC=IREAC+1
          INDXR=INDEX(C,'}') 
          !REMOVE COMMENTS
          CREACF(IREAC)=TRIM(ADJUSTL(C(INDXR+1:NREMX)))          
          !WRITE(*,*) CREACF(IREAC)
          DO WHILE( (INDEX(C(INDXR+1:NREMX),':').EQ.0)  )
           READ(ID,'(A)') C
           CREACF(IREAC)=TRIM(ADJUSTL(CREACF(IREAC)))//' '// &
                        TRIM(ADJUSTL(C))//' '
          ENDDO
          !WRITE(*,*) IREAC,TRIM(CREACF(IREAC))
         ENDIF         
        ENDIF
       ENDIF

      ENDDO      
      NREAC=IREAC

      !REMOVE RATE CONSTANTS
      DO I=1,NREAC
       C=CREACF(I)
       C(INDEX(C,':'):NREMX)=' '
       !REMOVE BLANCS AND FLUSH LEFT FIRST 
       !CALL FLSHL(LEN(C),C,LS)
       CALL RPLTXTE(C,'+',' + ',CWRK,NREMX)
       C=CWRK
       CALL RPLTXTE(C,'=',' = ',CWRK,NREMX)
       C=CWRK
       CALL RPLTXTE(C,'-',' - ',CWRK,NREMX)
       C=CWRK
       CALL RPLTXTE(C,'{',' { ',CWRK,NREMX)
       C=CWRK
       CALL RPLTXTE(C,'}',' } ',CWRK,NREMX)
       CREACS(I)=CWRK
      ENDDO


      END SUBROUTINE
!----------------------------------------------------------------------
!
      SUBROUTINE PARSE_SPEC(ID,NSPMX,NWRKC,NSPEC,CSPECS,CSPECF)
!
!     AUTHOR: ZACHARIAS M. NIKOLAOU
!
!     DESCRIPTION:
!     READS CHEM.DAT FILE, AND RETURNS NUMBER OF REACTIONS AND
!     REACTION NAMES
!    
!     INPUT: 
!     ID ~FILE ID
!     NWRKC ~CHARACTER ARRAY SIZE
!     
!     OUTPUT:
!     NSPEC ~NUMBER OF REACTIONS
!     CSPEC ~REACTION NAMES
 
      IMPLICIT NONE
      INTEGER :: ID,NSPMX,NWRKC
      !      
      CHARACTER(LEN=NSPMX) :: C,CSPECS(NWRKC),CSPECF(NWRKC)
      INTEGER :: I,STAT,ISPEC,NSPEC
      
      I=0
      STAT=0
      ISPEC=0
      DO WHILE(STAT.EQ.0)
       I=I+1
       READ(ID,'(A)',IOSTAT=STAT) C
       !WRITE(*,*) C
       IF(C(1:1).EQ.'#') THEN !DIRECTIVE IGNORE
        !WRITE(*,*) I, C
       ELSEIF(C(1:1).EQ.'{') THEN !COMMENT IGNORE                        
        !WRITE(*,*) I, C
       ELSEIF(C.EQ.' ') THEN !EMPTY SPACE IGNORE
        !WRITE(*,*) I, C
       ELSE
        !WRITE(*,*) I,TRIM(ADJUSTL(C))
         ISPEC=ISPEC+1 
         CSPECF(ISPEC)=TRIM(ADJUSTL(C))
         !WRITE(*,*) ISPEC         
       ENDIF
       !WRITE(*,*) I, STAT
      ENDDO
      NSPEC=ISPEC
      !WRITE(*,*) 'NUMBER OF SPECIES=',NSPEC
      !STOP
      !REMOVE EVERYTHING AFTER =
      DO I=1,NSPEC
       C=CSPECF(I)
       C(INDEX(C,'='):NSPMX)=' '
       CSPECS(I)=C
       !WRITE(*,*) I,CSPECS(I)
      ENDDO

      END
!
!-----------------------------------------------------------------------
!     
      SUBROUTINE PARSE_SPEC_KPP(NSPEC,RATEDIR,FLNM,CSPECS)
! 
      USE GLOBAL, ONLY : NSTR_SPMX
!
      IMPLICIT NONE
      INTEGER :: I,NSP,NSPEC
      CHARACTER(LEN=*) :: FLNM,RATEDIR
      CHARACTER(LEN=NSTR_SPMX) :: CSPECS(NSPEC)
      !      
      CHARACTER(LEN=NSTR_SPMX) :: C
!      
      OPEN(UNIT=1,FILE=RATEDIR//FLNM,STATUS='OLD')
      READ(1,*)
      READ(1,*) NSP
      IF(NSP.NE.NSPEC) THEN
       WRITE(*,*) 'PARSE_SPEC_KPP ERROR: MISMATCH'
       STOP
      ENDIF
      DO I=1,NSPEC
       READ(1,'(A)') C      
       CSPECS(I)=TRIM(ADJUSTL(C))
      ENDDO
    
      CLOSE(1)

      END SUBROUTINE
!
!-----------------------------------------------------------------------
!
      SUBROUTINE GET_SPEC_NU_MATRIX(NSPEC,NREAC,LCSPEC,LCREAC, &
                           CSPECNM,CREACNM,SPEC_NU_MATRIX)
!
!     AUTHOR: ZACHARIAS M. NIKOLAOU
!
!     DESCRIPTION: 
!     EXTRACTS SPECIES J INVOLVED IN REACTION I, RSPEC(I,J).
!
      IMPLICIT NONE
      INTEGER :: NSPEC,NREAC,LCSPEC,LCREAC
      CHARACTER(LEN=LCSPEC), DIMENSION(NSPEC) :: CSPECNM
      CHARACTER(LEN=LCREAC), DIMENSION(NREAC) :: CREACNM      
!
      INTEGER :: I,J,K,A,ILEN_SPEC(NSPEC),SPEC_NU_MATRIX(NREAC,NSPEC)
      INTEGER :: NOSTR!     
     
      DO J=1,NSPEC
       DO I=1,NREAC
        !CLEFT=GET_REACTANTS(NREAC,CREACNM(I))
        !CRIGHT=GET_PRODUCTS(NREAC,CREACNM(I))
        !NUL=NOSTR(CLEFT,CSPECNM(J))
        !NUR=NOSTR(CRIGH,CSPECNM(J))
        !SPEC_NU_MATRIX(I,J,1)=NUL
        !SPEC_NU_MATRIX(I,J,2)=NUR
        !SPEC_N_MATRIX(I,J,3)=NUR-NUL
       ENDDO
      ENDDO

      END SUBROUTINE
      !-----------------------------------------------------------------
      SUBROUTINE SEPARATE_REAC(CREAC,SEP,REAC,PROD)
      IMPLICIT NONE
      INTEGER ISEP
      CHARACTER(LEN=*) CREAC,SEP,REAC,PROD

      ISEP=INDEX(TRIM(ADJUSTL(CREAC)),SEP)
      REAC=CREAC(1:ISEP)
      PROD=CREAC(ISEP+1:)

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE GET_REACPROD(NSPEC,NREAC,LCSPEC,LCREAC, &
                              CSPECNM,CREACNM,I_REAC,I_PROD)

!     AUTHOR: ZACHARIAS M. NIKOLAOU
!
!     DESCRIPTION: 
!     I_REAC(NSPEC) : 1 SPEC IN REACTANTS
!     I_PROD(NSPEC) : 1 SPEC IN PRODUCTS
!
!     CREACNM: A + B = C + D
!
      IMPLICIT NONE
      INTEGER :: NSPEC,NREAC,LCSPEC,LCREAC
      CHARACTER(LEN=LCSPEC), DIMENSION(NSPEC) :: CSPECNM
      CHARACTER(LEN=LCREAC), DIMENSION(NREAC) :: CREACNM
!      
      INTEGER :: I,J,K,A,II,I_REAC(NREAC,NSPEC),I_PROD(NREAC,NSPEC), &
                 NOSTR
      CHARACTER(LEN=LCREAC) :: C
          
      DO J=1,NSPEC
       WRITE(*,*) TRIM(ADJUSTL(CSPECNM(J)))
       DO I=1,NREAC 
        WRITE(*,*) CREACNM(I)

        C=' '
        C=CREACNM(I)
        II=INDEX(C,'=>')
        IF(II.NE.0) THEN

         C(II:LCREAC)=' '
         WRITE(*,*) TRIM(ADJUSTL(C))
         A=NOSTR(C,CSPECNM(J)) 
         A=MIN(1,A)
         I_REAC(I,J)=A
         !
         C=' '
         C=CREACNM(I)
         C(1:II)=' '
         A=NOSTR(C,CSPECNM(J)) 
         A=MIN(1,A)
         I_PROD(I,J)=A
         WRITE(*,*) TRIM(ADJUSTL(C))

         IF(I_REAC(I,J).EQ.1) THEN 
          WRITE(*,*) TRIM(ADJUSTL(CSPECNM(J))),' ','IN REACTANTS'
         ENDIF
         IF(I_PROD(I,J).EQ.1) THEN 
          WRITE(*,*) TRIM(ADJUSTL(CSPECNM(J))),' ','IN PRODUCTS'
         ENDIF


        ENDIF

       ENDDO
      ENDDO

      !STOP

      END SUBROUTINE
      !-----------------------------------------------------------------
      SUBROUTINE SEPARATE_REACTIONS(N,CREAC,REAC,PROD)
      USE GLOBAL, ONLY : NSTR_REMX
      IMPLICIT NONE
      INTEGER :: N,I
      CHARACTER(LEN=*) :: CREAC(N),REAC(N),PROD(N)
      CHARACTER(LEN=NSTR_REMX) :: CLIST(2)

      DO I=1,N
       IF(INDEX(CREAC(I),'bolsig').EQ.0) THEN 
        CALL SEPARATE_STRING(CREAC(I),'=>',CLIST)
       ELSE
        CALL SEPARATE_STRING(CREAC(I),'->',CLIST)
       ENDIF
        REAC(I)=CLIST(1)
        PROD(I)=CLIST(2)
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE GET_REACTION_SPECIES(NSPEC,NREAC,CSPEC,CREAC, &
                                      REAC_SPEC,NOSPEC)
      !USE GLOBAL, ONLY: NSTR_SPMX,NSTR_REMX
      IMPLICIT NONE
      INTEGER :: NSPEC,NREAC,NOSPEC(NREAC),I,J,K,REAC_SPEC(NREAC,NSPEC)
      CHARACTER(LEN=*) :: CREAC(NREAC),CSPEC(NSPEC)
      LOGICAL :: IS_SPEC_IN_REACTION,IS_BOLSIG_REACTION, &
                 IS_ANY_NEUTRAL_REACTION,IS_CHARGED_SPECIES

      REAC_SPEC(1:NREAC,1:NSPEC)=0
      DO J=1,NREAC
       K=0
       DO I=1,NSPEC
        IF(IS_SPEC_IN_REACTION(CREAC(J),CSPEC(I))) THEN
         K=K+1
         REAC_SPEC(J,I)=1
         NOSPEC(J)=K
        ENDIF
       ENDDO
      ENDDO

      !ADD MANUALLY SPEC 'E' NOT PRESENT FOR BOLSIG REACTIONS
      DO I=1,NREAC
       IF(IS_BOLSIG_REACTION(CREAC(I))) THEN
        REAC_SPEC(I,1)=1 !TODO: THIS IS HARD-CODED HERE! 
        NOSPEC(I)=NOSPEC(I)+1
       ENDIF
      ENDDO

      !CHECK FOR 'ANY_NEUTRAL' REACTIONS
      DO I=1,NREAC
       IF(IS_ANY_NEUTRAL_REACTION(CREAC(I))) THEN
        DO J=1,NSPEC
         IF(.NOT.IS_CHARGED_SPECIES(CSPEC(J))) THEN
          REAC_SPEC(I,J)=1
         ENDIF
        ENDDO
       ENDIF
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
      FUNCTION IS_ANY_NEUTRAL_REACTION(REAC)
      IMPLICIT NONE
      CHARACTER(LEN=*) :: REAC,NEUTRAL
      LOGICAL :: IS_ANY_NEUTRAL_REACTION
      PARAMETER(NEUTRAL='ANY_NEUTRAL')

      IS_ANY_NEUTRAL_REACTION=INDEX(REAC,NEUTRAL).GT.0

      END FUNCTION
      !-----------------------------------------------------------------
      FUNCTION IS_BOLSIG_REACTION(REAC)
      IMPLICIT NONE
      LOGICAL :: IS_BOLSIG_REACTION
      CHARACTER(LEN=*) :: REAC
      CHARACTER(LEN=*) :: BOLSIG
      PARAMETER(BOLSIG='bolsig')
      
      IF(INDEX(REAC,BOLSIG).GT.0) THEN
       IS_BOLSIG_REACTION=.TRUE.
      ELSE
       IS_BOLSIG_REACTION=.FALSE.
      ENDIF

      END FUNCTION
      !-----------------------------------------------------------------
      FUNCTION IS_SPEC_IN_REACTION(REAC,SPEC)
      IMPLICIT NONE
      LOGICAL :: CASEA,CASEB,CASEC,IS_SPEC_IN_REACTION
      CHARACTER(LEN=*) :: REAC,SPEC
      CHARACTER(LEN=LEN(TRIM(ADJUSTL(SPEC)))+2) :: CA
      
      CA=' '//TRIM(ADJUSTL(SPEC))//' ' 
      IF(INDEX(REAC,CA).GT.0) THEN
       IS_SPEC_IN_REACTION=.TRUE.
      ELSE
       IS_SPEC_IN_REACTION=.FALSE.
      ENDIF

      END FUNCTION
      !-----------------------------------------------------------------
      SUBROUTINE GET_FORMATTED_REACTIONS(NREAC,CREAC,CREAC_F)
      !
      ! REFORMATS REACTIONS AS BELOW-INTENDED FOR PLASMA CHEMISTRY
      ! 
      ! AUTHOR: Z. NIKOLAOU (2025)     
      USE GLOBAL, ONLY : NSTR_REMX
      IMPLICIT NONE
      INTEGER :: NREAC,I
      CHARACTER(LEN=*) :: CREAC(NREAC)
      CHARACTER(LEN=*) :: CREAC_F(NREAC)
      CHARACTER(LEN=NSTR_REMX) :: C,CWRK

      DO I=1,NREAC
       C=CREAC(I)
       CALL RPLTXTE(C,'^+','^POS',CWRK,NSTR_REMX)
       C=CWRK
       CALL RPLTXTE(C,'^-','^NEG',CWRK,NSTR_REMX)
       C=CWRK
       CALL RPLTXTE(C,'+',' + ',CWRK,NSTR_REMX)
       C=CWRK
       CALL RPLTXTE(C,'=>',' => ',CWRK,NSTR_REMX)
       C=CWRK
       CALL RPLTXTE(C,'->',' -> ',CWRK,NSTR_REMX)
       C=CWRK
       CALL RPLTXTE(C,':',': ',CWRK,NSTR_REMX)
       C=CWRK
       CALL RPLTXTE(C,'(E-V)',' *(E-V)',CWRK,NSTR_REMX)
       C=CWRK
       CALL RPLTXTE(C,'ANY_NEUTRAL','  ANY_NEUTRAL ',CWRK,NSTR_REMX)
       C=CWRK
       CALL RPLTXTE(C,'POS','+',CWRK,NSTR_REMX)
       C=CWRK
       CALL RPLTXTE(C,'NEG','-',CWRK,NSTR_REMX)
       CREAC_F(I)='$ '//TRIM(ADJUSTL(CWRK))//' $'
      ENDDO
      
      RETURN
      END
      !-----------------------------------------------------------------
      FUNCTION IS_CHARGED_SPECIES(SPEC)
      IMPLICIT NONE
      CHARACTER(LEN=*) :: SPEC
      CHARACTER(LEN=LEN(TRIM(ADJUSTL(SPEC)))) :: C
      LOGICAL :: IS_CHARGED_SPECIES

      C=TRIM(ADJUSTL(SPEC))
      IS_CHARGED_SPECIES=(INDEX(C,'^+').GT.0).OR. &
                         (INDEX(C,'^-').GT.0).OR.(C.EQ.'E')
             
      END FUNCTION
      !-----------------------------------------------------------------

