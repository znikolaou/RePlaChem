      SUBROUTINE READ_CHEM(INDIR,RATEDIR,CHEMFILE,SPECFILE,DNUFILE, &
                           CREAC_S,CSPEC_S, &
                           CREAC_F,CSPEC_F, &
                           NREAC,NSPEC)
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
      INTEGER I,NREAC,NSPEC,I_REAC(NWRKC,NWRKC),I_PROD(NWRKC,NWRKC)
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
      STOP
      
      CLOSE(1)
      CLOSE(2)
      
      !REMOVE SPECIES APPEAR. ONLY AS PRODUCTS
      CALL GET_REACPROD(NWRKC,NWRKC,NSTR_SPMX,NSTR_REMX, &
                        CSPEC_S,CREAC_S,I_REAC,I_PROD)


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
      SUBROUTINE GET_RSPEC(NSPEC,NREAC,LCSPEC,LCREAC, &
                           CSPECNM,CREACNM,RSPEC)
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
      INTEGER :: I,J,K,A,ILEN_SPEC(NSPEC),RSPEC(NREAC,NSPEC)
      INTEGER :: NOSTR
!     
      !EXTRACT SPECIES STRING LENGTH
      !DO I=1,NSPEC
      ! ILEN_SPEC(I)=INT(LEN(TRIM(CSPECNM(I))))
      !ENDDO
     
      DO J=1,NSPEC
       DO I=1,NREAC
        A=NOSTR(CREACNM(I),CSPECNM(J)) 
        A=MIN(1,A)
        RSPEC(I,J)=A
       ENDDO
      ENDDO

      END SUBROUTINE
!
!-----------------------------------------------------------------------
 
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
        II=INDEX(C,'=')
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
!
!-----------------------------------------------------------------------
