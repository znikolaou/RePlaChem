      SUBROUTINE FLSHL(L,STRING,LS)
! 
!     AUTHOR: Z.M. NIKOLAOU
! 
!     DESCRIPTION:
!     REMOVES BLANCK SPACES IN STRING AND FLUSHES LEFT
! 
!     INPUT:
!     L ~STRING LENGTH
!     STRING ~ CHARACTER
!     OUTPUT:
!     STRING
!     LS     ~LENGTH OF FLUSHED STRING
!
      IMPLICIT NONE
!
      INTEGER :: L,LS
      CHARACTER(LEN=L) :: STRING,STRWRK
      !   
      INTEGER :: I,J
!       
      J=0
      STRWRK=' ' 
      DO I=1,L
       IF(STRING(I:I).NE.' ') THEN
        J=J+1
        STRWRK(J:J)=STRING(I:I)
       ENDIF 
      ENDDO
      STRING(1:J)=STRWRK 
      STRING(J+1:L)=' '
      LS=J
! 
      END SUBROUTINE
!----------------------------------------------------------------------
!
      INTEGER FUNCTION NOSTR(STRING,TEXT)
!
!     AUTHOR: ZACHARIAS M. NIKOLAOU     
!
!     DESCRIPTION: COUNTS NO OF OCCURANCES OF TEXT IN 
!     STRING, OF THE FORM TEXT' ' AND ' 'TEXT' '  
!
      IMPLICIT NONE
      CHARACTER(LEN=*) :: STRING
      CHARACTER(LEN=*) :: TEXT
      INTEGER :: I,NSTR,NTXTLR,NTXTR
!
      NSTR=LEN(STRING)
      NTXTLR=LEN(TRIM(TEXT))+2
      NTXTR=LEN(TRIM(TEXT))+1
      
      !CHECK FOR 1ST ELEMENT TEXT' '
      NOSTR=0
      IF( STRING(1:1+NTXTR-1).EQ.(TEXT(1:NTXTLR-2)//' ') ) THEN
       NOSTR=NOSTR+1
      ENDIF
      !CHECK FOR ' 'TEXT' '
      DO I=2,NSTR-NTXTLR+1
       IF( (STRING(I:I+NTXTLR-1)).EQ.(' '//TEXT(1:NTXTLR-2)//' ') ) THEN !TEXT IN STRING
        NOSTR=NOSTR+1
       ENDIF
      ENDDO
!
      END FUNCTION NOSTR
!-----------------------------------------------------------------------
!
      INTEGER FUNCTION INDX_TXT(STRING,TEXT)
!
!     AUTHOR: ZACHARIAS M. NIKOLAOU     
!
!     DESCRIPTION: 
!     RETURNS INDEX OF FIRST OCCURENCE OF 'TEXT' IN 
!     'STRING'.  
!
      IMPLICIT NONE
      CHARACTER(LEN=*) :: STRING
      CHARACTER(LEN=*) :: TEXT
      INTEGER :: I,NSTR,NTXT,NOSTR
!
      NSTR=LEN(TRIM(STRING))
      NTXT=LEN(TRIM(TEXT))
      
      NOSTR=0
      DO I=1,NSTR-NTXT+1
       IF( (STRING(I:I+NTXT-1)).EQ.(TEXT) ) THEN !TEXT IN STRING
        NOSTR=NOSTR+1
       ENDIF
      ENDDO
      INDX_TXT=NOSTR
!
      END FUNCTION INDX_TXT
!----------------------------------------------------------------------
!
      INTEGER FUNCTION NOTXT(STRING,TEXT)
!
!     AUTHOR: ZACHARIAS M. NIKOLAOU     
!
!     DESCRIPTION: 
!     COUNTS NO OF OCCURANCES OF TEXT IN 
!     STRING  
!
      IMPLICIT NONE
      CHARACTER(LEN=*) :: STRING
      CHARACTER(LEN=*) :: TEXT
      INTEGER :: I,NSTR,NTXT,N
!
      NSTR=LEN(TRIM(STRING))
      NTXT=LEN(TRIM(TEXT))
      
      !CHECK FOR 'TEXT'
      N=0
      DO I=1,NSTR-NTXT+1
       IF( (STRING(I:I+NTXT-1)).EQ.(TEXT) ) THEN !TEXT IN STRING
        N=N+1
       ENDIF
      ENDDO
      NOTXT=N
!
      END FUNCTION NOTXT
!----------------------------------------------------------------------
!
      FUNCTION STR2NUM(C)
!
!     AUTHOR: ZACHARIAS M. NIKOLAOU
!
!     DESCRIPTION:
!     RETURNS REAL FROM CHARACTER      
!
!      USE PRECIS, ONLY: DBL_P
!
      IMPLICIT NONE
!
      INTEGER :: LENC
      CHARACTER(LEN=*) :: C
!      REAL(KIND=DBL_P) :: CHAR2R
      REAL :: STR2NUM,A
!
!
      READ(C,'(G12.5)') A
     
      STR2NUM=A                

      END FUNCTION STR2NUM
!----------------------------------------------------------------------
      SUBROUTINE RPLTXT(STRING,TEXT,NEWTEXT)
!
!     AUTHOR: ZACHARIAS M. NIKOLAOU     
!
!     DESCRIPTION: 
!     REPLACES ALL OCCURANCES OF TEXT IN 
!     STRING WITH NEWTXT  
!
      IMPLICIT NONE
      CHARACTER(LEN=*) :: STRING
      CHARACTER(LEN=*) :: TEXT,NEWTEXT
      INTEGER :: I,NSTR,NTXT,NNEWTXT
!
      !CHECK LENGTH OF NEWTEXT
      IF( LEN(NEWTEXT).GT.LEN(TEXT) ) THEN
       WRITE(*,*) 'FUNCTION RPLTXT: ERROR'
       STOP
      ENDIF

      NSTR=LEN(TRIM(ADJUSTL(STRING)))
      NTXT=LEN(TRIM(ADJUSTL(TEXT)))
      NNEWTXT=LEN(TRIM(ADJUSTL(NEWTEXT)))
      
      !CHECK FOR 'TEXT'
      DO I=1,NSTR-NTXT+1
       IF( STRING(I:I+NTXT-1).EQ.TEXT ) THEN !TEXT IN STRING
        STRING(I:I+NTXT-1)=' '
        STRING(I:I+NNEWTXT-1)=NEWTEXT
       ENDIF
      ENDDO
!
      END SUBROUTINE
!----------------------------------------------------------------------
     SUBROUTINE RPLTXTE(STRING,TEXT,NEWTEXT,STRWRK,LWRK)
!
!     AUTHOR: ZACHARIAS M. NIKOLAOU     
!
!     DESCRIPTION: 
!     REPLACES ALL OCCURANCES OF TEXT IN 
!     STRING WITH NEWTXT  
!
      IMPLICIT NONE
      INTEGER :: LSTR,LTXT,LTXTN,LWRK       
      CHARACTER(LEN=*) :: STRING
      CHARACTER(LEN=*) :: TEXT
      CHARACTER(LEN=*) :: NEWTEXT
      CHARACTER(LEN=LWRK) :: STRWRK
!
      INTEGER :: I,IL,IR,IS,DIFF,NOCC,LENNEW
!
      LSTR=LEN(STRING)
      LTXT=LEN(TEXT)
      LTXTN=LEN(NEWTEXT)
      
      STRWRK(1:LWRK)=' '
      DIFF=LTXTN-LTXT       
      IS=LEN(TRIM(ADJUSTL(STRWRK))) 

!      WRITE(*,*) LSTR,LTXT,LTXTN,IS
!      WRITE(*,*) STRING

      NOCC=0
      IL=1
      DO I=1,LSTR-LTXT+1
       IF( (STRING(I:I+LTXT-1)).EQ.(TEXT) ) THEN !TEXT IN STRING
        NOCC=NOCC+1
        IR=I-1

!        WRITE(*,*) IS,IL,IR
!        WRITE(*,*) STRING(IL:IR)
        STRWRK(IS+1:)=STRING(IL:IR)//NEWTEXT
        IS=IS+(IR-IL+1)+LTXTN!LEN(TRIM(ADJUSTL(STRWRK))) 
        IL=I+LTXT         
!        WRITE(*,*) STRWRK

       ENDIF              
      ENDDO

      !NEW LENGTH
      LENNEW=LEN(TRIM(ADJUSTL(STRING)))+DIFF*NOCC
      IF(LENNEW.GT.LWRK) THEN
       WRITE(*,*) 'RPLTXTE: NOT ENOUGH WORSPACE PROVIDED'
       STOP
      ENDIF

      !LAST PART
      STRWRK(IS+1:)=STRING(IL:LSTR)
!
      END SUBROUTINE
!-----------------------------------------------------------------------
