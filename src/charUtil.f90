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
      INTEGER FUNCTION INDX_SPACE(STRING)
      IMPLICIT NONE
      CHARACTER(LEN=*) STRING
      INTEGER I

      DO I=1,LEN(STRING)
       !WRITE(*,*) I,STRING(I:I),STRING(I:I).EQ.''
       IF(STRING(I:I).EQ.'') THEN
        INDX_SPACE=I
        EXIT
       ENDIF
      ENDDO

      END FUNCTION
      !-----------------------------------------------------------------      
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
      !------------------------------------------------------------------
      LOGICAL FUNCTION ISEMPTY(STRING)
      IMPLICIT NONE
      CHARACTER(LEN=*) STRING
       
      ISEMPTY=STRING.EQ.' '

      END FUNCTION
      !-----------------------------------------------------------------
      LOGICAL FUNCTION ISCOMMENT(STRING)
      ! RETURNS TRUE IF FIRST CHAR IS ! OR #       
      IMPLICIT NONE
      CHARACTER(LEN=*) STRING

      ISCOMMENT=STRING(1:1).EQ.'!'.OR.STRING(1:1).EQ.'#'

      END FUNCTION
      !-----------------------------------------------------------------
      SUBROUTINE SPLIT_STRING(STRING,SEP,N,COLMS,NA)
      IMPLICIT NONE
      INTEGER :: N,NA,I,IS,IE,INEXT,ISNEW,INE,NTRIM,LSEP
      CHARACTER(LEN=*) :: STRING,SEP,COLMS(N)
      CHARACTER(LEN=LEN(TRIM(ADJUSTL(STRING)))) :: CSTR

      CSTR=TRIM(ADJUSTL(STRING))
      IS=INDEX(CSTR,SEP)
      IF(IS.EQ.0) THEN
       COLMS(1)=CSTR
       NA=0
      ELSE
       LSEP=LEN(SEP)
       IE=IS+LSEP-1
       COLMS(1)=CSTR(1:IS-1)
       I=1
       INEXT=IS
       DO WHILE(INEXT.GT.0)
        INEXT=INDEX(CSTR(IE+1:),SEP)
        ISNEW=INEXT+IE
        !WRITE(*,*) IS,IE,ISNEW,CSTR(IE+1:ISNEW-1)
        IF(ISNEW.GE.IE) THEN
         I=I+1
         COLMS(I)=CSTR(IE+1:ISNEW-1)
        ENDIF
        IE=ISNEW+LSEP-1
       ENDDO
        COLMS(I)=CSTR(ISNEW+1:)
        NA=I
      ENDIF

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE SEPARATE_STRING(STRING,SEP,STRLIST)
      !SEPARATES STRING WITH SINGLE SEPARATOR=SEP        
      IMPLICIT NONE
      CHARACTER(LEN=*) :: STRING,SEP,STRLIST(2)
      INTEGER ISEPS,ISEPE

      ISEPS=INDEX(STRING,SEP)
      ISEPE=ISEPS+LEN(SEP)-1
      STRLIST(1)=TRIM(ADJUSTL(STRING(1:ISEPS-1)))
      STRLIST(2)=TRIM(ADJUSTL(STRING(ISEPE+1:)))
      
      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE REMOVE_DUPLICATE_STRINGS(N,STRL,STRLF,NF)
      IMPLICIT NONE
      INTEGER :: N,NF,I,J
      CHARACTER(LEN=*) :: STRL(N),STRLF(N)
      LOGICAL :: IS_STRING_PRESENT

      STRLF(1:N)=' ' 
      J=0
      DO I=1,N
       IF(.NOT.IS_STRING_PRESENT(N,STRLF,STRL(I))) THEN
        J=J+1
        STRLF(J)=STRL(I)
       ENDIF
      ENDDO
      NF=J

      RETURN
      END
      !-----------------------------------------------------------------
      FUNCTION IS_STRING_PRESENT(N,STRL,STRING)
      IMPLICIT NONE
      INTEGER :: N,I
      CHARACTER(LEN=*) :: STRL(N),STRING
      LOGICAL IS_STRING_PRESENT

      IS_STRING_PRESENT=.FALSE.
      DO I=1,N
       IF(STRING.EQ.STRL(I)) THEN
        IS_STRING_PRESENT=.TRUE.
        EXIT
       ENDIF
      ENDDO

      END FUNCTION
      !-----------------------------------------------------------------
