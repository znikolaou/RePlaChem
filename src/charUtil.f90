      !-----------------------------------------------------------------
      LOGICAL FUNCTION IS_DOT(A)
      IMPLICIT NONE
      CHARACTER(LEN=*) :: A
      CHARACTER, PARAMETER :: CDOT='.'

      IS_DOT=A.EQ.CDOT

      END FUNCTION
      !-----------------------------------------------------------------
      LOGICAL FUNCTION IS_NUMBER(A)
      IMPLICIT NONE
      CHARACTER :: A
      INTEGER, PARAMETER :: ASCI_IL=48,ASCI_IH=57
      INTEGER :: IA

      IA=ICHAR(A)
      IS_NUMBER=.FALSE.
      IF((IA.GE.ASCI_IL).AND.(IA.LE.ASCI_IH)) THEN
       IS_NUMBER=.TRUE.
      ENDIF

      END FUNCTION
      !-----------------------------------------------------------------
      LOGICAL FUNCTION IS_AZ(A)
      IMPLICIT NONE
      LOGICAL :: IS_NUMBER,IS_DOT
      CHARACTER :: A
      INTEGER, PARAMETER :: ASCI_IL=65,ASCI_IH=90
      INTEGER :: IA

      IS_AZ=(.NOT.IS_NUMBER(A)).AND.(.NOT.IS_DOT(A))

      END FUNCTION
      !-----------------------------------------------------------------
      LOGICAL FUNCTION ISEMPTY(STRING)
      !
      !AUTHOR: Z. NIKOLAOU
      !
      IMPLICIT NONE
      CHARACTER(LEN=*) STRING
       
      ISEMPTY=STRING.EQ.' '

      END FUNCTION
      !-----------------------------------------------------------------
      LOGICAL FUNCTION ISCOMMENT(STRING)
      !
      !AUTHOR: Z. NIKOLAOU
      !    
      IMPLICIT NONE
      CHARACTER(LEN=*) STRING

      ISCOMMENT=STRING(1:1).EQ.'!'.OR.STRING(1:1).EQ.'#'
      ISCOMMENT=ISCOMMENT.OR.STRING(1:1).EQ.'$'

      END FUNCTION
      !-----------------------------------------------------------------
      LOGICAL FUNCTION IS_STRING_PRESENT(N,STRL,STRING)
      !
      !AUTHOR: Z. NIKOLAOU
      !      
      IMPLICIT NONE
      INTEGER :: N,I
      CHARACTER(LEN=*) :: STRL(N),STRING
   
      IS_STRING_PRESENT=.FALSE.
      DO I=1,N
       IF(STRING.EQ.STRL(I)) THEN
        IS_STRING_PRESENT=.TRUE.
        EXIT
       ENDIF
      ENDDO

      END FUNCTION  
      !-----------------------------------------------------------------
      SUBROUTINE SPLIT_NUM_AND_CHAR(STR,NS,NUM,CH)
      IMPLICIT NONE
      !
      !SPLITS A STRING OF THE FORM '123.313AFD' INTO 123.313 AND AFD. 
      !
      INTEGER I,IAZ,NS,GET_INDEX_FIRST_AZ
      CHARACTER(LEN=*) :: STR
      CHARACTER(LEN=NS) :: CH
      DOUBLE PRECISION :: NUM

      NUM=1.0e0
      CH=' '
      IAZ=GET_INDEX_FIRST_AZ(STR)
      IF(IAZ.EQ.0) THEN         
       READ(STR,*) NUM
      ELSEIF(IAZ.EQ.1) THEN
       CH=STR
      ELSE
       READ(STR(1:IAZ-1),*) NUM
       CH=STR(IAZ:)
      ENDIF

      RETURN
      END
      !-----------------------------------------------------------------
      FUNCTION GET_KEY_INDEX(KEY,NL,STRL,IS)
      !
      !GET INDEX OF (UNIQUE) KEY PRESENT IN LIST STRL
      !
      !AUTHOR: Z. NIKOLAOU
      !      
      IMPLICIT NONE
      CHARACTER(LEN=*) :: KEY,STRL(NL)
      INTEGER :: GET_KEY_INDEX,I,NL,IS

      GET_KEY_INDEX=0
      DO I=IS,NL
       IF(TRIM(ADJUSTL(STRL(I))).EQ.TRIM(ADJUSTL(KEY))) THEN
        GET_KEY_INDEX=I
        EXIT
       ENDIF
      ENDDO
    
      END FUNCTION
      !-----------------------------------------------------------------
      FUNCTION GET_INDEX_FIRST_AZ(STR)
      IMPLICIT NONE
      LOGICAL :: IS_AZ
      INTEGER :: I,GET_INDEX_FIRST_AZ
      CHARACTER(LEN=*) :: STR

      GET_INDEX_FIRST_AZ=0
      DO I=1,LEN(STR)
       IF(IS_AZ(STR(I:I))) THEN
        GET_INDEX_FIRST_AZ=I
        EXIT
       ENDIF
      ENDDO

      END FUNCTION
      !-----------------------------------------------------------------
      FUNCTION GET_INDEX_FIRST_CHAR(STR)
      !
      !AUTHOR: Z. NIKOLAOU
      !      
      IMPLICIT NONE
      CHARACTER(LEN=*) :: STR
      INTEGER :: GET_INDEX_FIRST_CHAR,I

      GET_INDEX_FIRST_CHAR=0
      DO I=1,LEN(STR)
       IF(STR(I:I).NE.' ') THEN
        GET_INDEX_FIRST_CHAR=I
        EXIT
       ENDIF
      ENDDO
      
      END FUNCTION
      !-----------------------------------------------------------------
      SUBROUTINE GET_INDEX(STR,TXT,N,IARR,NA)
      !         
      !GET ALL STARTING INDICES OF TXT IN STR
      ! 
      !AUTHOR: Z. NIKOLAOU
      !      
      IMPLICIT NONE
      CHARACTER(LEN=*) :: STR,TXT
      INTEGER :: N,NA,IARR(N,2),IS,IP,IE,J,LT

      IS=INDEX(STR,TXT)
      IP=IS
      LT=LEN(TXT)
      IE=IS+LT-1
      J=0
      DO WHILE(IP.NE.0)
       J=J+1
       IARR(J,1)=IS
       IARR(J,2)=IE
       IP=INDEX(STR(IE+1:),TXT)
       IS=IP+IE
       IE=IS+LT-1
      ENDDO
      NA=J

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE GET_SUBTEXT_LIST(STR,TXT,N,SLIST,NS)
      !
      !AUTHOR: Z. NIKOLAOU
      !
      IMPLICIT NONE
      INTEGER :: N,NS,I,IAR(N,2)
      CHARACTER(LEN=*) :: STR,TXT,SLIST(N)

      SLIST(1:N)=' '
      CALL GET_INDEX(STR,TXT,N,IAR,NS)
      DO I=1,NS
       SLIST(I)=STR(IAR(I,1):IAR(I,2))
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE GET_SUBTEXT_DISTINCT_LIST(STR,TXT,N,SLIST,ND)
      !
      !AUTHOR: Z. NIKOLAOU
      !
      IMPLICIT NONE
      INTEGER :: N,ND,I,LTXT,NA,IAR(N,2)
      CHARACTER(LEN=*) :: STR,TXT,SLIST(N)
      LOGICAL :: IS_STRING_PRESENT

      SLIST(1:N)=' '
      CALL GET_INDEX(STR,TXT,N,IAR,NA)
      ND=0
      DO I=1,NA
      IF(.NOT.IS_STRING_PRESENT(N,SLIST,STR(IAR(I,1):IAR(I,2)))) THEN
       SLIST(I)=STR(IAR(I,1):IAR(I,2))
       ND=ND+1
      ENDIF
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE GET_SUBTEXTP1_DISTINCT_LIST(STR,TXT,N,SLIST,ND)
      !
      !AUTHOR: Z. NIKOLAOU
      !
      IMPLICIT NONE
      INTEGER :: N,ND,I,LTXT,NA,IAR(N,2)
      CHARACTER(LEN=*) :: STR,TXT,SLIST(N)
      LOGICAL :: IS_STRING_PRESENT

      SLIST(1:N)=' '
      CALL GET_INDEX(STR,TXT,N,IAR,NA)
      ND=0
      DO I=1,NA
      IF(.NOT.IS_STRING_PRESENT(N,SLIST,STR(IAR(I,1):IAR(I,2)+1))) THEN
       SLIST(I)=STR(IAR(I,1):IAR(I,2)+1)
       ND=ND+1
      ENDIF
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
      FUNCTION GET_INDEX_FIRST_SPACE(STR)
      ! 
      !AUTHOR: Z. NIKOLAOU
      ! 
      IMPLICIT NONE
      CHARACTER(LEN=*) :: STR
      CHARACTER(LEN=*), PARAMETER :: SPACE=' '
      INTEGER :: GET_INDEX_FIRST_SPACE,I
      
      GET_INDEX_FIRST_SPACE=0
      DO I=1,LEN(STR)
       IF(STR(I:I).EQ.SPACE) THEN
        GET_INDEX_FIRST_SPACE=I
        EXIT
       ENDIF
      ENDDO

      END FUNCTION
      !-----------------------------------------------------------------
      FUNCTION GET_TEXT_WITH_SPACES_COUNT(STRING,TEXT)
      !
      ! AUTHOR: Z. NIKOLAOU     
      !
      ! DESCRIPTION: COUNTS NO OF OCCURANCES OF TEXT IN 
      ! STRING, OF THE FORM TEXT' ' AND ' 'TEXT' '  
      !
      IMPLICIT NONE
      CHARACTER(LEN=*) :: STRING,TEXT
      INTEGER :: I,NSTR,NTXTLR,NTXTR,GET_TEXT_WITH_SPACES_COUNT

      NSTR=LEN(STRING)
      NTXTLR=LEN(TRIM(TEXT))+2
      NTXTR=LEN(TRIM(TEXT))+1
      
      !CHECK FOR 1ST ELEMENT TEXT' '
      GET_TEXT_WITH_SPACES_COUNT=0
      IF( STRING(1:1+NTXTR-1).EQ.(TEXT(1:NTXTLR-2)//' ') ) THEN
       GET_TEXT_WITH_SPACES_COUNT=GET_TEXT_WITH_SPACES_COUNT+1
      ENDIF
      !CHECK FOR ' 'TEXT' '
      DO I=2,NSTR-NTXTLR+1
       IF( (STRING(I:I+NTXTLR-1)).EQ.(' '//TEXT(1:NTXTLR-2)//' ') ) THEN 
        GET_TEXT_WITH_SPACES_COUNT=GET_TEXT_WITH_SPACES_COUNT+1
       ENDIF
      ENDDO

      END FUNCTION 
      !----------------------------------------------------------------- 
      FUNCTION GET_TEXT_COUNT(STRING,TEXT)
      !
      ! AUTHOR: Z. NIKOLAOU     
      !
      IMPLICIT NONE
      CHARACTER(LEN=*) :: STRING,TEXT
      INTEGER :: I,NSTR,NTXT,GET_TEXT_COUNT

      NSTR=LEN(TRIM(STRING))
      NTXT=LEN(TRIM(TEXT))
      
      GET_TEXT_COUNT=0
      DO I=1,NSTR-NTXT+1
       IF( (STRING(I:I+NTXT-1)).EQ.(TEXT) ) THEN
        GET_TEXT_COUNT=GET_TEXT_COUNT+1
       ENDIF
      ENDDO
      
      END FUNCTION
      !-----------------------------------------------------------------
      FUNCTION STR2NUM(C)
      !
      ! AUTHOR: Z. NIKOLAOU
      !
      IMPLICIT NONE
      INTEGER :: LENC
      CHARACTER(LEN=*) :: C
      !REAL(KIND=DBL_P) :: CHAR2R
      REAL :: STR2NUM,A

      READ(C,'(G12.5)') A
     
      STR2NUM=A                

      END FUNCTION
      !-----------------------------------------------------------------
      SUBROUTINE REMOVE_SPACES(L,STRING,LS) 
      !
      ! AUTHOR: Z. NIKOLAOU
      !
      IMPLICIT NONE
      INTEGER :: L,LS,I,J
      CHARACTER(LEN=L) :: STRING,STRWRK
       
      J=0
      STRWRK=' ' 
      DO I=1,L
       IF(STRING(I:I).NE.' ') THEN
        J=J+1
        STRWRK(J:J)=STRING(I:I)
       ENDIF 
      ENDDO
      STRING(1:J)=STRWRK(1:J) 
      STRING(J+1:L)=' '
      LS=J

      END
      !-----------------------------------------------------------------
      SUBROUTINE REMOVE_DUPLICATE_STRINGS(N,STRL,STRLF,NF)
      !
      !AUTHOR: Z. NIKOLAOU
      !      
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
      SUBROUTINE REPLACE_TEXT_IN_STRING(STRING,TEXT,NEWTEXT)
      !
      ! TODO: SOME BUG HERE NEWTEXT OVERWRITES ORIGINAL STRING IF LENGTH NOT
      ! EQUAL TO TEXT!
      ! AUTHOR: Z. NIKOLAOU     
      ! 
      IMPLICIT NONE
      CHARACTER(LEN=*) :: STRING,TEXT,NEWTEXT
      INTEGER :: I,NSTR,NTXT,NNEWTXT

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

      END SUBROUTINE
      !-----------------------------------------------------------------
      SUBROUTINE REPLACE_TEXT(STRING,TEXT,NEWTEXT,STRWRK,LWRK)
      !
      ! AUTHOR: Z. NIKOLAOU     
      !
      IMPLICIT NONE
      INTEGER :: LSTR,LTXT,LTXTN,LWRK,I,IL,IR,IS,DIFF,NOCC,LENNEW
      CHARACTER(LEN=*) :: STRING,TEXT,NEWTEXT
      CHARACTER(LEN=LWRK) :: STRWRK
      
      LSTR=LEN(STRING)
      LTXT=LEN(TEXT)
      LTXTN=LEN(NEWTEXT)
      
      STRWRK(1:LWRK)=' '
      DIFF=LTXTN-LTXT       
      IS=LEN(TRIM(ADJUSTL(STRWRK))) 

      NOCC=0
      IL=1
      DO I=1,LSTR-LTXT+1
       IF( (STRING(I:I+LTXT-1)).EQ.(TEXT) ) THEN !TEXT IN STRING
        NOCC=NOCC+1
        IR=I-1
        STRWRK(IS+1:)=STRING(IL:IR)//NEWTEXT
        IS=IS+(IR-IL+1)+LTXTN!LEN(TRIM(ADJUSTL(STRWRK))) 
        IL=I+LTXT         
       ENDIF              
      ENDDO

      !NEW LENGTH
      LENNEW=LEN(TRIM(ADJUSTL(STRING)))+DIFF*NOCC
      IF(LENNEW.GT.LWRK) THEN
       WRITE(*,*) '*REPLACE_TEXT: NOT ENOUGH WORSPACE PROVIDED'
       STOP
      ENDIF

      !LAST PART
      STRWRK(IS+1:)=STRING(IL:LSTR)

      END
      !------------------------------------------------------------------
      SUBROUTINE SPLIT_STRING_WITH_SPACES(STRING,N,COLMS,NA)
      !
      !AUTHOR: Z. NIKOLAOU
      !       
      IMPLICIT NONE
      INTEGER :: N,NA,IS,IE,J,GET_INDEX_FIRST_CHAR, &
                 ISC,GET_INDEX_FIRST_SPACE
      CHARACTER(LEN=*) :: STRING,COLMS(N)
      CHARACTER(LEN=*), PARAMETER :: SPACE=' '
      LOGICAL :: ISEMPTY

      COLMS(1:N)=''
      NA=0
      IF(ISEMPTY(STRING)) THEN
       RETURN
      ELSE
       IF(INDEX(STRING,SPACE).EQ.0) THEN
        COLMS(1)=STRING(1:LEN(STRING))
        NA=1
        RETURN
       ENDIF
      ENDIF
     
      IS=1 
      J=0
      ISC=GET_INDEX_FIRST_CHAR(STRING(IS:))
      DO WHILE(ISC.GT.0)
       IS=ISC+IS-1
       IE=GET_INDEX_FIRST_SPACE(STRING(IS+1:))+IS
       IF(IE.EQ.IS) THEN
        J=J+1
        COLMS(J)=STRING(IS:LEN(STRING))
        EXIT
       ENDIF
       J=J+1
       COLMS(J)=STRING(IS:IE-1)
       IS=IE+1
       ISC=GET_INDEX_FIRST_CHAR(STRING(IS:))
      ENDDO
      NA=J
        
      END 
      !-----------------------------------------------------------------
      SUBROUTINE SPLIT_STRING(STRING,SEP,N,COLMS,NA)
      !
      !AUTHOR: Z. NIKOLAOU
      !      
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
      !
      !AUTHOR: Z. NIKOLAOU
      !      
      !SEPARATES STRING WITH SINGLE SEPARATOR=SEP    
      !    
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
      
