      SUBROUTINE GET_UNION_SET(NSPEC,NTRG,INDX_TRG,SPSET_TRGUP, &
                               SPSET_UNION)
      IMPLICIT NONE
      INTEGER :: I,J,NSPEC,NTRG,SPSET_UNION(NSPEC),INDX_TRG(NTRG), &
                 SPSET_TRGUP(NTRG,NSPEC)
     
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

      RETURN
      END
      !-----------------------------------------------------------------
