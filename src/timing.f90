      !-----------------------------------------------------------------
      !
      ! AUTHOR: Z. NIKOLAOU
      !
      !-----------------------------------------------------------------
      SUBROUTINE DATETIME(DATE,TIME,ZONE)
      IMPLICIT NONE
      CHARACTER(LEN=8)  :: DATE
      CHARACTER(LEN=10) :: TIME
      CHARACTER(LEN=5)  :: ZONE
      INTEGER :: VALUES(8)

      CALL DATE_AND_TIME(DATE,TIME,ZONE,VALUES)
            
      RETURN
      END
      !-----------------------------------------------------------------




