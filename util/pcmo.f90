      !-----------------------------------------------------------------
      !
      ! AUTHOR: Z. NIKOLAOU
      !
      !-----------------------------------------------------------------
      MODULE PCMO

      CONTAINS

      !----------------------------------------------------------------- 
      SUBROUTINE PCMO_WRITE_BIN(DIR,FLNAME,ICOUNT,TSTEP,TIME,N,ARRAY)
      IMPLICIT NONE
      INTEGER :: N,ICOUNT,TSTEP
      CHARACTER(LEN=*) :: DIR,FLNAME
      CHARACTER(LEN=LEN(DIR)+LEN(FLNAME)+14) :: FL
      DOUBLE PRECISION :: TIME,ARRAY(N)

      FL=BUILD_FLNAME(ICOUNT,DIR,FLNAME)
      OPEN(UNIT=1,FILE=FL,STATUS='REPLACE',FORM='UNFORMATTED')
      WRITE(1) ICOUNT,TSTEP,TIME,N,ARRAY
      CLOSE(1)
      
      RETURN
      END
      !-----------------------------------------------------------------
      FUNCTION BUILD_FLNAME(TSTEP,DIR,FLNAME)
      IMPLICIT NONE
      CHARACTER(LEN=*) :: DIR,FLNAME
      INTEGER :: TSTEP
      CHARACTER(LEN=8) :: CT
      CHARACTER(LEN=LEN(DIR)+LEN(FLNAME)+14) :: BUILD_FLNAME
    
      WRITE(CT,'(I8.8)') TSTEP
      BUILD_FLNAME=TRIM(ADJUSTL(DIR))//'/'// &
                   TRIM(ADJUSTL(FLNAME))//'_'//CT//'.dat'

      END FUNCTION
      !-----------------------------------------------------------------
      END MODULE PCMO
