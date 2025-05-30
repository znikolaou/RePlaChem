      MODULE PCMO
      !----------------------------------------------------------------- 
      SUBROUTINE PCMO_WRITE_BIN(DIR,FLNAME,ICOUNT,TSTEP,TIME,N,ARRAY)
      IMPLICIT NONE
      INTEGER, PARAMETER : IUNIT=1
      INTEGER :: N,ICOUNT,TSTEP
      CHARACTER(LEN=*) :: DIR,FLNAME
      DOUBLE PRECISION :: TIME,ARRAY(N)

      FL=BUILD_FLNAME(ICOUNT,DIR,FLNAME)
      OPEN(UNIT=IUNIT,FILE=FL,STATUS='REPLACE',FORM='UNFORMATTED')
      WRITE(IUNIT) ICOUNT,TSTEP,TIME,N,ARRAY
      CLOSE(IUNIT)
      
      RETURN
      END
      !-----------------------------------------------------------------
      FUNCTION BUILD_FLNAME(TSTEP,DIR,FLNAME)
      IMPLICIT NONE
      CHARACTER(LEN=*) :: DIR,FLNAME
      INTEGER :: TSTEP
      CHARACTER(LEN=8) :: CT
      CHARACTER(LEN=LEN(DIR)+LEN(FLNAME)+12) :: BUILD_FLNAME
    
      WRITE(CT,'(I8.8)') TSTEP
      BUILD_FLNAME=TRIM(ADJUSTL(DIR))//'/'// &
                   TRIM(ADJUSTL(FLNAME))//'_'//CT//'.dat'

      END FUNCTION
      !-----------------------------------------------------------------
      END MODULE PCMO
