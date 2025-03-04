      SUBROUTINE READ_RATES(IDAT,DIR,RNM,NSPEC,NREAC,WIJ,RR,JIJ)
!
!     AUTHOR: ZACHARIAS M. NIKOLAOU.
!
!     DESCRIPTION: READS IN DATA.
!     NSP=NUMBER OF SPECIES.
!     NRE=NUMBER OF REACTIONS.
!     TIME=TIME (NOT USED).
!     THETA=SOLAR AZIMUTH ANGLE (NOT USED).
!     WIJ=SPECIES I RATE FROM REACTION J.
!     RR=RATE OF REACTION.
!     JIJ=JACOBIAN-BASED DIC VALUES AS CALCULATED IN
!     Y. Chen, J.Y Chen. Combust. Flame 174 (2016) 77-84. 
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
!
      USE GLOBAL, ONLY :NWRK
      USE PRECIS, ONLY : DBL_P
!
      IMPLICIT NONE
      INTEGER I,J,IDAT,NSPEC,NREAC
      CHARACTER(LEN=*) :: DIR,RNM
      REAL(KIND=DBL_P) :: WIJ(NREAC,NSPEC),RR(NREAC),JIJ(NSPEC,NSPEC)
      !
      CHARACTER(LEN=4) :: CDAT
      CHARACTER(LEN=NWRK) :: FLNM
      INTEGER :: NSP,NRE
      REAL(KIND=DBL_P) :: TIME,THETA
! 
      !BUILD FILENAMES
      WRITE(CDAT,'(I4.4)') IDAT     
      FLNM=DIR//RNM//'_'//CDAT//'.dat'

      OPEN(UNIT=1,FILE=FLNM,STATUS='OLD',FORM='UNFORMATTED')
      READ(1) NSP,NRE
      CLOSE(1)
      
      IF(NSP.NE.NSPEC.OR.NRE.NE.NREAC) THEN
       WRITE(*,*) 'READ_RATES:ERROR, MISMATCH',NSP,NSPEC,NRE,NREAC
       STOP
      ELSE
       OPEN(UNIT=1,FILE=FLNM,STATUS='OLD',FORM='UNFORMATTED')
       READ(1) NSP,NRE,TIME,THETA,WIJ,RR,JIJ
       CLOSE(1)
      ENDIF
!           
      END SUBROUTINE
