      MODULE GLOBAL
      
      INTEGER, PARAMETER :: NREMX=10000,NSPMX=500,NSMX=500, &
                            NLINEMX=100000,NSFLMX=1000
      CHARACTER(LEN=*), PARAMETER :: INDIR='./input/', &
                                     OUTDIR='./output/', &
                                     CONTROL_FL='control.txt', &
                                     RATE_FL='rates.txt', &
                                     CHEMRED_FL='reducedChemistry.txt'
      INTEGER :: NREAC,NSPEC,NSPEC_BOLSIG,NELEM
      CHARACTER(LEN=NSMX) :: ELEM(NSPMX),SPEC(NSPMX),REAC(NREMX), &
                             SPEC_CHARGE(NSPMX),REACF(NREMX), &
                             SPEC_BOLSIG(NSPMX), REAC_CONST(NREMX), &
                             REAC_SPEC(NREMX,NSPMX)
      LOGICAL :: IS_SPEC_IN_REAC(NSPMX,NREMX),IS_SPEC_CHARGED(NSPMX), &
                 IS_BOLSIG_REAC(NREMX)
      INTEGER :: RSPEC(NREMX,NSPMX)

      END MODULE GLOBAL
