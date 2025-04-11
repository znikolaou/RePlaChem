      !-----------------------------------------------------------------
      !
      ! AUTHOR: Z. NIKOLAOU
      !
      !-----------------------------------------------------------------
      MODULE GLOBAL
      CHARACTER(LEN=*), PARAMETER :: &
                        INDIR='./input/',   &
                        OUTDIR='./output/', &
                        RATE_DIR='rates/',  &
                        CONTROL_FL='control.txt',      &
                        CHEMRED_FL='reducedChemistry.txt', &
                        RATE_FL='speciesMeanRatesMatrix.dat', &
                        REAC_RATE_FL='reactionRates'

      INTEGER, PARAMETER :: NREMX=6000,    &
                            NSPMX=300,     &
                            NSMX=500,      &
                            NLINEMX=10000, &
                            NSFLMX=1000 
      DOUBLE PRECISION, PARAMETER :: ZERO=0.0E0, &
                                     ONE=1.0E0 
      !                             
      INTEGER :: NELEM,NSPEC,NREAC,NSPEC_BOLSIG,NBOLS_SET,NREAC_DOLLAR
      CHARACTER(LEN=NSMX) :: ELEM(NSPMX),SPEC(NSPMX),REAC(NREMX), &
                             SPEC_BOLSIG(NSPMX), REAC_CONST(NREMX), &
                             BOLSIG_SEC_SET_LIST(NLINEMX), &
                             REAC_SEC_DOLLAR_LIST(NLINEMX), &
                             THIRD_BODY_SPEC(NREMX,NSPMX)
      LOGICAL :: IS_SPEC_CHARGED(NSPMX), &
                 IS_BOLSIG_REAC(NREMX), &
                 IS_ANY_NEUTRAL_REAC(NREMX), &
                 IS_ANY_ION_POS_REAC(NREMX), &
                 IS_ANY_ION_NEG_REAC(NREMX), &
                 IS_ANY_SPEC_REAC(NREMX), &
                 IS_THIRD_BODY_REAC(NREMX)
      INTEGER :: RSPEC(NREMX,NSPMX)
      DOUBLE PRECISION :: NUR(NREMX,NSPMX),NUP(NREMX,NSPMX), &
                          DELTANU(NREMX,NSPMX),SPEC_CHARGE(NSPMX)
      END MODULE GLOBAL
      !-----------------------------------------------------------------
