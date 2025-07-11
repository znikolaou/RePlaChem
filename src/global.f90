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
       CONTROL_FL='control.txt', &
       CHEMRED_FL='skeletalChemistry.txt', &
       REAC_RATE_FL='reactionRates', &
       STATS_FL='summary.txt', &
       VERSION='25.05', &
       LOGO=  '# RePlaChem_'//VERSION//'     #', &
       AUTHOR='# AUTHOR: Z. NIKOLAOU #', &
       DASH=  '#---------------------#'
      INTEGER, PARAMETER :: & 
       NREMX=20000, &
       NSPMX=5000, &
       NSMX=300, &
       NSMX_TB=30, &
       NLINEMX=80000, &
       NSFLMX=1200 
      DOUBLE PRECISION, PARAMETER :: ZERO=0.0E0, ONE=1.0E0 
      !
      INTEGER :: NELEM,NSPEC,NREAC,NREAC_BOLSIG,NSPEC_BOLSIG, &
                 NBOLS_SET,NTRG,NCASE,NDATA,NREAC_DOLLAR, &
                 RSPEC(NREMX,NSPMX),INDX_TRG(NSPMX)
      CHARACTER(LEN=NSFLMX) :: CHEMFL
      CHARACTER(LEN=NSMX) :: ELEM(NSPMX),SPEC(NSPMX),REAC(NREMX), &
       SPEC_BOLSIG(NSPMX),REAC_CONST(NREMX), &
       BOLSIG_SEC_SET_LIST(NLINEMX),REAC_SEC_DOLLAR_LIST(NLINEMX)
      CHARACTER(LEN=NSMX_TB) :: THIRD_BODY_SPEC(NREMX,NSPMX)
      LOGICAL :: IS_SPEC_CHARGED(NSPMX),IS_BOLSIG_REAC(NREMX), &
       IS_ANY_NEUTRAL_REAC(NREMX),IS_ANY_ION_POS_REAC(NREMX), &
       IS_ANY_ION_NEG_REAC(NREMX),IS_ANY_SPEC_REAC(NREMX), &
       IS_THIRD_BODY_REAC(NREMX)
      DOUBLE PRECISION :: NUR(NREMX,NSPMX),NUP(NREMX,NSPMX), &
       DELTANU(NREMX,NSPMX),SPEC_CHARGE(NSPMX),ETOL(NSPMX)
      END MODULE GLOBAL
      !-----------------------------------------------------------------
