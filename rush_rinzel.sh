#!/bin/bash
# $1 - Job Number
# $2 - Intervals
LOWMU=0
HIGHMU=0.3
LOWSIGMA=0
HIGHSIGMA=38
INDEX=$1
INTERVAL=`echo $(($2))`
INTERVALM1=`echo $(($2-1))`
# matlab linspace(d1,d2,n): y = [d1+(0:n-2)*(d2-d1)/(floor(n)-1) d2];
MUINDEX=`echo $(($INDEX / $INTERVAL))`
SIGMAINDEX=`echo $(($INDEX % $INTERVAL))`
MU=`echo $LOWMU + $MUINDEX*\($HIGHMU-\($LOWMU\)\)/$INTERVALM1 | bc -lq`
SIGMA=`echo $LOWSIGMA + $SIGMAINDEX*\($HIGHSIGMA-$LOWSIGMA\)/$INTERVALM1 | bc -lq`
FILE=`printf 'stats/%s_%s_mu=%.2f_sig=%.2f_c=%.2f_dt=%.3f_tmax=%.0f.dat' 'rr' 'stats' $MU $SIGMA .5 .01 800000`

mkdir -p code/corr/data
#echo "$MU $SIGMA" | ./rush_rinzel_simulate 0.5 0.01 800000 0 > code/corr/data/$1.simout
#./analysis $MU $SIGMA 0.5 0 0 0 0 0 0.01 800000 $1.simout > $FILE
echo would generate $FILE
# optionally, upload to traviscj.com's research database(password removed)
# cat $1.analysis | python stat_to_mysql.py $MU $SIGMA 0.5 | mysql -h traviscj.com -u research -pPASSWORD research
