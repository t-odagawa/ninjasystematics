#!/bin/sh

cd ${HOME}/NinjaSystematics/build/src

#nominalfilename=/hsm/nu/ninja/pra_tmp/CC0pi_20221213/detsys/Nominal/output/output_mode$1.root
nominalfilename=/hsm/nu/ninja/pra_tmp/CC0pi_20221213/signal/output/output_mode$1.root # Hit threshold
#nominalfilename=/hsm/nu/ninja/pra_tmp/CC0pi_20221213/detsys/MPPCNoise/nominal/output/output_mode$1.root # MPPC noise

dial=$2

detdirname=/hsm/nu/ninja/pra_tmp/CC0pi_20221213/detsys/${dial}

histname=$3

if [ ! -d /hsm/nu/ninja/pra_tmp/CC0pi_20221213/detsys/covariance_matrix/${dial} ];then
    mkdir -p /hsm/nu/ninja/pra_tmp/CC0pi_20221213/detsys/covariance_matrix/${dial}
fi

outputfilename=/hsm/nu/ninja/pra_tmp/CC0pi_20221213/detsys/covariance_matrix/${dial}/${histname}_covmat.root

./DetectorCovariance ${nominalfilename} ${detdirname} ${dial} ${histname} ${outputfilename}
