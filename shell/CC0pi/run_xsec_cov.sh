#!/bin/sh

cd ${HOME}/NinjaSystematics/build/src

nominalfilename=/hsm/nu/ninja/pra_tmp/CC0pi_20221213/signal/output/output_mode$1.root

xsecdirname=/hsm/nu/ninja/pra_tmp/CC0pi_20221213/xsecsys

dial=/$2

histname=$3

if [ ! -d /hsm/nu/ninja/pra_tmp/CC0pi_20221213/xsecsys/covariance_matrix/${dial} ];then
    mkdir -p /hsm/nu/ninja/pra_tmp/CC0pi_20221213/xsecsys/covariance_matrix/${dial}
fi

outputfilename=/hsm/nu/ninja/pra_tmp/CC0pi_20221213/xsecsys/covariance_matrix/${dial}/${histname}_covmat.root

./InteractionModelCovariance ${nominalfilename} ${xsecdirname} ${dial} ${histname} ${outputfilename}
