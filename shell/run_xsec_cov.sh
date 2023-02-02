#!/bin/sh

cd ${HOME}/NinjaSystematics/build/src

nominalfilename=/group/nu/ninja/work/odagawa/20221020-phd-thesis-preliminary/signal/output/output_mode$1.root

xsecdirname=/hsm/nu/ninja/pra_tmp/xsecsys

dial=/$2

histname=$3

outputfilename=/hsm/nu/ninja/pra_tmp/xsecsys/covariance_matrix/${dial}/${histname}_covmat.root

./InteractionModelCovariance ${nominalfilename} ${xsecdirname} ${dial} ${histname} ${outputfilename}
