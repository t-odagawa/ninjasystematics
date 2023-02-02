#!/bin/sh

cd ${HOME}/NinjaSystematics/build/src

nominalfilename=/hsm/nu/ninja/pra_tmp/detsys/Nominal/output/output_mode$1.root
#nominalfilename=/hsm/nu/ninja/pra_tmp/detsys/MPPCNoise/nominal/output/output_mode$1.root # MPPC noise
#nominalfilename=/group/nu/ninja/work/odagawa/20221020-phd-thesis-preliminary/signal/output/output_mode$1.root # Hit threshold

dial=$2

detdirname=/hsm/nu/ninja/pra_tmp/detsys/${dial}

histname=$3

outputfilename=/hsm/nu/ninja/pra_tmp/detsys/covariance_matrix/${dial}/${histname}_covmat.root

./DetectorCovariance ${nominalfilename} ${detdirname} ${dial} ${histname} ${outputfilename}
