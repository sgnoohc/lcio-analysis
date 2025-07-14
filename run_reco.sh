# git clone git@github.com:madbaron/MyBIBUtils.git
# cd MyBIBUtils/
# mkdir build
# cd build/
# cmake ..
# make install
# export MARLIN_DLL=/home/users/phchang/work/muc/workdir/MyBIBUtils/lib/libMyBIBUtils.so:$MARLIN_DLL

k4run \
    SteeringMacros/k4Reco/steer_reco.py \
    --TypeEvent muonGun_pT_0_50
