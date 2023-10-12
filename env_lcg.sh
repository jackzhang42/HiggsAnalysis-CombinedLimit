. /cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-centos7-gcc11-opt/setup.sh
# . /cvmfs/sft.cern.ch/lcg/views/LCG_98/x86_64-centos7-gcc9-opt/setup.sh
export PATH=${PATH}:${PWD}/build/bin
export LD_LIBRARY_PATH=${PWD}/build/lib:${LD_LIBRARY_PATH}
# export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PWD}/build/lib
export PYTHONPATH=${PYTHONPATH}:${PWD}/build/lib/python

# dir="/afs/cern.ch/user/y/yangfan/eos/analysis/projects/couplings/cms"
# export PATH=${PATH}:${dir}/build/bin
# export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${dir}/build/lib
# export PYTHONPATH=${PYTHONPATH}:${dir}/build/lib/python
