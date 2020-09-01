Bootstrap: docker
From: ubuntu:bionic
Include: apt

%help
    This container is a portable installation of the FortBPSF program

%post
    echo "Installing gcc toolchain"
    apt update
    apt install -y build-essential gfortran libcurl4-gnutls-dev git 
    apt update
    apt install -y libhdf5-dev
    apt install -y rapidjson-dev

    cd /opt
    echo "Pulling fort-bpsf from Git"
    git clone https://github.com/lin-compchem/fort-bpsf.git
    echo "Making fort-bpsf"
    cd fort-bpsf/src/hdf5-1.12.0/
    bash fortbpsf-build.sh
    cd ../
    make

%runscript
    export LD_LIBRARY_PATH=/opt/fort-bpsf/lib/${LD_LIBRARY_PATH}
    exec /opt/fort-bpsf/bin/gen_symfuncs_parallel "$@"
