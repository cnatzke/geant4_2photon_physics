#!/bin/bash

geant_install_dir="/data_fast/cnatzke/simulations/two-photon-emission/development/geant4/src/geant4-v10.7.4"

rsync -avzh source/ ${geant_install_dir}/source/
