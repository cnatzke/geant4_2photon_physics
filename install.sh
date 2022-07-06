#!/bin/bash

geant_install_dir="/data/cnatzke/temp/cern/geant4.10.07.p03"

rsync -avzh source/ ${geant_install_dir}/source/
