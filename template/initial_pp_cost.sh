#!/bin/bash

currentdir=$(pwd)

cd "main/WT"
./pp_plus_cost.sh
cd $currentdir

cd "main/C"
./pp_plus_cost.sh
cd $currentdir

cd "main/mix"
./pp_plus_cost.sh
cd $currentdir
