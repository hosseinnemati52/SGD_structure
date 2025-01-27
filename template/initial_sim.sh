#!/bin/bash

currentdir=$(pwd)

cd "main/WT"
./clear_copy_init_run_pp.sh
./terminal_opener.sh
cd $currentdir

cd "main/C"
./clear_copy_init_run_pp.sh
./terminal_opener.sh
cd $currentdir

cd "main/mix"
./clear_copy_init_run_pp.sh
./terminal_opener.sh
cd $currentdir
