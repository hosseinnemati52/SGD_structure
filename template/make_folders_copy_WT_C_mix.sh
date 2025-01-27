#!/bin/bash

mkdir main
mkdir plus_dh
mkdir minus_dh

# Copy all contents from source into the run_$i folders
cp -r source_WT_C_mix/* "main"
cp -r source_WT_C_mix/* "plus_dh"
cp -r source_WT_C_mix/* "minus_dh"

