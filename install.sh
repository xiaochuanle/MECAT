#!/bin/bash

cd src
make -j4
cd ..

cd mecat2canu/src
make -j4
cd ../..

cd extract_sequences
make -j4
cd ..

