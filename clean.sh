#!/bin/bash

# build pw
cd mecat2pw
make clean
cd ..

# build ref
cd mecat2ref
make clean
cd ..

# build cns
cd mecat2cns
make clean
cd ..

# build canu
cd mecat2canu/src
make clean
cd ../..

# build extract sequences
cd extract_sequences
make clean
cd ..
