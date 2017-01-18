PWD := $(shell pwd)
OS_TYPE		:= $(shell echo `uname`)
MACHINE_TYPE	:= $(shell echo `uname -m`)
ifeq (${MACHINE_TYPE}, x86_64)
	MACHINE_TYPE = amd64
endif
BUILD_DIR	:= ${PWD}/${OS_TYPE}-${MACHINE_TYPE}/bin

all: dextract extractSequences mecatCanu mecat

aux_tools/hdf5/lib/libhdf5_cpp.so:
	cd aux_tools/hdf5 && make

dextract: aux_tools/hdf5/lib/libhdf5_cpp.so
	echo ${PWD}
	cd aux_tools/dextractor && make PATH_HDF5=${PWD}/aux_tools/hdf5 BUILD_DIR=${BUILD_DIR}

extractSequences: 
	cd extract_sequences && make

mecatCanu:
	cd mecat2canu/src && make

mecat:
	cd src && make

.PHONY: clean
clean:
	cd extract_sequences && make clean
	cd mecat2canu/src && make clean
	cd src && make clean
	
