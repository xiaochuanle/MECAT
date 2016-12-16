ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := mecat2cns
SOURCES  := main.cpp argument.cpp MECAT_AlnGraphBoost.C mecat_correction.cpp overlaps_partition.cpp reads_correction_aux.cpp reads_correction_can.cpp reads_correction_m4.cpp reads_correction_option.cpp

SRC_INCDIRS  := . libboost

TGT_LDFLAGS := -L${TARGET_DIR}
TGT_LDLIBS  := -lcanu
TGT_PREREQS := libcanu.a

SUBMAKEFILES :=
