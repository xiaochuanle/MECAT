ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := mecat2cns
SOURCES  := main.cpp \
	dw.cpp \
	MECAT_AlnGraphBoost.C \
	mecat_correction.cpp \
	options.cpp \
	overlaps_partition.cpp \
	reads_correction_aux.cpp \
	reads_correction_can.cpp \
	reads_correction_m4.cpp \
	packed_db.cpp \

SRC_INCDIRS  := . libboost

# make sure large files are okay (and off_t is 8 bytes);
# requires c++11 or higher for <mutex> headers
TGT_CXXFLAGS := -static -D _FILE_OFFSET_BITS=64 -std=c++11 -O3 -mfpmath=sse -march=native -flto -fno-fat-lto-objects -fno-builtin -mmmx -msse -msse2 -mssse3 -msse4.1 -msse4.2 -mavx -maes -mpopcnt -mfxsr -mxsave -mxsaveopt
TGT_LDFLAGS  := -static -L${TARGET_DIR} ${TGT_CXXFLAGS}
TGT_LDLIBS   := -lmecat
TGT_PREREQS  := libmecat.a

SUBMAKEFILES :=
