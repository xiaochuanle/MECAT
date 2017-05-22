ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := filter_reads
SOURCES  := filter_reads.cpp 

SRC_INCDIRS  := ../common .

TGT_LDFLAGS := -L${TARGET_DIR}
TGT_LDLIBS  := -lmecat
TGT_PREREQS := libmecat.a

SUBMAKEFILES :=
