ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := mecat2ref
SOURCES  := mecat2ref.cpp mecat2ref_impl_large.cpp output.cpp mecat2ref_aux.cpp

SRC_INCDIRS  := . 

TGT_LDFLAGS := -L${TARGET_DIR}
TGT_LDLIBS  := -lmecat
TGT_PREREQS := libmecat.a

SUBMAKEFILES :=
