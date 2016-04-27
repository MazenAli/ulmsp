
# --- System configuration file for Linux with gcc

AR = ar
ARFLAGS = cru
RM = rm
CC = gcc
GCC = gcc

CFLAGS = -Wall -Wfatal-errors -Wpedantic -Wextra
LDFLAGS =
LIBS =

#
# Index options
#
ifdef BASE
  CFLAGS += -DINDEX_BASE=$(BASE)
endif

#
# Verbose setting
#
ifdef VERBOSE
  CFLAGS += -DVERBOSE
endif

#
# Debug settings
#
ifdef DEBUG
  CFLAGS += -g
endif

#
# Optimization settings
#
ifdef OPT
  CFLAGS += -O3 -march=native -funroll-loops -funswitch-loops -DNDEBUG
endif

#
# BLAS
#
ifdef USE_BLAS
  CFLAGS += -DUSE_BLAS -DADD_LIBS += -lblas -lgfortran
  LIBS += -lblas -lgfortran
endif

ifdef USE_CBLAS
  CFLAGS += -DUSE_CBLAS
endif

ifdef USE_OPENBLAS
  CFLAGS += -DUSE_OPENBLAS
  LIBS += -lopenblas
endif

#
# Math lib
#
LIBS += -lm
