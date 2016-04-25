
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
