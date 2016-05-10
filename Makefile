
# ------------------------------------------------------------
# Components of the main library
# ------------------------------------------------------------

ULMSP_STORAGE = \
	Lib/storage/coo.c \
	Lib/storage/crs.c \
	Lib/storage/realvector.c \
	Lib/storage/indexvector.c \
	Lib/storage/realmatrix.c \
	Lib/storage/indexmatrix.c \

ULMSP_OPS= \
	Lib/ops/gecoomv.c \
	Lib/ops/gecrsmv.c \

ULMSP_SOLVERS= \
	Lib/solvers/cgcrs.c \
	Lib/solvers/cgcoo.c \
	Lib/solvers/gausscrs.c \

ULMSP_MESH= \
	Lib/mesh/mesh.c \

SOURCES_libulmsp := \
	$(ULMSP_STORAGE) \
	$(ULMSP_OPS) \
	$(ULMSP_SOLVERS) \
	$(ULMSP_MESH) \

HEADERS_libulmsp := $(SOURCES_libulmsp:.c=.h)

OBJECTS_libulmsp := $(SOURCES_libulmsp:.c=.o)

DEPENDENCIES_libulmsp := $(SOURCES_libulmsp:.c=.d)

# ------------------------------------------------------------
# Test programs
# ------------------------------------------------------------

SOURCES_stable := \
	Tests/setup_coo.c \
	Tests/setup_crs.c \
	Tests/setup_realvector.c \
	Tests/laplace_crs.c \
	Tests/coo2crs.c \
	Tests/testSolveGauss.c \
	Tests/refinetria.c \
	Tests/solveLaplace.c \
	Tests/solveLaplaceL.c \
	Tests/solveLaplace_CG_BPX.c \

SOURCES_tests = $(SOURCES_stable) \

OBJECTS_tests := \
	$(SOURCES_tests:.c=.o)

DEPENDENCIES_tests := \
	$(SOURCES_tests:.c=.d)

PROGRAMS_tests := \
	$(SOURCES_tests:.c=)

# ------------------------------------------------------------
# All files
# ------------------------------------------------------------

SOURCES := \
	$(SOURCES_libulmsp) \
	$(SOURCES_tests)

HEADERS := \
	$(HEADER_libulmsp)

OBJECTS := \
	$(OBJECTS_libulmsp) \
	$(OBJECTS_tests)

DEPENDENCIES := \
	$(DEPENDENCIES_libulmsp) \
	$(DEPENDENCIES_tests)

PROGRAMS := \
	$(PROGRAMS_tests)

# ------------------------------------------------------------
# Standard target
# ------------------------------------------------------------

all: programs

# ------------------------------------------------------------
# Build configuration
# ------------------------------------------------------------

$(OBJECTS): options.inc
include options.inc

# ------------------------------------------------------------
# System-dependent parameters (e.g., name of compiler)
# ------------------------------------------------------------

$(OBJECTS): system.inc
include system.inc

# ------------------------------------------------------------
# System-independent configuration (e.g., variants of algorithms)
# ------------------------------------------------------------


# ------------------------------------------------------------
# Rules for test programs
# ------------------------------------------------------------

programs: $(PROGRAMS_tests)

$(PROGRAMS_tests): %: %.o
ifdef BRIEF_OUTPUT
	@echo Linking $@
	@$(CC) $(LDFLAGS) -Wl,-L,.,-rpath,. $< -o $@.out -lulmsp $(LIBS)
else
	$(CC) $(LDFLAGS) -Wl,-L,.,-rpath,. $< -o $@.out -lulmsp $(LIBS)
endif

$(PROGRAMS_tests) $(PROGRAMS_tools): libulmsp.a

$(OBJECTS_tests): %.o: %.c
ifdef BRIEF_OUTPUT
	@echo Compiling $<
	@$(GCC) -MT $@ -MM -I Lib -I Lib/storage -I Lib/ops -I Lib/solvers \
			-I Lib/mesh \
			$< > $(<:%.c=%.d)
	@$(CC) $(CFLAGS) -I Lib -I Lib/storage -I Lib/ops -I Lib/solvers \
			-I Lib/mesh \
			-c $< -o $@
else
	@$(GCC) -MT $@ -MM -I Lib -I Lib/storage -I Lib/ops -I Lib/solvers \
			-I Lib/mesh \
			$< > $(<:%.c=%.d)
	$(CC) $(CFLAGS) -I Lib -I Lib/storage -I Lib/ops -I Lib/solvers \
			-I Lib/mesh \
			-c $< -o $@
endif

-include $(DEPENDENCIES_tests) $(DEPENDENCIES_tools)
$(OBJECTS_tests): Makefile

# ------------------------------------------------------------
# Rules for the Doxygen documentation
# ------------------------------------------------------------

doc:
	doxygen Doc/doxyfile

# ------------------------------------------------------------
# Rules for the main library
# ------------------------------------------------------------

libulmsp.a: $(OBJECTS_libulmsp)
ifdef BRIEF_OUTPUT
	@echo Building $@
	@$(AR) $(ARFLAGS) $@ $(OBJECTS_libulmsp)
else
	$(AR) $(ARFLAGS) $@ $(OBJECTS_libulmsp)
endif

$(OBJECTS_libulmsp): %.o: %.c
ifdef BRIEF_OUTPUT
	@echo Compiling $<
	@$(GCC) -MT $@ -MM -I Lib -I Lib/storage -I Lib/ops -I Lib/solvers \
			-I Lib/mesh \
			$< > $(<:%.c=%.d)
	@$(CC) $(CFLAGS) -I Lib -I Lib/storage -I Lib/ops -I Lib/solvers \
			-I Lib/mesh \
			-c $< -o $@
else
	@$(GCC) -MT $@ -MM -I Lib -I Lib/storage -I Lib/ops -I Lib/solvers \
			-I Lib/mesh \
			$< > $(<:%.c=%.d)
	$(CC) $(CFLAGS) -I Lib -I Lib/storage -I Lib/ops -I Lib/solvers \
			-I Lib/mesh \
			-c $< -o $@
endif

-include $(DEPENDENCIES_libulmsp)
$(OBJECTS_libulmsp): Makefile

# ------------------------------------------------------------
# Useful additions
# ------------------------------------------------------------

.PHONY: clean cleandoc programs indent

clean:
	$(RM) -f $(OBJECTS) $(DEPENDENCIES) $(PROGRAMS) libulmsp.a

cleandoc:
	$(RM) -rf Doc/html Doc/latex
