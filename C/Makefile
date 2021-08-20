# Makefile for MVGC mex files. Run 'make' from the command line
# in this directory. If build is successful, run 'mextest' from
# the Matlab prompt.
#
# Tested on gcc (Linux 64-bit) and mingw-w64 (Windows 64-bit);
# note that MSVC (still) does not fully support the C99
# standard, so you may have problems building with it. Should
# also build okay on MacOS with gcc (not sure about clang/LLVM).
#
# NOTE: use dynamic linking! (why wouldn't you...?)

ifdef ComSpec
	OS = WIN
else
	UNAME = $(shell uname)
	ifeq ($(UNAME),Linux)
		OS = LINUX
	else
		ifeq ($(UNAME),Darwin)
			OS = MAC
		else
			OS = UNKNOWN
		endif
	endif
endif

MEXEXT =  $(shell mexext)
MEXDIR = ../mex

MSRC   = findin_mex.c mvfilter_mex.c voulag_mex.c slidlyap_mex.c slidare_mex.c
# MSRC   = findin_mex.c mvfilter_mex.c voulag_mex.c
MBIN   = $(patsubst %.c,$(MEXDIR)/%.$(MEXEXT),$(MSRC))

MFLAGS  = -std=c99 -march=native
WFLAGS  = -Wall -Wextra -Wconversion -Winline -pedantic-errors
DFLAGS  = -D_POSIX_C_SOURCE=199309L -D_DEFAULT_SOURCE -D$(OS)
LDFLAGS = -lmwblas -lmwlapack -lmwslicot
ifeq ($(OS),WIN)
	OFLAGS    = -O3
#	MATLAB () say to use 'CFLAGS' for MinGW, but 'COMPFLAGS' seems to work
#   See: https://uk.mathworks.com/help/matlab/matlab_external/compiling-c-mex-files-with-mingw.html
	MEXCFLAGS = COMPFLAGS
	RM        = del /F /Q
	RMBIN     = $(subst /,\,$(MBIN))
else
	OFLAGS    = -O3 -flto
	MEXCFLAGS = CFLAGS
	RM        = rm -f
	RMBIN     = $(MBIN)
endif
GCCFLAGS = $(MFLAGS) $(OFLAGS) $(WFLAGS) $(DFLAGS)
MEXFLAGS = -O -R2018a

CC  = gcc
MEX = mex

.PHONY: all clean diag

all: $(MBIN)

clean:
	$(RM) $(RMBIN)

diag:
	@echo "*** MSRC   = " $(MSRC)
	@echo "*** MBIN   = " $(MBIN)
	@echo "*** MEXEXT = " $(MEXEXT)

$(MEXDIR)/%.$(MEXEXT): %.c
	@echo
	@echo Building $@
	$(MEX) $(MEXFLAGS) $(MEXCFLAGS)="\$$CFLAGS $(GCCFLAGS)" $(LDFLAGS) -outdir $(MEXDIR) $<
	@echo If build succeeded, you should run \'mextest\' in Matlab
