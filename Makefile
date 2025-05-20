FC      = gfortran
FFLAGS  = -ffree-line-length-none -fbacktrace -Wall -Wextra -O2
LDLIBS  = -llapack -lblas

SRC     = numerics.f90 source.f90 fdtd.f90 main.f90
OBJDIR  = obj
MODDIR  = mod
BINDIR  = bin
DATADIR = data
OBJ     = $(patsubst %.f90,$(OBJDIR)/%.o,$(SRC))


$(OBJDIR)/%.o: %.f90 | $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -J$(MODDIR) -c $< -o $@


exec: $(OBJ) | $(BINDIR)
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ $^ $(LDLIBS)

all : $(DATADIR) exec

all: $(BINDIR) exec


clean:
	rm -f $(OBJDIR)/*.o $(MODDIR)/*.mod $(BINDIR)/exec $(DATADIR)/* .txt
	rm -f frames/*

.PHONY: all clean