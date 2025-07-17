# Configuration MUMPS
topdir = /home/emin/Documents/MUMPS_5.8.0
libdir = $(topdir)/lib
includedir = $(topdir)/include

# Inclusion du makefile MUMPS qui définit les variables essentielles
include $(topdir)/Makefile.inc

# Structure des répertoires
SRC     = numerics.f90 source.f90 fdtd.f90 main.f90
OBJDIR  = obj
MODDIR  = mod
BINDIR  = bin
DATADIR = data
OBJ     = $(patsubst %.f90,$(OBJDIR)/%.o,$(SRC))

# Options de compilation supplémentaires
FFLAGS  += -ffree-line-length-none -fbacktrace -Wall -Wextra -O2 -I$(includedir) -fcheck=all

# Définition des bibliothèques MUMPS 
LIBSDMUMPS = -L$(libdir) -lsmumps$(PLAT) -ldmumps$(PLAT) -lmumps_common$(PLAT)

# Règles de compilation
$(OBJDIR)/%.o: %.f90 | $(OBJDIR) $(MODDIR)
	$(FC) $(OPTF) $(FFLAGS) -I. -I$(includedir) -I$(topdir)/src $(INCS) -J$(MODDIR) -c $< -o $@

exec: $(OBJ) | $(BINDIR)
	$(FL) -o $(BINDIR)/$@ $(OPTL) $^ $(LIBSDMUMPS) $(LORDERINGS) $(LIBS) $(RPATH_OPT) $(LIBBLAS) $(LIBOTHERS)

all: directories exec

directories: $(BINDIR) $(OBJDIR) $(MODDIR) $(DATADIR)

$(BINDIR) $(OBJDIR) $(MODDIR) $(DATADIR):
	mkdir -p $@

clean:
	rm -f $(OBJDIR)/*.o $(MODDIR)/*.mod $(BINDIR)/exec $(DATADIR)/*.txt
	rm -f frames/*

.PHONY: all clean directories



