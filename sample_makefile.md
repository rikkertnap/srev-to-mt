# Define the compiler and flags
FC = gfortran
FFLAGS = -O2 -Wall

# List all object files needed to link the final executable
OBJ = globals.o precision_definition.o chains_spheriod.o solver_mod.o main.o

# The final target executable
program.exe: $(OBJ)
	$(FC) $(FFLAGS) -o program.exe $(OBJ)

# -----------------------------------------------------------------
# THE DEPENDENCY ORDER (This is what handles your 'use' visibility)
# -----------------------------------------------------------------

# 1. Base level modules (they don't 'use' anything else)
globals.o: globals.f90
precision_definition.o: precision_definition.f90

# 2. Middle tier: chains_spheriod 'uses' globals and precision
chains_spheriod.o: chains_spheriod.f90 globals.o precision_definition.o
	$(FC) $(FFLAGS) -c chains_spheriod.f90

# 3. Higher tier: solver_mod 'uses' chains_spheriod
solver_mod.o: solver_mod.f90 chains_spheriod.o
	$(FC) $(FFLAGS) -c solver_mod.f90

# 4. Top tier: main 'uses' globals, chains_spheriod, and solver_mod
main.o: main.f90 globals.o chains_spheriod.o solver_mod.o
	$(FC) $(FFLAGS) -c main.f90

# Standard pattern rule for the basic independent files
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

# Clean rule to wipe binaries and module maps
clean:
	rm -f *.o *.mod program.exe
