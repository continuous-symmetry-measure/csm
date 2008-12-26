#fortran
FORTRAN = gfortran

# the compiler flags
FLAGS = -O3 -fno-range-check -fconvert=little-endian

# the first target. All executables will be made when 'make' is
# invoked with no target
all: mir inv c2 

c2: c2_3d_symmetry.f
	$(FORTRAN) $(FLAGS) $^ -o $@

inv: inv3d_symmetry.f
	$(FORTRAN) $(FLAGS) $^ -o $@

mir: mir3d_symmetry.f
	$(FORTRAN) $(FLAGS) $^  -o $@

