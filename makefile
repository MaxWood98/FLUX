#build arguments  #-Wall -fbounds-check -ffpe-trap=underflow,zero -fopt-info-optimized=$@_opt.dat
# buildargs = -O0 -Wall -fbounds-check -ffpe-trap=zero,overflow,underflow,invalid,denormal -fbacktrace -fcheck=all -fopenmp
# buildargs = -O0 -fopenmp
buildargs = -O2 -fopenmp
# buildargs = -O2 -fopt-info-optimized=$@_opt.dat

#build settings
buildsettings = -J obj

#directories
OBJDIR = obj/

#list object file names -> 1 for each .f90 in the correct compilation order
OBJS = $(addprefix $(OBJDIR), \
		io_utilities_module.o\
		flux_data_methods.o\
		flux_io.o\
		flux_edge_flux.o\
		flux_solve.o\
		)

#object patturn rules -> for every file in "dir"/*.f90, make the file *.o in $(OBJDIR) from it
$(OBJDIR)%.o : src/%.f90
	gfortran $(buildsettings) $(buildargs) -c $< -o $@

$(OBJDIR)%.o : io_utilities/%.f90
	gfortran $(buildsettings) $(buildargs) -c $< -o $@

#main build procedure
build: flux_link

#linking procedure
flux_link: $(OBJS) $(addprefix $(OBJDIR), flux_main.o)
	gfortran -o flux $^ $(buildsettings) -I obj $(buildargs) 

#clean procedure 
clean: 
	rm obj/*.mod
	rm obj/*.o 
	rm obj/*.dat