#build arguments  #-Wall -fbounds-check -ffpe-trap=underflow,zero -fopt-info-optimized=$@_opt.dat
buildargs = -O0 -Wall -fbounds-check -ffpe-trap=zero,overflow,underflow,invalid,denormal -fbacktrace -fcheck=all
# buildargs = -O2
# buildargs = -O2 -fopt-info-optimized=$@_opt.dat

#build settings
buildsettings = -J obj

#directories
OBJDIR = obj/

#list object file names -> 1 for each .f90 in the correct compilation order
OBJS = $(addprefix $(OBJDIR), \
		io_utilities_module.o\
		facevertex_module.o\
		halfedge_module.o\
		geometry_kdtree_module.o\
		hex_data_methods_module.o\
		hex_geometry.o\
		hex_utilities.o\
		hex_postprocess.o\
		hex_tritree.o\
		hex_io.o\
		hex_2d.o\
		hex_gradient.o\
		)

#object patturn rules -> for every file in "dir"/*.f90, make the file *.o in $(OBJDIR) from it
$(OBJDIR)%.o : src/%.f90
	gfortran $(buildsettings) $(buildargs) -c $< -o $@

$(OBJDIR)%.o : io_utilities/%.f90
	gfortran $(buildsettings) $(buildargs) -c $< -o $@

$(OBJDIR)%.o : halfedge/%.f90
	gfortran $(buildsettings) $(buildargs) -c $< -o $@

$(OBJDIR)%.o : kdtree/%.f90
	gfortran $(buildsettings) $(buildargs) -c $< -o $@

#main build procedure
build: hex_link

#linking procedure
hex_link: $(OBJS) $(addprefix $(OBJDIR), hex_main.o)
	gfortran -o hex $^ $(buildsettings) -I obj $(buildargs) 

#clean procedure 
clean: 
	rm obj/*.mod
	rm obj/*.o 
	rm obj/*.dat