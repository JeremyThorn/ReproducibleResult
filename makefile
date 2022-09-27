run_info_directory = synthetic_data
make:
	gfortran -Wall -fcheck=all rmd2.3.f95 -o rmd2
	./rmd2 $(run_info_directory)/rmd2_infile.dat $(run_info_directory)/rmd2_qfile.dat 
