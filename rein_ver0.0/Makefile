CC = gcc
CFLAGS = -O3 -Wall -lm 
gsl_lib = -I$(GSL_ROOT)/include -L$(GSL_ROOT)/lib -lgsl -lgslcblas
gsl_include = -I$(GSL_ROOT)/include
addinclude = -I./supp2 -I./supp2/funcs
sacpref =/uufs/chpc.utah.edu/sys/installdir/SAC/101.6a 
SACLIB=`sac-config --prefix=$(sacpref) --libs sac`
SACIO=`sac-config --prefix=$(sacpref)  --libs sacio`
SACFLAGS = `sac-config --cflags`

default: all

all: 
	$(CC) -c $(CFLAGS) supp2/funcs/sacioH.c
	$(CC) -c $(CFLAGS) $(gsl_lib) $(gsl_include) supp2/funcs/complex_pack-2.c
	$(CC) -c $(CFLAGS) $(gsl_lib) $(gsl_include) supp2/funcs/matrix_pack.c
	$(CC) -c $(CFLAGS) $(gsl_lib) $(gsl_include) supp2/funcs/fft_pack.c
	$(CC) -c $(CFLAGS) $(gsl_lib) $(gsl_include) supp2/funcs/calculation_w.c
	$(CC) -c $(CFLAGS) $(gsl_lib) $(gsl_include) supp2/funcs/chris_only.c
	$(CC) -c $(CFLAGS) $(gsl_lib) $(gsl_include) supp2/funcs/acf_functions.c
	$(CC) -c $(CFLAGS) $(gsl_lib) $(gsl_include) $(addinclude) supp2/func_models.c
	$(CC) -c -ld $(CFLAGS) $(gsl_lib) $(gsl_include) $(addinclude) supp2/func_synthetic_iso-3.c
	$(CC) -c $(CFLAGS) $(gsl_lib) $(gsl_include) $(addinclude) supp2/func_readwrite-2.c
	$(CC) -c $(CFLAGS) $(gsl_lib) $(gsl_include) $(addinclude) supp2/func_likelihood-2adv.c
	mpicc -o trinversion_noise_ver10.2-chol3 $(CFLAGS) $(SACFLAGS) main_inversion_mpi_ver9-chol.c func_synthetic_iso-3.o $(SACLIB) $(SACIO) $(gsl_lib) $(gsl_include) $(addinclude) sacioH.o func_likelihood-2adv.o func_models.o func_readwrite-2.o complex_pack-2.o matrix_pack.o fft_pack.o calculation_w.o chris_only.o acf_functions.o

clean:
	rm -rf *.o 
