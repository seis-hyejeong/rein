CC = gcc
CFLAGS = -O2 -Wall -lm
gsl_lib = -L/home/hyejeong/gsl/lib -lgsl -lgslcblas
gsl_include = -I/home/hyejeong/gsl/include

default: all

all: 
	$(CC) -c $(CFLAGS) sacioH.c
	$(CC) -c $(CFLAGS) $(gsl_lib) $(gsl_include) readwrite_w.c
	$(CC) -c $(CFLAGS) $(gsl_lib) $(gsl_include) complex_pack.c
	$(CC) -c $(CFLAGS) $(gsl_lib) $(gsl_include) matrix_pack.c
	$(CC) -c $(CFLAGS) $(gsl_lib) $(gsl_include) fft_pack.c
	$(CC) -c $(CFLAGS) $(gsl_lib) $(gsl_include) calculation_w.c
	$(CC) -c $(CFLAGS) $(gsl_lib) $(gsl_include) chris_only.c
	$(CC) -c $(CFLAGS) $(gsl_lib) $(gsl_include) acf_functions.c
	$(CC) -o paps_water $(CFLAGS) $(gsl_lib) $(gsl_include) paps_water_main.c readwrite_w.o complex_pack.o matrix_pack.o fft_pack.o calculation_w.o chris_only.o sacioH.o
	$(CC) -o acf_water $(CFLAGS) $(gsl_lib) $(gsl_include) acf_tele_main.c readwrite_w.o complex_pack.o matrix_pack.o fft_pack.o calculation_w.o chris_only.o sacioH.o acf_functions.o
	$(CC) -o acf_water_smooth $(CFLAGS) $(gsl_lib) $(gsl_include) acf_tele_smoothing_main.c readwrite_w.o complex_pack.o matrix_pack.o fft_pack.o calculation_w.o chris_only.o sacioH.o acf_functions.o
	$(CC) -o acf_hann_smooth $(CFLAGS) $(gsl_lib) $(gsl_include) acf_tele_hanning_main.c readwrite_w.o complex_pack.o matrix_pack.o fft_pack.o calculation_w.o chris_only.o sacioH.o acf_functions.o
	$(CC) -o hv_water $(CFLAGS) $(gsl_lib) $(gsl_include) hv_tele_main.c readwrite_w.o complex_pack.o matrix_pack.o fft_pack.o calculation_w.o chris_only.o sacioH.o acf_functions.o
	cp paps_water acf_water hv_water acf_hann_smooth acf_water_smooth ../bin/
clean:
	rm -rf *.o
