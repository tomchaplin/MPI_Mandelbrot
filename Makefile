CC = mpicc
LIBS = lib/libbmp.c
OPTS = -lm -O3

all: mandelbrot mandelbrot-serial

mandelbrot: mandelbrot.c
	$(CC) $(OPTS) mandelbrot.c $(LIBS) -o mandelbrot.exe

mandelbrot-serial: mandelbrot-serial.c
	$(CC) $(OPTS) mandelbrot-serial.c $(LIBS) -o mandelbrot-serial.exe
