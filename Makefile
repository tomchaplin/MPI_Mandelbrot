CC = mpicc
LIBS = lib/libbmp.c
OPTS = -lm

mandelbrot: mandelbrot.c
	$(CC) $(OPTS) mandelbrot.c $(LIBS) -o mandelbrot

mandelbrot-serial: mandelbrot-serial.c
	$(CC) $(OPTS) mandelbrot-serial.c $(LIBS) -o mandelbrot-serial
