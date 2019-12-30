#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "lib/libbmp.h"

int MAX_ITER = 250;
int width = 1920;
int height = 1080;
char filename[9] = "file1.bmp";

typedef struct Compl 
{
	float re, im;
} Compl;

typedef struct Pixel
{
	int r, g, b;
} Pixel;

typedef struct MandelOpts
{
	int p_width, p_height;
	int max_iterations;
	float c_width, c_height;
	float escape_radius;
	Compl startZ;
} MandelOpts;

typedef struct Payload
{
	int row;
	Pixel pixel_arr[1920];
} Payload;

int computePoint(Compl c, MandelOpts* opts) {
	int iterations = 0;
	Compl z;
	float z_length, temp;
	z.re = 0.0; z.im = 0.0;
	do {
		temp = z.re * z.re - z.im * z.im + c.re;
		z.im = 2.0 * z.re * z.im + c.im;
		z.re = temp;
		z_length = z.re * z.re + z.im * z.im;
		iterations++;
	} while ( (z_length < opts->escape_radius) && (iterations < opts->max_iterations) );
	return iterations;
}

void colourPixel(int iterations, int max_iterations, Pixel* pixel) {
	int colour = (int)floor( iterations * 255.0 / max_iterations );
	pixel->r = colour;
	pixel->g = colour;
	pixel->b = colour;
}

void computeMandelbrot(bmp_img* img, MandelOpts* opts) {
	Compl c = opts->startZ;
	float re_jump = opts->c_width / (opts->p_width - 1);
	float im_jump = opts->c_height / (opts->p_height - 1);
	int pixelIter;
	Pixel pixel;
	for(int y = 0, x; y < opts->p_height; y++) {
		for(x = 0; x < opts->p_width; x++) {
			pixelIter = computePoint(c, opts);
			colourPixel(pixelIter, opts->max_iterations, &pixel);
			bmp_pixel_init(&img->img_pixels[y][x], pixel.r, pixel.g, pixel.b);
			c.re += re_jump;
		}
		c.re = opts->startZ.re;
		c.im -= im_jump;
	}
}

int main(int argc, char** argv) {
	Compl z;
	z.re = -2.0;
	z.im = 1.0;

	MandelOpts opts;
	opts.p_width = width;
	opts.p_height = height;
	opts.c_width = 3.0;
	opts.c_height = 2.0;
	opts.max_iterations = MAX_ITER;
	opts.escape_radius = 4.0;
	opts.startZ = z;

	bmp_img img;
	bmp_img_init_df(&img, width, height);

	computeMandelbrot(&img, &opts);

	bmp_img_write(&img, "file1.bmp");
	bmp_img_free(&img);
	return 0;
}
