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

void computeRow(MandelOpts* opts, Payload* load) {
	Compl c = opts->startZ;
	c.im = c.im - (load->row - 1) * opts->c_height / (opts->p_height -1);
	float re_jump = opts->c_width / (opts->p_width - 1);
	int pixelIter;
	Pixel pixel;
	for(int x = 0; x < opts->p_width; x++) {
		pixelIter = computePoint(c, opts);
		colourPixel(pixelIter, opts->max_iterations, &pixel);
		load->pixel_arr[x] = pixel;
		c.re += re_jump;
	}
}

void writePayload(bmp_img* img, Payload* payload) {
	for(int x = 0; x < width; x++) {
		bmp_pixel_init(
				&img->img_pixels[payload->row][x], 
				payload->pixel_arr[x].r,
				payload->pixel_arr[x].g,
				payload->pixel_arr[x].b);
	}
}

int main(int argc, char** argv) {
	/* Setup options for our mandelbot */
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
	/* Init MPI */
	int size, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	/* Setup MPI types */

	/* First Compl */
	MPI_Datatype MPI_Compl;
	{
	int nblocks = 2;
	int blocklengths[2] = {1, 1};
	MPI_Datatype types[2] = {MPI_DOUBLE, MPI_DOUBLE};
	MPI_Aint offsets[2];
	offsets[0] = offsetof(Compl, re);
	offsets[1] = offsetof(Compl, im);
	MPI_Type_create_struct(nblocks, blocklengths, offsets, types, &MPI_Compl);
	MPI_Type_commit(&MPI_Compl);
	}

	/* Now for MandelOpts */
	MPI_Datatype MPI_MandelOpts;
	{
	int nblocks = 3;
	int blocklengths[3] = {3,3,1};
	MPI_Datatype types[3] = {MPI_INT, MPI_DOUBLE, MPI_Compl};
	MPI_Aint offsets[3];
	offsets[0] = offsetof(MandelOpts, p_width);
	offsets[1] = offsetof(MandelOpts, c_width);
	offsets[2] = offsetof(MandelOpts, startZ);
	MPI_Type_create_struct(nblocks, blocklengths, offsets, types, &MPI_MandelOpts);
	MPI_Type_commit(&MPI_MandelOpts);
	}

	/* Now for Pixel */
	MPI_Datatype MPI_Pixel;
	{
	int nblocks = 1;
	int blocklengths[1] = {3};
	MPI_Datatype types[1] = {MPI_INT};
	MPI_Aint offsets[1];
	offsets[0] = offsetof(Pixel, r);
	MPI_Type_create_struct(nblocks, blocklengths, offsets, types, &MPI_Pixel);
	MPI_Type_commit(&MPI_Pixel);
	}

	/* Now for Payload */
	MPI_Datatype MPI_Payload;
	{
	int nblocks = 2;
	int blocklengths[2] = {1,1920};
	MPI_Datatype types[3] = {MPI_INT, MPI_Pixel};
	MPI_Aint offsets[2];
	offsets[0] = offsetof(Payload, row);
	offsets[1] = offsetof(Payload, pixel_arr);
	MPI_Type_create_struct(nblocks, blocklengths, offsets, types, &MPI_Payload);
	MPI_Type_commit(&MPI_Payload);
	}

	/* Get MPI rank */
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	/* Broadcast the context to our workers */
	MPI_Bcast(&opts, 1, MPI_MandelOpts, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		/* Open up the image file */
		bmp_img img;
		bmp_img_init_df(&img, width, height);
		int targetThread;
		MPI_Status status;
		Payload payload;
		/* Do the actual work */
		for(int currentRow = 0; currentRow < opts.p_height; currentRow++) {
			if(currentRow  < size - 1) {
				/* We send out the initial work */
				MPI_Send(&currentRow, 1, MPI_INT, currentRow + 1, 0, MPI_COMM_WORLD);
			} else {
				/* We listen for results and distribute work */
				MPI_Recv(&payload, 1, MPI_Payload, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
				MPI_Send(&currentRow, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
				writePayload(&img, &payload);
			}
		}
		/* Now we have to collect last set of results and tell threads to finish */
		int currentRow = -1;
		for(int thread = 0; thread < size - 1; thread++) {
			MPI_Recv(&payload, 1, MPI_Payload, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
			MPI_Send(&currentRow, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
			writePayload(&img, &payload);
		}
		/* Write to file */
		bmp_img_write(&img, "file1.bmp");
		bmp_img_free(&img);
	} else {
		MPI_Status status;
		Payload payload;
		int targetRow;
		MPI_Recv(&targetRow, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		do {
			payload.row = targetRow;
			computeRow(&opts, &payload);
			MPI_Send(&payload, 1, MPI_Payload, 0, 0, MPI_COMM_WORLD);
			MPI_Recv(&targetRow, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		} while ( targetRow != -1 );
	}
	/* Close MPI */
	MPI_Finalize();
	return 0;
}
