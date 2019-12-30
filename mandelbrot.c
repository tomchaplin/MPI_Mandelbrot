#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "lib/libbmp.h"

MPI_Datatype MPI_Compl;
MPI_Datatype MPI_MandelOpts;
MPI_Datatype MPI_Payload;

int MAX_ITER = 250;
int width = 1920;
int height = 1080;

const int AWAIT_OPTIONS = -2;
const int POISON_PILL = -1;

typedef struct Compl 
{
	double re, im;
} Compl;

typedef struct Pixel
{
	int r, g, b;
} Pixel;

typedef struct MandelOpts
{
	int p_width, p_height;
	int max_iterations;
	double c_width, c_height;
	double escape_radius;
	Compl startZ;
} MandelOpts;

typedef struct Payload
{
	int row;
	int iterations_arr[1920];
} Payload;

void setupTypes(void) {
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
	{
	int nblocks = 2;
	int blocklengths[2] = {1,width};
	MPI_Datatype types[2] = {MPI_INT, MPI_INT};
	MPI_Aint offsets[2];
	offsets[0] = offsetof(Payload, row);
	offsets[1] = offsetof(Payload, iterations_arr);
	MPI_Type_create_struct(nblocks, blocklengths, offsets, types, &MPI_Payload);
	MPI_Type_commit(&MPI_Payload);
	}
}

int computePoint(Compl c, MandelOpts* opts) {
	int iterations = 0;
	Compl z;
	double z_length, temp;
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
	double re_jump = opts->c_width / (opts->p_width - 1);
	for(int x = 0; x < opts->p_width; x++) {
		load->iterations_arr[x] = computePoint(c, opts);
		c.re += re_jump;
	}
}

void writePayload(bmp_img* img, Payload* payload, MandelOpts* opts) {
	Pixel pixel;
	for(int x = 0; x < width; x++) {
		colourPixel(payload->iterations_arr[x], opts->max_iterations, &pixel);
		bmp_pixel_init(
				&img->img_pixels[payload->row][x], 
				pixel.r,
				pixel.g,
				pixel.b);
	}
}

void computeMandelbrot(bmp_img* img, MandelOpts* opts) {
	int size;
	int targetThread;
	MPI_Status status;
	Payload payload;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	for(int currentRow = 0; currentRow < opts->p_height; currentRow++) {
		if(currentRow  < size - 1) {
			// Send the options
			MPI_Send(opts, 1, MPI_MandelOpts, currentRow + 1, 0, MPI_COMM_WORLD);
			/* We send out the initial work */
			MPI_Send(&currentRow, 1, MPI_INT, currentRow + 1, 0, MPI_COMM_WORLD);
		} else {
			/* We listen for results and distribute work */
			MPI_Recv(&payload, 1, MPI_Payload, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
			MPI_Send(&currentRow, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
			writePayload(img, &payload, opts);
		}
	}
	/* Now we have to collect last set of results and tell threads to finish */
	int message = POISON_PILL;
	for(int thread = 0; thread < size - 1; thread++) {
		MPI_Recv(&payload, 1, MPI_Payload, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
		MPI_Send(&message, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
		writePayload(img, &payload, opts);
	}
}

void worker(void) {
	MandelOpts opts;
	MPI_Status status;
	Payload payload;
	int targetRow;
	// Get the first options
	MPI_Recv(&opts, 1, MPI_MandelOpts, 0, 0, MPI_COMM_WORLD, &status);
	// Get the first row
	MPI_Recv(&targetRow, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
	while ( targetRow != POISON_PILL ) {
		payload.row = targetRow;
		computeRow(&opts, &payload);
		MPI_Send(&payload, 1, MPI_Payload, 0, 0, MPI_COMM_WORLD);
		MPI_Recv(&targetRow, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		if(targetRow == AWAIT_OPTIONS) {
			// Get the new options
			MPI_Recv(&opts, 1, MPI_MandelOpts, 0, 0, MPI_COMM_WORLD, &status);
			// Get the first job on the new options
			MPI_Recv(&targetRow, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		}
	}
}

int main(int argc, char** argv) {
	/* Init MPI */
	MPI_Init(&argc, &argv);
	/* Setup MPI types */
	setupTypes();
	/* Get rank */
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0) {
		/* Setup options for our mandelbot */
		Compl z;
		z.re = -2.5;
		MandelOpts opts;
		opts.p_width = width;
		opts.p_height = height;
		opts.c_width = 4.0;
		opts.c_height = opts.c_width / width * height;
		z.im = opts.c_height / 2.0;
		opts.max_iterations = MAX_ITER;
		opts.escape_radius = 4.0;
		opts.startZ = z;

		/* Open up the image file */
		bmp_img img;
		bmp_img_init_df(&img, width, height);
		/* Do work */
		computeMandelbrot(&img, &opts);
		/* Write to file */
		bmp_img_write(&img, "file1.bmp");
		bmp_img_free(&img);
	} else {
		worker();
	}
	/* Close MPI */
	MPI_Finalize();
	return 0;
}
