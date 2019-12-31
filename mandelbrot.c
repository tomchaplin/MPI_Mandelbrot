#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include "lib/libbmp.h"

MPI_Datatype MPI_Compl;
MPI_Datatype MPI_MandelOpts;
MPI_Datatype MPI_Payload;

int MAX_ITER = 1000;
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

typedef struct Palette
{
	int size;
	Pixel* swatch;
} Palette;

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

void colourPixel(double hue, Pixel* pixel, Palette* palette) {
	int colourNumber = floor( hue * ((double)palette->size - 1.0) );
	pixel->r = palette->swatch[colourNumber].r;
	pixel->g = palette->swatch[colourNumber].g;
	pixel->b = palette->swatch[colourNumber].b;
}

void colourFrame(int* iterationGrid, bmp_img* img, MandelOpts* opts, Palette* palette) {
	Pixel pixel;
	int histogram[opts->max_iterations + 1];
	for(int y = 0; y < opts->p_height; y++) {
		for(int x = 0; x < opts->p_width; x++) {
			histogram[ iterationGrid[y*opts->p_width + x] ]++;
		}
	}
	int total = 0;
	for(int i = 0; i <= opts->max_iterations; i++) {
		total += histogram[i];
	}
	double hueList[opts->max_iterations + 1];	
	double currentHue = 0;
	for(int i = 0; i <= opts->max_iterations; i++) {
		hueList[i] = currentHue;
		currentHue += (double)histogram[i] / total;
	}
	for(int y = 0; y < opts->p_height; y++) {
		for(int x = 0; x < opts->p_width; x++) {
			//colourPixel(iterationGrid[y*opts->p_width + x], opts->max_iterations, &pixel);
			currentHue = hueList[ iterationGrid[y*opts->p_width + x] ];
			colourPixel(currentHue, &pixel, palette);
			bmp_pixel_init(
					&img->img_pixels[y][x], 
					pixel.r, pixel.g, pixel.b);
		}
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
void computeRow(MandelOpts* opts, Payload* load) {
	Compl c = opts->startZ;
	c.im = c.im - (load->row - 1) * opts->c_height / (opts->p_height -1);
	double re_jump = opts->c_width / (opts->p_width - 1);
	for(int x = 0; x < opts->p_width; x++) {
		load->iterations_arr[x] = computePoint(c, opts);
		c.re += re_jump;
	}
}

void writePayload(int* iterationGrid, Payload* payload, MandelOpts* opts) {
	int start_point = payload->row * opts->p_width;
	for(int x = 0; x < width; x++) {
		iterationGrid[x + start_point] = payload->iterations_arr[x];
	}
}

void computeMandelbrot(int* iterationGrid, MandelOpts* opts, int message) {
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
			writePayload(iterationGrid, &payload, opts);
		}
	}
	/* Now we have to collect last set of results and tell threads to finish */
	for(int thread = 0; thread < size - 1; thread++) {
		MPI_Recv(&payload, 1, MPI_Payload, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
		MPI_Send(&message, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
		writePayload(iterationGrid, &payload, opts);
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

void setupDirectory(char* directory) {
	time_t rawtime;
	struct tm * timeinfo;
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	sprintf(directory,  "output/%d_%02d_%02d_%02d%02d%02d/",
			timeinfo->tm_year + 1900,
			timeinfo->tm_mon + 1,
			timeinfo->tm_mday,
			timeinfo->tm_hour,
			timeinfo->tm_min,
			timeinfo->tm_sec
			);
	struct stat st = {0};
	if (stat(directory, &st) == -1) {
			mkdir(directory, 0700);
	}
}

double mod1(double x) {
	while (x >= 1.0) {
		x--;
	}
	while (x < 0.0) {
		x++;
	}
	return x;
}

double hsl_tests(double x, double tmp_1, double tmp_2) {
	if(6.0*x < 1.0) {
		return tmp_2 + (tmp_1 - tmp_2) * 6.0 * x;
	} else if(2.0*x < 1.0) {
		return tmp_1;
	} else if(3.0*x < 2.0) {
		return tmp_2 + (tmp_1 - tmp_2) * (2.0/3.0 - x) * 6.0;
	} else {
		return tmp_2;
	}
}

// Accepts HSL
// Hue in (0,360), saturation in (0,1), lightness in (0,1)
void HSL_2_RGB(double *hsl, int *rgb) {
	double hue = hsl[0];
	double saturation = hsl[1];
	double luminance = hsl[2];

	if(saturation == 0) {
		// The colour is grey
		rgb[0] = (int) luminance * 255.0;
		rgb[1] = (int) luminance * 255.0;
		rgb[2] = (int) luminance * 255.0;
	} else {
		// We need to calculate some temporary variables
	double tmp_1, tmp_2;
	if(luminance < 0.5) {
		tmp_1 = luminance * (1.0 + saturation);
	} else {
		tmp_1 = luminance + saturation - luminance*saturation;
	}
		tmp_2 = 2.0*luminance - tmp_1;
		hue = hue/360.0;
		rgb[0] = (int) ( hsl_tests ( mod1(hue + (1.0/3.0)), tmp_1, tmp_2 ) * 255 );
		rgb[1] = (int) ( hsl_tests ( mod1(hue), tmp_1, tmp_2 ) * 255 );
		rgb[2] = (int) ( hsl_tests ( mod1(hue - (1.0/3.0)), tmp_1, tmp_2 ) * 255 );
	}
}


void getPaletteOne(Palette* palette) {
	int rgb[3];
	double hsl[3];
	hsl[0] = 287.0;
	hsl[1] = 1.0;
	for(int i = 0; i <= 17; i++) {
		hsl[2] = (double) i / 100.0;
		HSL_2_RGB(hsl, rgb);
		palette->swatch[i].r = rgb[0];
		palette->swatch[i].g = rgb[1];
		palette->swatch[i].b = rgb[2];
	}

	for(int i=18; i<=70; i++) {
		hsl[0]++;
		hsl[2] = hsl[2] + (0.19 / 53.0);
		HSL_2_RGB(hsl, rgb);
		palette->swatch[i].r = rgb[0];
		palette->swatch[i].g = rgb[1];
		palette->swatch[i].b = rgb[2];
	}

	for(int i=71; i<=90; i++) {
		hsl[0]++;
		hsl[2] = hsl[2] + (0.19 / 64.0);
		HSL_2_RGB(hsl, rgb);
		palette->swatch[i].r = rgb[0];
		palette->swatch[i].g = rgb[1];
		palette->swatch[i].b = rgb[2];
	}

	hsl[0] = 0.0;

	for(int i=91; i<=134; i++) {
		hsl[0]++;
		hsl[2] = hsl[2] + (0.19 / 64.0);
		HSL_2_RGB(hsl, rgb);
		palette->swatch[i].r = rgb[0];
		palette->swatch[i].g = rgb[1];
		palette->swatch[i].b = rgb[2];
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
		MandelOpts opts;
		opts.p_width = width;
		opts.p_height = height;
		opts.max_iterations = MAX_ITER;
		opts.c_width = 1.0;
		opts.escape_radius = 4.0;
		Compl centerZ;
		centerZ.re = -0.8332313247446712;
		centerZ.im = 0.20542799512012438;
		double zoom_factor = 1.05;
		int frames = 600;

		/* Setup palette */
		Palette palette;
		Pixel palette_swatch[135];
		palette.swatch = palette_swatch;
		palette.size = 135;
		getPaletteOne(&palette);

		/* Setup directory */
		char directory[128];
		setupDirectory(directory);

		/* Open up the image file */
		bmp_img img;
		char filename[128];
		/* Do work */
		for(int currentFrame = 0; currentFrame < frames; currentFrame++) {
			bmp_img_init_df(&img, width, height);
			// Figure out the options for this frame
			opts.c_height = opts.c_width * ((double)opts.p_height / (double)opts.p_width);
			opts.startZ = centerZ;
			opts.startZ.re -= 0.5*opts.c_width;
			opts.startZ.im += 0.5*opts.c_height;
			// Work out the mandelbrot
			int message = AWAIT_OPTIONS;
			if(currentFrame == frames -1) {
				message = POISON_PILL;
			};
			int iterationGrid[width * height];
			computeMandelbrot(iterationGrid, &opts, message);
			colourFrame(iterationGrid, &img, &opts, &palette);
			// Figure out the filename
			sprintf(filename, "%sfile%04d.bmp", directory, currentFrame);
			/* Write to file */
			bmp_img_write(&img, filename);
			bmp_img_free(&img);
			printf("Finished frame %04d\n", currentFrame);
			fflush(stdout);
			// Zooooom
			opts.c_width = opts.c_width / zoom_factor;
		}
	} else {
		worker();
	}
	/* Close MPI */
	MPI_Finalize();
	return 0;
}
