/*
	Rhea Pottathuparambil
	1001551132
*/


#include "bitmap.h"

#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <sys/time.h>
#include <pthread.h>
#include <math.h>

//making a struct to pass multiple arguments to the pthread function
//this struct contains details needed in the pthread function
struct data
{
	int start;
	int end;
	double xmin;
	double xmax;
	double ymin;
	double ymax;
	int max;
};

int iteration_to_color(int i, int max);
int iterations_at_point(double x, double y, int max);
void *compute_image(void *t); 
//declared this as global as it was harder to pass it as an argument to fucntions everytim
struct bitmap *bm; 


void show_help()
{
	printf("Use: mandel [options]\n");
	printf("Where options are:\n");
	printf("-m <max>    The maximum number of iterations per point. (default=1000)\n");
	printf("-x <coord>  X coordinate of image center point. (default=0)\n");
	printf("-y <coord>  Y coordinate of image center point. (default=0)\n");
	printf("-s <scale>  Scale of the image in Mandlebrot coordinates. (default=4)\n");
	printf("-W <pixels> Width of the image in pixels. (default=500)\n");
	printf("-H <pixels> Height of the image in pixels. (default=500)\n");
	printf("-o <file>   Set output file. (default=mandel.bmp)\n");
	printf("-h          Show this help text.\n");
	printf("-n <thread> Number of threads\n");
	printf("\nSome examples are:\n");
	printf("mandel -x -0.5 -y -0.5 -s 0.2\n");
	printf("mandel -x -.38 -y -.665 -s .05 -m 100\n");
	printf("mandel -x 0.286932 -y 0.014287 -s .0005 -m 1000 -n 50\n\n");
}

int main(int argc, char *argv[])
{	
	char c;

	// These are the default configuration values used
	// if no command line arguments are given.

	const char *outfile = "mandel.bmp";
	double xcenter = 0;
	double ycenter = 0;
	double scale = 4;
	int image_width = 500;
	int image_height = 500;
	int max = 1000;
	int n = 1; //Initial number of threads is 1

	// For each command line argument given,
	// override the appropriate configuration value.

	while ((c = getopt(argc, argv, "x:y:s:W:H:m:o:n:h")) != -1)
	{
		switch (c)
		{
		case 'x':
			xcenter = atof(optarg);
			break;
		case 'y':
			ycenter = atof(optarg);
			break;
		case 's':
			scale = atof(optarg);
			break;
		case 'W':
			image_width = atoi(optarg);
			break;
		case 'H':
			image_height = atoi(optarg);
			break;
		case 'm':
			max = atoi(optarg);
			break;
		case 'o':
			outfile = optarg;
			break;
		case 'n': 
			n = atoi(optarg); //getting the number of threads from arguments
			break;
		case 'h':
			show_help();
			exit(1);
			break;
		}
	}

	// Display the configuration of the image.
	printf("mandel: x=%lf y=%lf scale=%lf max=%d outfile=%s num_threads=%d Height= %d Width=%d\n", xcenter, ycenter, scale, max, outfile, n,image_height,image_width);

	struct data temp[n]; //making n elements for the struct as I have to pass it as an argument everytime

	// Create a bitmap of the appropriate size.
	bm = bitmap_create(image_width, image_height);

	// Fill it with a dark blue, for debugging
	bitmap_reset(bm, MAKE_RGBA(0, 0, 255, 0));
	


	int i = 0; 
	pthread_t thread[n]; //declaring n threads (the threads are not created yet)
	int threshold = ceil(image_height/n);
	for (i = 0; i < n; i++)
	{
		//storing all values in the element of struct array.
		temp[i].start = i*threshold; //calculating the start height and end height for each thread
		temp[i].end = (i+1)*threshold;
		if(temp[i].end > image_height)
			temp[i].end = image_height;
		temp[i].xmin = xcenter - scale;
		temp[i].xmax = xcenter + scale;
		temp[i].ymin = ycenter - scale;
		temp[i].ymax = ycenter + scale;
		temp[i].max = max;

		if(pthread_create(&thread[i], NULL, compute_image, (void *)&temp[i])) //creating threads 
		{
			perror("Error while creating threads.");
			exit(1);
		}			
	}
	
	for (i = 0; i < n; i++)
	{
		if (pthread_join(thread[i], NULL)) //joining each thread.
		{
			perror("Error while joining threads.");
			exit(EXIT_FAILURE);
		}
	}
	
	// Save the image in the stated file.
	if (!bitmap_save(bm, outfile))
	{
		fprintf(stderr, "mandel: couldn't write to %s: %s\n", outfile, strerror(errno));
		return 1;
	}												

	return 0;
}

/*
Compute an entire Mandelbrot image, writing each point to the given bitmap.
Scale the image to the range (xmin-xmax,ymin-ymax), limiting iterations to "max"
*/

void *compute_image(void *temp)
{

	struct data *t = (struct data *)temp; //casting the void * pointer
	//getting all values and storing them in seperate variables
	double xmin = t->xmin;
	double xmax = t->xmax;
	double ymin = t->ymin;
	double ymax = t->ymax;
	int start = t->start;
	int end = t->end;
	int max = t->max;

	int i, j;

	int width = bitmap_width(bm);
	int height = bitmap_height(bm);

	
	//printf("start = %d, end = %d\n",start,end);
	// For every pixel in the image...
	for (j = start; j < end; j++) 
	{
		//printf("j = %d\n", j);
		for (i = 0; i < width; i++)
		{

			// Determine the point in x,y space for that pixel.
			double x = xmin + i * (xmax - xmin) / width;
			double y = ymin + j * (ymax - ymin) / height;

			// Compute the iterations at that point.
			int iters = iterations_at_point(x, y, max);

			// Set the pixel in the bitmap.
			bitmap_set(bm, i, j, iters);
		}
	}
	
	//printf("start = %d, end = %d, ended.\n",start,end);
	pthread_exit(NULL); //exiting the thread
}

/*
Return the number of iterations at point x, y
in the Mandelbrot space, up to a maximum of max.
*/

int iterations_at_point(double x, double y, int max) 
{
	double x0 = x;
	double y0 = y;

	int iter = 0;

	while ((x * x + y * y <= 4) && iter < max)
	{

		double xt = x * x - y * y + x0;
		double yt = 2 * x * y + y0;

		x = xt;
		y = yt;

		iter++;
	}

	return iteration_to_color(iter, max);
}

/*
Convert a iteration number to an RGBA color.
Here, we just scale to gray with a maximum of imax.
Modify this function to make more interesting colors.
*/

int iteration_to_color(int i, int max)
{
	int gray = 255 * i / max;
	return MAKE_RGBA(gray, gray, gray, 0);
}
