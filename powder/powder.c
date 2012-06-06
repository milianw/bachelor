/*
 * powder.c - Generate powder spectrum from a set of discrete EPR spectrum files
 *
 * Copyright (C) 2012 John Paul Adrian Glaubitz <glaubitz@physik.fu-berlin.de>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Library General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program; if not, write to the
 * Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fftw3.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>

static int debug = 0;

/* typedef struct s_lebedev { */

/*   double phi; */
/*   double theta; */
/*   double weight; */
/* } lebedev; */

typedef struct s_orientation {

  double x;
  double y;
  double z;
  double weight;
} orientation;

typedef struct s_epr_spectrum {

  fftw_complex * B; /* pointer to an array of B field values */
  fftw_complex * I; /* pointer to an array of I field values */
  orientation * O; /* pointer to an array of orientation values */
  int size; /* amount of data points in spectrum */

} epr_spectrum;

/** 
 * alloc_epr_spectrum: allocate memory for an EPR spectrum
 * 
 * @param spectrum - pointer to the EPR spectrum struct which is to be allocated
 * @param size - size of the EPR spectrum (number of data points)
 * 
 * @return - returns size of EPR spectrum if successful, zero otherwise
 */
int alloc_epr_spectrum (epr_spectrum * spectrum, int size) {

  if (!spectrum)
    return 0;

  if (!(spectrum->B = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size)))
    return 0;

  if (!(spectrum->I = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size))) {
    fftw_free (spectrum->B);
    return 0;
  }

  if (!(spectrum->O = (orientation *) malloc (sizeof(orientation) * size))) {
    fftw_free (spectrum->B);
    fftw_free (spectrum->I);
    return 0;
  }

  spectrum->size = size;

  return size;
}

/** 
 * dealloc_epr_spectrum: deallocate memory for an EPR spectrum
 * 
 * @param spectrum - pointer to the EPR spectrum struct to be deallocated
 */
void dealloc_epr_spectrum (epr_spectrum * spectrum) {

  if (!spectrum)
    return;

  if (spectrum->B)
    free(spectrum->B);

  if (spectrum->I)
    free(spectrum->I);

  if (spectrum->O)
    free(spectrum->O);

  spectrum->size = 0;
}

/** 
 * determine_line_count: read number of lines containing at least one double value in a text file
 * 
 * @param datafile - file handle of the data file to be line-counted
 * 
 * @return - returns number of lines with at least one double number, zero otherwise
 */
unsigned int determine_line_count (FILE * datafile) {

  unsigned int lines;
  char buffer [256];
  double tmp;

  lines = 0;

  if (datafile)
    while(fgets(buffer, sizeof(buffer), datafile)) {
      if(buffer[0] && buffer[strlen(buffer)-1] == '\n')
	if (sscanf (buffer, "%lg", &tmp) == 1) lines++;
    }
  else
    return 0;

  return lines;
}

/** 
 * read_input_epr_spectrum: read EPR spectrum from data file into an array of struct epr_spectrum
 * 
 * @param path - absolute path to the data file containing the EPR spectrum data
 * @param spectrum - uninitialized pointer to an array which will be allocated and filled with EPR spectrum data
 * 
 * @return - returns the number of EPR spectrum data points successfully parsed, zero otherwise
 */
int read_input_epr_spectrum (FILE * spectrum_file, epr_spectrum * spectrum) {

  int lines;
  int i;
  char buffer [256];

  lines = 0;
  i = 0;

  /* determine size of EPR spectrum array to be allocated */
  if ((lines = determine_line_count (spectrum_file)) > 0) {
    if (!alloc_epr_spectrum (spectrum, lines))
      return 0;
  }
  else
    return 0;

  if(debug)
    printf ("The Lebedev file has %d lines\n", lines);

  /* the file has to be rewound to the beginning, since it is currently standing
     at the end from the previous call to determine_line_count
  */
  rewind (spectrum_file);

  /* parse EPR spectrum input file */
  while(fgets(buffer, sizeof(buffer), spectrum_file)) {
    if (sscanf(buffer, "%lg %lg %lg %lg %lg %lg[^\n]\n", &spectrum->B[i][0], &spectrum->I[i][0], &spectrum->O[i].x, &spectrum->O[i].y, &spectrum->O[i].z, &spectrum->O[i].weight) == 6) {
      if (debug)
	printf ("i: %d B: %lg I: %lg x: %lg y: %lg z: %lg weight: %lg\n", i, spectrum->B[i][0], spectrum->I[i][0], spectrum->O[i].x, spectrum->O[i].y, spectrum->O[i].z, spectrum->O[i].weight);
      i++;
    }
  }

  return lines;
}

/** 
 * fill_spectrum_equidistant: create spectrum with equidistant data points from
 *                            input spectrum, so the data can be used for FFT later
 * 
 * @param spectrum - pointer to an array with spectrum data which will be replaced with a new array,
 *                   the memory for the input array will be deallocated
 * @param accuracy - accuracy of the new grid to be generated (i.e. smallest step size)
 * 
 * @return - returns the size of the new, intercalated array, zero otherwise
 */
int fill_spectrum_equidistant (epr_spectrum * spectrum, double accuracy) {

  int i, j;
  int new_size;
  double B_min;
  double B_max;

  epr_spectrum new_spectrum;

  if (spectrum && spectrum->size >= 2 && accuracy > 0) {

    i = 0;
    B_min = spectrum->B[0][0];
    B_max = spectrum->B[0][0];

    /* find minimum and maximum field values */

    for (i = 0; i < spectrum->size; i++) {
      if (spectrum->B[i][0] < B_min)
        B_min = spectrum->B[i][0];
      if (spectrum->B[i][0] > B_max)
        B_max = spectrum->B[i][0];
    }

    /* allocate new, larger array for spectrum data with accuracy provided */

    new_size = (int) ((B_max - B_min) / accuracy) + 1;

    if (!alloc_epr_spectrum(&new_spectrum, new_size))
      return 0;
    else {
      /* fill B field values of spectrum */
      for (i = 0; i < new_size; i++) {
	new_spectrum.B[i][0] = B_min + i * accuracy;
	new_spectrum.B[i][1] = 0;
	new_spectrum.I[i][0] = 0;
	new_spectrum.I[i][1] = 0;
	new_spectrum.O[i] = (orientation) {0, 0, 0, 0};
      }

      for (i = 0; i < spectrum->size; i++) {
	j = (int) floor ((spectrum->B[i][0] - B_min) / accuracy);
	new_spectrum.B[j][0] = spectrum->B[i][0];
	new_spectrum.B[j][1] = spectrum->B[i][1];
	new_spectrum.I[j][0] = spectrum->I[i][0];
	new_spectrum.I[j][1] = spectrum->I[i][1];
	new_spectrum.O[j] = spectrum->O[i];
      }
    }

    dealloc_epr_spectrum(spectrum);
    (* spectrum) = new_spectrum;
    return new_size;
  }
  else
    return 0;
}

/** 
 * fft_epr_spectrum: perform fast fourier transformation of EPR spectrum (both forward and reverse)
 * 
 * @param spectrum - pointer to an array with the EPR spectrum input data
 * @param size - size of the EPR spectrum input array
 * @param sign - sign (determines direction of FFT, either FFTW_FORWARD (-1) or FFTW_BACKWARD (+1) from fftw3)
 * 
 * @return - return the size of the transformed EPR spectrum array
 */
int fft_epr_spectrum(epr_spectrum * spectrum, int sign) {

  fftw_complex * out;
  fftw_plan plan;

  if(!spectrum)
    if (spectrum->size == 0)
      return 0;

  /* transformed spectrum values will be in the out array */
  out = (fftw_complex *) fftw_malloc (sizeof(fftw_complex) * spectrum->size);

  /* the next three calls perform the actual discrete FFT */
  plan = fftw_plan_dft_1d (spectrum->size, spectrum->I, out, sign, FFTW_ESTIMATE);
  fftw_execute (plan);
  fftw_destroy_plan (plan);

  /* free the memory for the input intensity values */
  fftw_free (spectrum->I);

  /* set the intensity values to point at the transformed intensity data */
  spectrum->I = out;

  return -1;
}

/** 
 * broaden_spectrum: introduce line shape broadening to a spectrum using a given decay
 * 
 * @param spectrum - pointer to the struct of the EPR spectrum to be broadenend
 * @param decay - decay constant to be used for the exponential decay factor
 * 
 * @return - returns size of spectrum if successful, zero otherwise
 */
int broaden_spectrum (epr_spectrum * spectrum, double decay) {

  int i;

  if (spectrum->size > 0) {
    /* since the input spectrum is in frequency domain, do a reverse transform into time domain */
    fft_epr_spectrum (spectrum, FFTW_BACKWARD);

    /* multiply intensity values with Exp decay function */
    for (i = 0; i < spectrum->size; i++) {
      spectrum->I[i][0] *= exp(-(spectrum->B[i][0] * decay));
      spectrum->I[i][1] *= exp(-(spectrum->B[i][0] * decay));
    }

    /* transform the spectrum back into frequency domain */
    fft_epr_spectrum (spectrum, FFTW_FORWARD);

    return spectrum->size;
  }
  else
    return 0;
}

/** 
 * usage: print usage information
 * 
 * @param cmdname - pointer to a string containing the name of the binary (argv[0])
 */
void usage(char * cmdname)
{
  /* printf("Usage: %s --fill --input-spectrum [spectrum file] --output-file [output file] --lebedev-average [lebedev file] --broadening [decay constant] --debug.\n\n\ */
  printf("Usage: %s --fill [accuracy] --input-spectrum [spectrum file] --output-file [output file] --broadening [decay constant] --debug.\n\n\
          Where:\n\
          \n\
          fill - insert additional data points with [accuracy] to EPR to enable discrete FFT\n\
          input-spectrum - input file with EPR spectrum data\n\
          output-file - write output to [output file]\n\
          broadening - perform line broadening through FFT using the supplied [decay constant]\n\
          debug - output extended debug data\n\n", cmdname);

          /* lebedev-average - perform Lebedev average with Lebedev data from [lebedev file]\n\ */
}

int main (int argc, char** argv) {

  int i, c, option_index;
  /* option flags for averaging, fill and broadening; default = 0 => OFF */
  /* int lebedev_average = 0; */
  int fill = 0, broadening = 0;
  double decay;
  double accuracy;
  /* lebedev * lebedev_grid; */
  epr_spectrum spectrum;
  FILE * input_spectrum_file;
  /* FILE * lebedev_file; */
  FILE * output_file;
  char input_spectrum_path [256];
  /* char lebedev_file_path [256]; */
  char output_path [256];

  memset (input_spectrum_path, 0, 256);
  /* memset (lebedev_file_path, 0, 256); */
  memset (output_path, 0, 256);

  static struct option options [] = {
    {"input-spectrum", required_argument, NULL, 'i'},
    {"output-file", required_argument, NULL, 'o'},
    {"fill", required_argument, NULL, 'f'},
    /* {"lebedev-average", required_argument, NULL, 'l'}, */
    {"broadening", required_argument, NULL, 'b'},
    {"debug", no_argument, NULL, 'd'},
    {0, 0, 0, 0}
  };

  /* parse command line options using getopt_long */
  while ((c = getopt_long(argc, argv, "i:o:f:b:d", options, &option_index)) != -1) {

    switch (c) {

    case 'i':
      strncpy (input_spectrum_path, optarg, 255);
      break;

    case 'o':
      strncpy (output_path, optarg, 255);
      break;

    case 'f':
      fill = 1;
      accuracy = atof (optarg);
      break;

    /* case 'l': */
    /*   strncpy (lebedev_file_path, optarg, 255); */
    /*   lebedev_average = 1; */
    /*   break; */

    case 'b':
      broadening = 1;
      decay = atof (optarg);
      break;

    case 'd':
      debug = 1;
      break;

    default:
      usage (argv[0]);
      return 0;
    }
  }

  /* open input files for EPR spectrum and Lebedev grid data, determine file for output or use stdout */
  if (!(input_spectrum_file = fopen (input_spectrum_path, "r"))) {
    printf ("Error opening file for input or no input file specified.\n");
    usage(argv[0]);
    return 0;
  }

  if (!(output_file = fopen (output_path, "w")))
    output_file = stdout;

  /* read input files for EPR spectrum and Lebedev grid data */
  if (!(read_input_epr_spectrum (input_spectrum_file, &spectrum) > 0)) {
    printf ("Error parsing input spectrum, no data points found.\n");
    return -1;
  }

  /* generate equidistant grid if fill parameter is set */
  if (fill)
    fill_spectrum_equidistant (&spectrum, accuracy);

  if (broadening)
    broaden_spectrum (&spectrum, decay);

  /* output new spectrum to stdout/output_file */
  for (i = 0; i < spectrum.size; i++)
    fprintf (output_file, "%lg %lg %lg %lg %lg %lg\n", spectrum.B[i][0], spectrum.I[i][0], spectrum.O[i].x, spectrum.O[i].y, spectrum.O[i].z, spectrum.O[i].weight);

  /* free memory for EPR spectrum */
  dealloc_epr_spectrum (&spectrum);

  return 0;
}
