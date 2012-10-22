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
#include <dirent.h>

static int debug = 0;

/* struct to hold orientation (x,y,z)
 * of a spin plus its weighting
 */
typedef struct s_orientation {

  double x;
  double y;
  double z;
  double weight;
} orientation;

/* struct to hold a complete EPR
 * spectrum including B field, intensity,
 * orientation (using s_orientation);
 * since the size of the spectrum is dynamic,
 * the array size is stored in 'size'
 */
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
 * @param have_orientations - set whether allocate memory for orientation data
 * 
 * @return - returns size of EPR spectrum if successful, zero otherwise
 */
int alloc_epr_spectrum (epr_spectrum * spectrum, int size, int have_orientations) {

  if (!spectrum)
    return 0;

  /* fftw provides its own function to allocate data for complex numbers,
   * so we use this one rather than the plain libc malloc; fftw_malloc
   * makes sure that the array boundaries are SIMD-aligned (MMX, SSE etc)
   * which increases performance on SIMD-enabled systems; complex numbers
   * for the B field were chosen for new particular reason, the imaginary
   * part is currently always set to zero
   */
  if (!(spectrum->B = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size)))
    return 0;

  /* allocate array for the intensities the same way we do it for the B field,
   * in the case of the intensity we are actually using both the real and imaginary
   * part of the complex numbers
   */
  if (!(spectrum->I = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size))) {
    fftw_free (spectrum->B);
    return 0;
  }

  /* allocate memory for the orientations only when explicitly requested; this
   * is useful for EPR spectra where we've already summed over all orientations
   * and don't need to keep track of the orientation of a single spin anymore
   */
  if (have_orientations) {
    if (!(spectrum->O = (orientation *) malloc (sizeof(orientation) * size))) {
      fftw_free (spectrum->B);
      fftw_free (spectrum->I);
      return 0;
    }
  }
  else
    spectrum->O = NULL;

  /* return the size of the allocated EPR spectrum,
   * zero otherwise (see above); this way we have
   * a useful return code to test the success of this
   * function call
   */
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
 * determine_line_count: read number of lines containing a valid line in a spectrum file
 * 
 * @param datafile - file handle of the data file to be line-counted
 * 
 * @return - returns number of lines with at least one double number, zero otherwise
 */
unsigned int determine_line_count (FILE * datafile) {

  unsigned int lines;
  char buffer [256];
  double tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;

  lines = 0;

  if (datafile)
    while(fgets(buffer, sizeof(buffer), datafile)) {
      if(buffer[0] && buffer[strlen(buffer)-1] == '\n')
	if (sscanf (buffer, "%lg %lg %lg %lg %lg %lg[^\n]\n", &tmp1, &tmp2, &tmp3, &tmp4, &tmp5, &tmp6) == 6) {
	  lines++;
	  buffer[0] = '\0';
	}
    }
  else
    return 0;

  /* the file has to be rewound to the beginning before it can
     be used for further reading in
  */
  rewind (datafile);

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
  /* double B; */

  lines = 0;
  i = 0;

  /* determine size of EPR spectrum array to be allocated */
  if ((lines = determine_line_count (spectrum_file)) > 0) {
    if (!alloc_epr_spectrum (spectrum, lines, 1))
      return 0;
  }
  else
    return 0;

  if(debug)
    printf ("The read spectrum file has %d lines\n", lines);

  /* parse EPR spectrum input file */
  while(fgets(buffer, sizeof(buffer), spectrum_file)) {
    if (sscanf(buffer, "%lg %lg %lg %lg %lg %lg[^\n]\n", &spectrum->B[i][0], &spectrum->I[i][0], &spectrum->O[i].x, &spectrum->O[i].y, &spectrum->O[i].z, &spectrum->O[i].weight) == 6) {
      if (debug)
	printf ("i: %d B: %lg I: %lg x: %lg y: %lg z: %lg weight: %lg\n", i, spectrum->B[i][0], spectrum->I[i][0], spectrum->O[i].x, spectrum->O[i].y, spectrum->O[i].z, spectrum->O[i].weight);
      i++;
      buffer[0] = '\0';
    }
  }

  return i;
}

/** 
 * create_uniform_grid: create uniform grid with equidistant B field values
 * 
 * @param spectrum - pointer to an array with spectrum data which will be replaced with a new array,
 *                   the memory for the input array will be deallocated
 * @param accuracy - accuracy of the new grid to be generated (i.e. smallest step size)
 * 
 * @return - returns the size of the new, intercalated array, zero otherwise
 */
int create_uniform_grid (epr_spectrum * spectrum, double accuracy) {

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
      /* if (spectrum->B[i][0] < B_min) */
      /*   B_min = spectrum->B[i][0]; */
      if (spectrum->B[i][0] > B_max)
        B_max = spectrum->B[i][0];
    }

    /* allocate new, larger array for spectrum data with accuracy provided */

    new_size = (int) ((B_max - B_min) / accuracy) + 1;

    if (!alloc_epr_spectrum(&new_spectrum, new_size, 1))
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
	if (!isnan(spectrum->I[i][0]))
	  new_spectrum.I[j][0] += spectrum->I[i][0] * spectrum->O[i].weight;
	if (!isnan(spectrum->I[i][1]))
	  new_spectrum.I[j][1] += spectrum->I[i][1] * spectrum->O[i].weight;
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
 * average_multiple_spectra - average over multiple spectra
 * 
 * @param input_spectra - an array containing pointers to EPR spectra
 * @param averaged_spectrum - pointer to the averaged output spectrum
 * @param spectrum_count - number of spectra in the input spectrum array
 * @param accuracy - accuracy of the uniform grid of the input spectra
 * 
 * @return - returns size of the averaged spectrum on success, zero in case of error
 */
int average_multiple_spectra(epr_spectrum * input_spectra, epr_spectrum * averaged_spectrum, int spectrum_count, double accuracy) {

  int averaged_size;
  double B_min, B_max;
  int i, j, k;

  B_min = input_spectra[0].B[0][0];
  B_max = input_spectra[0].B[input_spectra[0].size - 1][0];

  /* find minimum and maximum B field for all input spectra */

  for (i = 1; i < spectrum_count; i++) {
    if (input_spectra[i].B[0][0] < B_min)
      B_min = input_spectra[i].B[0][0];

    if (input_spectra[i].B[input_spectra[i].size - 1][0] > B_max)
      B_max = input_spectra[i].B[input_spectra[i].size - 1][0];
  }

  /* allocate spectrum large enough to accomodate the range of
     all input spectra */

  averaged_size = (int) ((B_max - B_min) / accuracy) + 1;

  if (!alloc_epr_spectrum(averaged_spectrum, averaged_size, 0))
    return 0;
  else {
    for (i = 0; i < averaged_spectrum->size; i++) {
      averaged_spectrum->B[i][0] = B_min + i * accuracy;
      averaged_spectrum->B[i][1] = 0;
      averaged_spectrum->I[i][0] = 0;
      averaged_spectrum->I[i][1] = 0;
    }
    
    for (i = 0; i < spectrum_count; i++)
      for (j = 0; j < input_spectra[i].size; j++) {
	k = (int) floor ((input_spectra[i].B[j][0] - B_min) / accuracy);
	//if (!isnan(input_spectra[j]->I[i][0]))
	averaged_spectrum->I[k][0] += input_spectra[i].I[j][0];
	//if (!isnan(input_spectra[j]->I[i][1]))
	averaged_spectrum->I[k][1] += input_spectra[i].I[j][1];
      }
  }
  
  return averaged_size;
}

/** 
 * usage: print usage information
 * 
 * @param cmdname - pointer to a string containing the name of the binary (argv[0])
 */
void usage(char * cmdname)
{
  /* printf("Usage: %s --fill --input-spectrum [spectrum file] --output-file [output file] --lebedev-average [lebedev file] --broadening [decay constant] --debug.\n\n\ */
  printf("Usage: %s --uniform, -u [accuracy] --input-directory, -i [input spectrum directory] --output-directory, -o [output spectrum directory] --broadening, -b [decay constant] --average, -a --debug, -d.\n\n\
          Where:\n\
          \n\
          uniform - insert additional data points with [accuracy] to EPR to create a uniform grid\n\
          input-directory - directory containing input EPR spectra\n\
          output-directory - directory for the output EPR spectra\n\
          broadening - perform line broadening through FFT using the supplied [decay constant]\n\
          average - perform average over all processed spectra\n\
          debug - output extended debug data\n\n", cmdname);
}

int main (int argc, char** argv) {

  int i, c, option_index, spectrum_count, spectrum_index;
  /* option flags for uniform, broadening, average; default = 0 => OFF */
  int uniform = 0, broadening = 0, average = 0, output = 0;
  double decay;
  double accuracy;
  epr_spectrum * input_spectra;
  epr_spectrum averaged_spectrum;
  DIR * input_directory;
  DIR * output_directory;
  struct dirent * input_directory_entry;
  struct dirent * output_directory_entry;
  FILE * input_spectrum_file;
  FILE * output_spectrum_file;
  FILE * averaged_spectrum_file;
  char input_directory_path [256];
  char input_file_path [512];
  char output_directory_path [256];
  char output_file_path [512];
  char average_file_path [256];

  memset (input_directory_path, 0, 256);
  memset (input_file_path, 0, 512);
  memset (output_directory_path, 0, 256);
  memset (output_file_path, 0, 512);
  memset (average_file_path, 0, 256);

  static struct option options [] = {
    {"input-directory", required_argument, NULL, 'i'},
    {"output-directory", required_argument, NULL, 'o'},
    {"uniform", required_argument, NULL, 'u'},
    {"broadening", required_argument, NULL, 'b'},
    {"average", required_argument, NULL, 'a'},
    {"debug", no_argument, NULL, 'd'},
    {0, 0, 0, 0}
  };

  /* parse command line options using getopt_long */
  while ((c = getopt_long(argc, argv, "i:o:u:b:a:d", options, &option_index)) != -1) {

    switch (c) {

    case 'i':
      strncpy (input_directory_path, optarg, 255);
      break;

    case 'o':
      output = 1;
      strncpy (output_directory_path, optarg, 255);
      break;

    case 'u':
      uniform = 1;
      accuracy = atof (optarg);
      break;

    case 'b':
      broadening = 1;
      decay = atof (optarg);
      break;

    case 'a':
      average = 1;
      strncpy (average_file_path, optarg, 255);
      break;

    case 'd':
      debug = 1;
      break;

    default:
      usage (argv[0]);
      return 0;
    }
  }

  if (average && !uniform) {
    printf ("Error: Averaging over multiple spectra requires all grids to be of uniform size. Abort.\n");
    return 0;
  }

  if (broadening && !uniform) {
    printf ("Error: Broadening spectra requires all grids to be of uniform size. Abort.\n");
    return 0;
  }

  /* try to open input directory containing the EPR spectrum data */
  if (!(input_directory = opendir (input_directory_path))) {
    printf ("Error opening directory for input or no input directory specified.\n");
    usage(argv[0]);
    return 0;
  }

  /* try to open output directory; this directory will also be used
   * as the input directory for the averaging process, if applicable
   */
  if (output && !(output_directory = opendir (output_directory_path))) {
    printf ("Error opening directory for output or no input directory specified.\n");
    usage(argv[0]);
    return 0;
  }

  spectrum_count = 0;
  spectrum_index = 0;

  /* determine number of input spectrum files in directory
   * and allocate appropriate amount of memory for FILE handlers
   */

  while (input_directory_entry = readdir (input_directory)) {
    if (input_directory_entry->d_type == DT_REG)
      spectrum_count++;
  }

  rewinddir (input_directory);

  /* allocate memory for array to hold the input spectra structs */
  input_spectra = malloc (spectrum_count * sizeof (epr_spectrum));

  /* main loop for reading and processing input spectra */
  while (input_directory_entry = readdir (input_directory)) {
    
    if (input_directory_entry->d_type == DT_REG) {
      strncpy (input_file_path, input_directory_path, 256);
      strncat (input_file_path, input_directory_entry->d_name, 256);
      
      strncpy (output_file_path, output_directory_path, 256);
      strncat (output_file_path, input_directory_entry->d_name, 256);
      
      if (!(input_spectrum_file = fopen (input_file_path, "r"))) {
	printf ("Error opening input spectrum file %s, skipping.\n", input_directory_entry->d_name);
	continue;
      }
      
      /* read input EPR spectrum into memory */
      if (!(read_input_epr_spectrum (input_spectrum_file, &input_spectra[spectrum_index]) > 0)) {
	printf ("Error parsing input spectrum file %s, no data points found. Skipping.\n", input_directory_entry->d_name);
	continue;
      }

      if (output && !(output_spectrum_file = fopen (output_file_path, "w"))) {
	printf ("Error opening output spectrum file %s, skipping.\n", input_directory_entry->d_name);
	continue;
      }
      
      /* generate equidistant grid if uniform flag is set */
      if (uniform)
	create_uniform_grid (&input_spectra[spectrum_index], accuracy);
      
      if (broadening)
	broaden_spectrum (&input_spectra[spectrum_index], decay);
      
      /* output new spectrum to stdout/output_file */
      if (output)
	for (i = 0; i < input_spectra[spectrum_index].size; i++)
	  fprintf (output_spectrum_file, "%lg %lg %lg %lg %lg %lg\n", input_spectra[spectrum_index].B[i][0], input_spectra[spectrum_index].I[i][0], input_spectra[spectrum_index].O[i].x, input_spectra[spectrum_index].O[i].y, input_spectra[spectrum_index].O[i].z, input_spectra[spectrum_index].O[i].weight);

      if (debug)
	printf ("The line count of %s is: %i\n", output_file_path, input_spectra[spectrum_index].size);
      
      fclose (input_spectrum_file);
      if (output)
	fclose (output_spectrum_file);
      spectrum_index++;
    }
  }

  if(average) {
    //strncpy (output_file_path, output_directory_path, 256);
    //strncat (output_file_path, "averaged-spectrum.data", 512);
    if (!(averaged_spectrum_file = fopen (average_file_path, "w"))) {
      printf ("Error opening output file for averaged spectrum. Abort.");
      return -1;
    }

    if (!average_multiple_spectra(input_spectra, &averaged_spectrum, spectrum_count, accuracy))
      printf ("Error averaging input spectra.\n");
    else {
      if (debug)
	printf ("The size of the averaged spectrum is: %i\n", averaged_spectrum.size);
      
      for (i = 0; i < averaged_spectrum.size; i++)
	fprintf (averaged_spectrum_file, "%lg %lg\n", averaged_spectrum.B[i][0], averaged_spectrum.I[i][0]);
    }

    dealloc_epr_spectrum (&averaged_spectrum);
  }

  for (i = 0; i < spectrum_count; i++)
    dealloc_epr_spectrum (&input_spectra[i]);

  return 0;
}
