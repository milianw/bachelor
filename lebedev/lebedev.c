/*
 * lebedev.c - Convert Lebedev data to cartesian and vice versa
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
#include <math.h>
#include <unistd.h>
#include <getopt.h>

static int debug = 0;

typedef struct s_lebedev {

  double phi;
  double theta;
  double weight;
} lebedev;

typedef struct s_cartesian {

  double x;
  double y;
  double z;
  double weight;
} cartesian;

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
 * read_lebedev_data: read Lebedev grid data from file into an array of struct lebedev_grid
 * 
 * @param lebedev_file - file handle to the data file containing the Lebedev grid data
 * @param lebedev_data - uninitialized pointer to an array which will be allocated and filled with grid data
 * 
 * @return - returns the number of Lebedev grid points successfully parsed, zero otherwise
 */
int read_lebedev_data (FILE * lebedev_file, lebedev ** lebedev_data) {

  int lines;
  int i;
  char buffer [256];

  lines = 0;
  i = 0;

  /* determine size of Lebedev array to be allocated */
  if ((lines = determine_line_count (lebedev_file)) > 0) {
    if (!(*lebedev_data = malloc (sizeof(lebedev)*lines)))
      return 0;
  }
  else
    return 0;

  if(debug)
    printf ("The input file has %d lines\n", lines);

  /* the file has to be rewound to the beginning, since it is currently standing
     at the end from the previous call to determine_line_count
  */
  rewind (lebedev_file);

  /* parse Lebedev grid data from file, use exponential form for floating point data if applicable */
  while(fgets(buffer, sizeof(buffer), lebedev_file)) {
    if (sscanf(buffer, "%lf %lf %lf[^\n]\n", &(* lebedev_data)[i].phi, &(* lebedev_data)[i].theta, &(* lebedev_data)[i].weight) == 3) {
      if (debug)
	printf ("i: %d phi: %.15lf theta: %.15lf weight: %.15lf\n", i, (* lebedev_data)[i].phi, (* lebedev_data)[i].theta, (* lebedev_data)[i].weight);
      i++;
    }
  }

  return lines;
}

/** 
 * read_cartesian_data: read Lebedev grid data from file into an array of struct lebedev_grid
 * 
 * @param cartesian_file - file handle to the data file containing the cartesian data
 * @param cartesian_data - uninitialized pointer to an array which will be allocated and filled with cartesian data
 * 
 * @return - returns the number of cartesian data points successfully parsed, zero otherwise
 */
int read_cartesian_data (FILE * cartesian_file, cartesian ** cartesian_data) {

  int lines;
  int i;
  char buffer [256];

  lines = 0;
  i = 0;

  /* determine size of Lebedev array to be allocated */
  if ((lines = determine_line_count (cartesian_file)) > 0) {
    if (!(*cartesian_data = malloc (sizeof(cartesian)*lines)))
      return 0;
  }
  else
    return 0;

  if(debug)
    printf ("The input file has %d lines\n", lines);

  /* the file has to be rewound to the beginning, since it is currently standing
     at the end from the previous call to determine_line_count
  */
  rewind (cartesian_file);

  /* parse Lebedev grid data from file, use exponential form for floating point data if applicable */
  while(fgets(buffer, sizeof(buffer), cartesian_file)) {
    if (sscanf(buffer, "%lf %lf %lf[^\n]\n", &(* cartesian_data)[i].x, &(* cartesian_data)[i].y, &(* cartesian_data)[i].z, &(* cartesian_data)[i].weight) == 4) {
      if (debug)
	printf ("i: %d x: %.15lf y: %.15lf z: %.15lf weight: %.15lf\n", i, (* cartesian_data)[i].x, (* cartesian_data)[i].y, (* cartesian_data)[i].z, (* cartesian_data)[i].weight);
      i++;
    }
  }

  return lines;
}

int lebedev_to_cartesian (lebedev ** in, cartesian ** out, int npoints) {
 
  int i;
  
  if (npoints > 0) {
    if (!(*out = malloc (sizeof(cartesian)*npoints)))
      return 0;
  }
  else
    return 0;
  
  for (i = 0; i < npoints; i++) {

    (* out)[i].x = sin ((* in)[i].theta) * cos ((* in)[i].phi);
    (* out)[i].y = sin ((* in)[i].theta) * sin ((* in)[i].phi);;
    (* out)[i].z = cos ((* in)[i].theta);
    (* out)[i].weight = (* in)[i].weight;
  }

  return npoints;
}

int cartesian_to_lebedev (cartesian ** in, lebedev ** out, int npoints) {

  int i;
  
  if (npoints > 0) {
    if (!(*out = malloc (sizeof(lebedev)*npoints)))
      return 0;
  }
  else
    return 0;
  
  for (i = 0; i < npoints; i++) {

    if (((* in)[i].x) > 0)
      (* out)[i].phi = atan(((* in)[i].y) / ((* in)[i].x));
    else if (((* in)[i].x) == 0)
      (((* in)[i].y) / ((* in)[i].y)) * (M_PI / 2);
    else if (((* in)[i].x) < 0 && ((* in)[i].y) >= 0)
      (* out)[i].phi = atan(((* in)[i].y) / ((* in)[i].x)) + M_PI;
    else if (((* in)[i].x) < 0 && ((* in)[i].y) < 0)
      (* out)[i].phi = atan(((* in)[i].y) / ((* in)[i].x)) - M_PI;

    (* out)[i].theta = acos(((* in)[i].z) / sqrt(((* in)[i].x) * ((* in)[i].x) + ((* in)[i].y) * ((* in)[i].y) + ((* in)[i].z) * ((* in)[i].z)));
    (* out)[i].weight = (* in)[i].weight;
  }

  return npoints;
}

/** 
 * usage: print usage information
 * 
 * @param cmdname - pointer to a string containing the name of the binary (argv[0])
 */
void usage(char * cmdname)
{
  printf("Usage: %s --input-file [input file] --output-file [output file] --to-cartesian --to-lebedev --debug.\n\n\
          Where:\n\
          \n\
          input-file - read input coordinate data from [input file]\n\
          output-file - write output coordinate data to [output file]\n\
          to-cartesian - treat input data as lebedev and convert to cartesian\n\
          to-lebedev - treat input data as cartesian and convert to lebedev\n\
          debug - output extended debug data\n\n", cmdname);
}

int main (int argc, char** argv) {

  int i, lines;
  int c, option_index;
  int to_cartesian, to_lebedev;
  cartesian * cartesian_data;
  lebedev * lebedev_data;
  FILE * input_file;
  FILE * output_file;
  char input_path [256];
  char output_path [256];

  memset (input_path, 0, 256);
  memset (output_path, 0, 256);
  
  lines = 0;
  to_cartesian = 0;
  to_lebedev = 0;

  static struct option options [] = {
    {"input-file", required_argument, NULL, 'i'},
    {"output-file", required_argument, NULL, 'o'},
    {"to-cartesian", no_argument, NULL, 'c'},
    {"to-lebedev", no_argument, NULL, 'l'},
    {"debug", no_argument, NULL, 'd'},
    {0, 0, 0, 0}
  };

  /* parse command line options using getopt_long */
  while ((c = getopt_long(argc, argv, "i:o:cld", options, &option_index)) != -1) {

    switch (c) {

    case 'i':
      strncpy (input_path, optarg, 255);
      break;

    case 'o':
      strncpy (output_path, optarg, 255);
      break;
      
    case 'c':
      if (!to_lebedev)
	to_cartesian= 1;
      break;
      
    case 'l':
      if (!to_cartesian)
	to_lebedev = 1;
      break;

    case 'd':
      debug = 1;
      break;

    default:
      usage (argv[0]);
      return 0;
    }
  }

  if ((strlen(input_path) > 0)){
    input_file = fopen (input_path, "r");
    
    if (!input_file) {
      printf ("Error opening file for input.\n");
      usage(argv[0]);
      return 0;
    }
  }
  else
    input_file = stdin;
  
  if ((strlen(output_path) > 0)){
    output_file = fopen (output_path, "w");
    
    if (!output_file) {
      printf ("Error opening file for output.\n");
      usage(argv[0]);
      return 0;
    }
  }
  else
    output_file = stdout;
  
  if(to_cartesian) {
    if (!(lines = read_lebedev_data(input_file, &lebedev_data))) {
      printf("Error parsing Lebedev input data.\n");
      return 0;
    }
    lebedev_to_cartesian(&lebedev_data, &cartesian_data, lines);
    for (i = 0; i < lines; i++)
      fprintf (output_file, "%lg %lg %lg %lg\n", cartesian_data[i].x, cartesian_data[i].y, cartesian_data[i].z, cartesian_data[i].weight);
  }
  else if (to_lebedev) {
    if (!(lines = read_cartesian_data(input_file, &cartesian_data))) {
      printf("Error parsing cartesian input data.\n");
      return 0;
    }
    cartesian_to_lebedev(&cartesian_data, &lebedev_data, lines);
      for (i = 0; i < lines; i++)
	fprintf (output_file, "%lg %lg %lg\n", lebedev_data[i].phi, lebedev_data[i].theta, lebedev_data[i].weight);
  }
  else {
    printf ("Either --to-cartesian or --to-lebedev must specified.\n");
    usage(argv[0]);
    return 0;
  }

  return 0;
}
