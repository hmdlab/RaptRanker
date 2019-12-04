/*
 * Main.cpp
 * Copyright (c) 2011 Yasuo Tabei All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE and * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include "SketchSort.hpp"

#include <iostream>
#include <cstdlib>

/* Globals */
void usage();
void version();
void parse_parameters (int argc, char **argv);

char *fname, *oname;
int   hamDist      = 1;
int   numblocks    = 4;
int   numchunks    = 3;
float cosDist      = 0.01;
bool  autoFlag     = false;
float missingratio = 0.0001;
bool  centering    = false;

int main(int argc, char **argv) 
{
  version();

  parse_parameters(argc, argv);

  SketchSort sketchsort;
  sketchsort.run(fname, oname, numblocks, hamDist, cosDist, numchunks, autoFlag, missingratio, centering);

  return 0;
}

void version(){
  std::cerr << "SketchSort version 0.0.8" << std::endl
            << "Written by Yasuo Tabei" << std::endl << std::endl;
}

void usage(){
  std::cerr << std::endl
       << "Usage: sketchsort [OPTION]... INFILE OUTFILE" << std::endl << std::endl
       << "       where [OPTION]...  is a list of zero or more optional arguments" << std::endl
       << "             INFILE       is the name of an input file" << std::endl
       << "             OUTFILE      is the name of an output file" << std::endl << std::endl
       << "Additional arguments (input and output files may be specified):" << std::endl
       << "       -hamdist [maximum hamming distance]" << std::endl
       << "       (default: " << hamDist << ")" << std::endl
       << "       -numblocks [the number of blocks]" << std::endl
       << "       (default: " << numblocks << ")" << std::endl
       << "       -cosdist [maximum cosine distance]" << std::endl
       << "       (default: " << cosDist << ")" << std::endl
       << "       -numchunks [the number of chunks]" << std::endl
       << "       (default: " << numchunks << ")" << std::endl
       << "       -auto " << std::endl
       << "       -missingratio " << std::endl
       << "       (default: " << missingratio << ")" << std::endl
       << "       -centering" << std::endl
       << std::endl;
  exit(0);
}

void parse_parameters (int argc, char **argv){
  if (argc == 1) usage();
  int argno;
  for (argno = 1; argno < argc; argno++){
    if (argv[argno][0] == '-'){
      if      (!strcmp (argv[argno], "-version")){
	version();
      }
      else if (!strcmp (argv[argno], "-auto")) {
	autoFlag = true;
      }
      else if (!strcmp (argv[argno], "-centering")) {
	centering = true;
      }
      else if (!strcmp (argv[argno], "-numblocks")) {
	if (argno == argc - 1) std::cerr << "Must specify minimum support after -numblocks" << std::endl;
	numblocks = atoi(argv[++argno]);
      }
      else if (!strcmp (argv[argno], "-hamdist")) {
	if (argno == argc - 1) std::cerr << "Must specify hamming distance threshold after -hamdist" << std::endl;
	hamDist = atoi(argv[++argno]);
      }
      else if (!strcmp (argv[argno], "-cosdist")) {
	if (argno == argc - 1) std::cerr << "Must specify cosine distance threshold size after -cosdist" << std::endl;
	cosDist = atof(argv[++argno]);
      }
      else if (!strcmp (argv[argno], "-numchunks")) {
	if (argno == argc - 1) std::cerr << "Must specify number of chunks after -numchunks" << std::endl;
	numchunks = atoi(argv[++argno]);
      }
      else if (!strcmp (argv[argno], "-missingratio")) {
	if (argno == argc - 1) std::cerr << "Must specify missing edge ratio after -missingratio" << std::endl;
	missingratio = atof(argv[++argno]);
      }
      else {
	usage();
      }
    } else {
      break;
    }
  }
  if (argno > argc)
    usage();

  fname = argv[argno];
  oname = argv[argno + 1];
}
