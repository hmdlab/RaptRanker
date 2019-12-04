/*
 * SketchSort.hpp
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

#ifndef SKETCHSORT_HPP
#define SKETCHSORT_HPP

// for boost
#include <boost/pool/pool.hpp>
#include <boost/pool/object_pool.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/hermitian.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>


#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <cmath>
#include <cstring>
#include <iterator>
#include <fstream>
#include <strstream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <set>
#include <stdint.h>
#include <time.h>

#include <boost/random.hpp>

using namespace boost::numeric::ublas;  // boost::numeric::ublas
using namespace boost::numeric;         // boost::numeric

#ifndef M_PI
#define M_PI 3.14159265358979323846/* pi */
#define M_PI_2 1.57079632679489661923/* pi/2 */
#define M_PI_4 0.78539816339744830962/* pi/4 */
#endif

struct pstat {
  int start;
  int end;
};

struct params {
  unsigned int numblocks;
  unsigned int numchunks;
  unsigned int projectDim;
  unsigned int chunk_dist;
  unsigned int chunks;
  unsigned int num_seq;
  unsigned int seq_len;
  unsigned int chunk_len;
  unsigned int start_chunk;
  unsigned int end_chunk;
  unsigned int cchunk;
  unsigned int *counter;
  pstat        *pos;
  pstat        *pchunks;
  float         cosDist;
  std::vector<unsigned int> ids;
  std::vector<int> blocks;
  std::ostream *os;
};

class SketchSort {
  boost::pool<> *p;
  std::vector<boost::numeric::ublas::vector<float> > fvs;
  std::vector<float> norms;
  uint8_t num_char;
  unsigned int dim;

  uint64_t numSort;
  uint64_t numHamDist;
  uint64_t numCosDist;

  void readFeature(const char *fname);
  void centeringData();
  void preComputeNorms();
  int projectVectors(unsigned int projectDim, std::vector<uint8_t*> &sig, params &param);
  void report(std::vector<uint8_t*> &sig, int l, int r, params &param);
  void sort(std::vector<uint8_t*> &sig, int spos, int epos, int l, int r, params &param);
  void radixsort(std::vector<uint8_t*> &sig, int spos, int epos, int l, int r, params &param);
  void insertionSort(std::vector<uint8_t*> &sig, int spos, int epos, int l, int r, params &param);
  bool calc_chunk_hamdist(uint8_t *seq1, uint8_t *seq2, const params &param);
  bool calc_hamdist(uint8_t *seq1, uint8_t *seq2, const params &param);
  bool check_canonical(uint8_t *seq1, uint8_t *seq2, const params &param);
  bool check_chunk_canonical(uint8_t *seq1, uint8_t *seq2, const params &param);
  float checkCos(unsigned int id1, unsigned int id2);
  double calcMissingEdgeRatio(params &param);
  void classify(std::vector<uint8_t*> &sig, int spos, int epos, int l, int r, int bpos, params &param);
  void multi_classification(std::vector<uint8_t*> &sig, int maxind, int l, int r, params &param);
  void refinement();
  void insertKnnList(unsigned int from_id, unsigned int to_id, float cosDist);
  void decideParameters(float _missingratio, params &param);
 public:
  SketchSort() {};
  void run(const char *fname, const char *oname,
	   unsigned int _numblocks, 
           unsigned int _dist,
	   float        _cosDist,
	   unsigned int _numchunks, 
	   bool         _autoFlag,
	   float        _missingratio, 
	   bool         _centering);

};

#endif // SKETCHSORT_HPP
