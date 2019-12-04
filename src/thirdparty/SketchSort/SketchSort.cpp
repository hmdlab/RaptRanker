/*
 * SketchSort.cpp
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

template<class T>
inline uint8_t sign(T val) {
  if (val > 0)
    return 1;
  return 0;
}

template<class T>
inline T max(T a1, T a2) {
  if (a1 > a2)
    return a1;
  return a2;
}

bool cmp(const std::pair<int, float> &p1, const std::pair<int, float> &p2) {
  return p1.second < p2.second;
}

void SketchSort::readFeature(const char *fname) {
  std::ifstream ifs(fname);

  if (!ifs) {
    std::cerr << "can not open " << fname << std::endl;
    exit(0);
  }

  dim             = 0;
  float val       = 0.f;
  std::string line;
  while (std::getline(ifs, line)) {
    fvs.resize(fvs.size() + 1);
    boost::numeric::ublas::vector<float> &fv = fvs[fvs.size() - 1];
    uint32_t counter = 0;
    std::istringstream is(line);
    if (dim != 0) {
      fv.resize(dim);
      while (is >> val) {
	fv[counter++]= val;
      }
      if (counter !=  dim) {
	std::cerr << "dimesions of the input vector should be same!" << std::endl;
	std::cerr << line << std::endl;
	std::cerr << "dim:" << dim << " dim:" << counter << std::endl;
	exit(1);
      }
    }
    else {
      while (is >> val) {
	fv.resize(counter + 1);
	fv[counter] = val;
	counter++;
      }
      dim = counter;
    }
  }
}

void SketchSort::centeringData() {
  size_t dim     = fvs[0].size();
  size_t numData = fvs.size();
  float  mean;
  for (size_t i = 0; i < dim; i++) {
    mean = 0.f;
    for (size_t j = 0; j < numData; j++) {
      mean += fvs[j][i];
    }
    mean /= (float)numData;
    for (size_t j = 0; j < numData; j++) {
      fvs[j][i] -= mean;
    }
  }
}

/* sparce random projection
int SketchSort::projectVectors(unsigned int projectDim, std::vector<uint8_t*> &sig, params &param) {

  p = new boost::pool<>(sizeof(uint8_t));
  sig.resize(fvs.size());
  param.ids.resize(fvs.size());
  for (size_t i = 0; i < sig.size(); i++) {
    //    sig[i]    = new uint32_t[projectDim + 1];                                      
    sig[i]    = (uint8_t*)p->ordered_malloc(projectDim + 1);
    param.ids[i] = i;
  }

  boost::mt19937 gen(static_cast<unsigned long>(time(0)));
  boost::uniform_real<> dst(0.f, 1.f);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> > rand(gen, dst);
  //  double tiny = 1.0/1.79e+308;                                           
  std::vector<std::pair<int, float> > randMat;
  float s = sqrt(float(dim));
  //  float s     = dim/log(dim);
  float thr   = 1.f/(2*s);
  float coff  = sqrt(s);
  for (size_t i = 0; i < projectDim; i++) {
    randMat.clear();
    for (size_t j = 0; j < dim; j++) {
      float r   = rand();
      if       (r < thr) {
        randMat.push_back(std::make_pair(j, coff));
      } else if (r < 2*thr) {	
        randMat.push_back(std::make_pair(j, -coff));
      }
    }

    for (size_t j = 0; j < fvs.size(); j++) {
      boost::numeric::ublas::vector<float> &fv  = fvs[j];
      double proc = 0.f;
      for (size_t k = 0; k < randMat.size(); k++) {
        proc += fv[randMat[k].first] * randMat[k].second;
      }
      sig[j][i+1] = sign(proc);
    }
  }
  param.seq_len = projectDim;
  param.num_seq = fvs.size();

  return 1;
}
*/

int SketchSort::projectVectors(unsigned int projectDim, std::vector<uint8_t*> &sig, params &param) {
  std::vector<float> randMat;
  p = new boost::pool<>(sizeof(uint8_t));
  sig.resize(fvs.size());
  param.ids.resize(fvs.size());
  for (size_t i = 0; i < sig.size(); i++) {
    //    sig[i]    = new uint32_t[projectDim + 1];
    sig[i]    = (uint8_t*)p->ordered_malloc(projectDim + 1);
    param.ids[i] = i;
  }

  boost::mt19937 gen(static_cast<unsigned long>(time(0)));
  boost::normal_distribution<> dst(0.f, 1.f);
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > rand(gen, dst);

  //  double tiny = 1.0/1.79e+308;
  randMat.resize(dim + 1);
  for (size_t i = 0; i < projectDim; i++) {
    for (size_t j = 0; j <= dim; j++) {
      randMat[j] = rand();
    }

    for (size_t j = 0; j < fvs.size(); j++) {
      boost::numeric::ublas::vector<float> &fv  = fvs[j];
      double proc = 0.f;
      for (size_t k = 0; k < fv.size(); k++)
        proc += fv[k] * randMat[k];

      sig[j][i+1] = sign(proc);
    }
  }
  param.seq_len = projectDim;
  param.num_seq = fvs.size();

  return 1;
}

inline float SketchSort::checkCos(unsigned int id1, unsigned int id2) {
  ++numCosDist;
  boost::numeric::ublas::vector<float> &fv_1 = fvs[id1];
  boost::numeric::ublas::vector<float> &fv_2 = fvs[id2];
  float sum = boost::numeric::ublas::inner_prod(fv_1, fv_2);

  return (1.f - sum*(norms[id1]*norms[id2]));
}

inline void SketchSort::sort(std::vector<uint8_t*> &sig, int spos, int epos, int l, int r, params &param) {
   if (r - l + 1 > 50) radixsort(sig, spos, epos, l, r, param);
  else                insertionSort(sig, spos, epos, l, r, param);
}

inline void SketchSort::radixsort(std::vector<uint8_t*> &sig, int spos, int epos, int l, int r, params &param) {
  unsigned int *c      = param.counter;
  std::vector<unsigned int> &ids      = param.ids;
  std::vector<uint8_t*> newsig(r - l + 1);
  std::vector<unsigned int> newids(r - l + 1);
  unsigned int tmp;
  int tpos = spos - 1;
  while (++tpos <= epos) {
    for (int i = 0; i < num_char; i++) *(c + i) = 0;
    for (int i = l; i <= r; i++) c[sig[i][tpos]]++;
    for (int i = 1; i < num_char; i++) *(c + i) += *(c + i - 1);
    for (int i = r; i >= l; --i) {
      tmp = --c[sig[i][tpos]] + l;
      newids[tmp - l] = ids[i];
      newsig[tmp - l] = sig[i];
    }
    if (++tpos <= epos) {
      for (int i = 0; i < num_char; i++) *(c + i) = 0;
      for (int i = l; i <= r; i++) c[newsig[i - l][tpos]]++;
      for (int i = 1; i < num_char; i++) *(c + i) += *(c + i - 1);
      for (int i = r; i >= l; --i) {
	tmp = --c[newsig[i - l][tpos]] + l;
	ids[tmp] = newids[i - l];
	sig[tmp] = newsig[i - l];
      }
    }
    else {
      for (int i = l; i <= r; i++) {
	ids[i] = newids[i - l];
	sig[i] = newsig[i - l];
      }
      return;
    }
  }
}

inline void SketchSort::insertionSort(std::vector<uint8_t*> &sig, int spos, int epos, int l, int r, params &param) {
  int i, j;
  uint8_t *pivot, pval;
  unsigned int pid;
  std::vector<unsigned int> &ids = param.ids;
  for (int tpos = spos; tpos <= epos; tpos++) {
    for (i = l + 1; i <= r; i++) {
      pivot = sig[i]; pval = sig[i][tpos]; pid = ids[i];
      for (j = i; j > l && sig[j-1][tpos] > pval; j--) {
	sig[j]       = sig[j-1];
	ids[j]       = ids[j-1];
      }
      sig[j] = pivot;
      ids[j] = pid;
    }
  }
}

inline void SketchSort::classify(std::vector<uint8_t*> &sig, int spos, int epos, int l, int r, int bpos, params &param) {
  int n_l = l, n_r = r;
  for (int iter = l + 1; iter <= r; iter++) {
    if (!std::equal(sig[n_l] + spos, sig[n_l] + epos + 1, sig[iter] + spos)) {
      n_r = iter - 1;
      if (n_r - n_l >= 1)
	multi_classification(sig, bpos + 1, n_l, n_r, param);
      n_l = iter;
    }
  }
  if (r - n_l >= 1)
    multi_classification(sig, bpos + 1, n_l, r, param);
}

inline bool SketchSort::calc_chunk_hamdist(uint8_t *seq1, uint8_t *seq2, const params &param) {
  ++numHamDist;
  unsigned int d = 0;
  for (size_t i = 1;  i <= param.chunk_len; i++) 
    if (*seq1++ != *seq2++ && ++d > param.chunk_dist) return false;
  return true;
}

inline bool SketchSort::check_chunk_canonical(uint8_t *seq1, uint8_t *seq2, const params &param) {
  unsigned int d = 0;
  int end        = param.pchunks[param.cchunk].start - 1;
  int j          = 1;
  int tend       = param.pchunks[j].end;
  int i          = 0;
  
  while (++i <= end) {
    if ((d += abs(seq1[i] - seq2[i])) > param.chunk_dist) {
      while (++i <= tend) d += abs(seq1[i] - seq2[i]);
			    //	if (seq1[i] != seq2[i]) ++d;
      d = 0;
      tend = param.pchunks[++j].end;
      i    = param.pchunks[j].start - 1;
      continue;
    }
    if (tend == i)
      return false;
  }
  return true;
}

inline bool SketchSort::check_canonical(uint8_t *seq1, uint8_t *seq2, const params &param) {
  size_t sb = 1, eb = 1;
  size_t b;
  for (size_t i = 0, size = param.blocks.size(); i < size; i++) {
    eb = param.blocks[i];
    for (b = sb; b < eb; b++) {
      if (std::equal(seq1 + param.pos[b].start, seq1 + param.pos[b].end + 1, seq2 + param.pos[b].start))
	  return false;
    }
    sb = param.blocks[i] + 1;
  }
  return true;
}


inline void SketchSort::report(std::vector<uint8_t*> &sig, int l, int r, params &param) {
  //  std::cout << "report" << std::endl;
  float cosDist;
  for (int i = l; i < r; i++) {
    for (int j = i + 1; j <= r; j++) {
      if (check_canonical(sig[i], sig[j], param) && 
	  calc_chunk_hamdist(sig[i] + param.start_chunk, sig[j] + param.start_chunk, param) && 
	  check_chunk_canonical(sig[i], sig[j], param) &&
	  ((cosDist = checkCos(param.ids[i], param.ids[j])) <= param.cosDist)) {
	(*param.os) << param.ids[i] << " " << param.ids[j] << " " << cosDist << std::endl;
      }
    }
  }
}

void SketchSort::multi_classification(std::vector<uint8_t*> &sig, int maxind, int l, int r, params &param) {
  if (param.blocks.size() == param.numblocks - param.chunk_dist) {
    report(sig, l, r, param);
    return;
  }

  for (int bpos = maxind; bpos <= (int)param.numblocks; bpos++) {
    if (param.blocks.size() + (param.numblocks - bpos + 1) < param.numblocks - param.chunk_dist) { // pruning
      //      std::cerr << "return " << std::endl;
      return;
    }
    param.blocks.push_back(bpos);
    sort(sig, param.pos[bpos].start, param.pos[bpos].end, l, r, param);
    classify(sig, param.pos[bpos].start, param.pos[bpos].end, l, r, bpos, param);
    param.blocks.pop_back();
  }
}

double combination(int n, int m) {
  double sum = 1.0;
  for (int i = 0; i < m; i++) {
    sum *= (n-i)/(m-i);
  }
  return sum;
}

double SketchSort::calcMissingEdgeRatio(params &param) {
  double sum = 0.f;
  double prob = acos(1.0 - param.cosDist)/M_PI;
  for (unsigned int k = 0; k <= param.chunk_dist; k++) {
    sum += (combination(param.projectDim, k) * pow(prob, k) * pow(1 - prob, param.projectDim - k));
  }
  return pow(1.0 - sum, param.numchunks);
}

void SketchSort::preComputeNorms() {
  norms.resize(fvs.size());
  float sum;
  for (size_t i = 0; i < fvs.size(); i++) {
    boost::numeric::ublas::vector<float> &fv = fvs[i];
    sum = 0.f;
    for (size_t j = 0; j < fv.size(); j++) {
      sum += pow(fv[j], 2);
    }
    norms[i] = 1.f/sqrt(sum);
  }
}

void SketchSort::decideParameters(float _missingratio, params &param) {
  unsigned int hamDist   = 1;
  unsigned int numBlocks = hamDist + 3;
  unsigned int numchunks = 0;

  do {
    if (numchunks > 30) {
      hamDist   += 1;
      numBlocks  = hamDist + 3;
      numchunks  = 0;
    }
    numchunks += 1;
    param.chunk_dist = hamDist;
    param.numblocks  = numBlocks;
    param.numchunks  = numchunks;
  } while (calcMissingEdgeRatio(param) >= _missingratio);
}

void SketchSort::run(const char *fname, const char *oname, 
		  unsigned int _numblocks, 
		  unsigned int _dist, 
		  float        _cosDist,
		  unsigned int _numchunks, 
		  bool         _autoFlag,
		  float        _missingratio, 
		  bool         _centering)
{
  params param;
  param.numblocks  = _numblocks;
  param.numchunks  = _numchunks;
  param.chunk_dist = _dist;
  param.cosDist    = _cosDist;
  num_char         = 2;
  param.projectDim = 32;

  numSort    = 0;
  numCosDist = 0;
  numHamDist = 0;
  
  if (_autoFlag) {
    std::cerr << "deciding parameters such that the missing edge ratio is no more than " << _missingratio << std::endl;
    decideParameters(_missingratio, param);
    std::cout << "decided parameters:" << std::endl;
    std::cout << "hamming distance threshold: " << param.chunk_dist << std::endl;
    std::cout << "number of blocks: " << param.numblocks << std::endl;
    std::cout << "number of chunks: "  << param.numchunks << std::endl;
    std::cout << std::endl;
  }

  std::ofstream ofs(oname);
  param.os         = &ofs;

  std::cout << "missing edge ratio:" << calcMissingEdgeRatio(param) << std::endl;

  std::cerr << "start reading" << std::endl;
  double readstart = clock();
  readFeature(fname);
  double readend   = clock();
  std::cerr << "end reading" << std::endl;  
  std::cout << "readtime:" << (readend - readstart)/(double)CLOCKS_PER_SEC << std::endl;

  if (_centering) {
    std::cerr << "start making input-data centered at 0" << std::endl;
    double centeringstart = clock();
    centeringData();
    double centeringend = clock();
    std::cerr << "end making input-data centered at 0" << std::endl;
    std::cout << "centering time:" << (centeringend - centeringstart)/(double)CLOCKS_PER_SEC << std::endl;

  }


  double totalstart = clock(); 
  preComputeNorms();
  //param.projectDim = 2*(int)log(dim);

  param.counter = new unsigned int[num_char];
  
  std::cout << "number of data:" << fvs.size() << std::endl;
  std::cout << "data dimension:" << dim << std::endl;
  std::cout << "projected dimension:" << param.projectDim << std::endl;
  std::cout << "length of strings:" << param.projectDim * param.numchunks << std::endl;
  std::cout << "number of chunks:" << param.numchunks << std::endl;

  double projectstart = clock();
  std::cerr << "start projection" << std::endl;
  std::vector<uint8_t*> sig;
  projectVectors(param.projectDim * param.numchunks, sig, param);
  //read(fname, sig, param);
  std::cerr << "end projection" << std::endl;
  double projectend = clock();
  std::cout << "projecttime:" << (projectend - projectstart)/(double)CLOCKS_PER_SEC << std::endl;

  param.pchunks = new pstat[param.numchunks + 1];
  for (int i = 1; i <= (int)param.numchunks; i++) {
    param.pchunks[i].start = (int)ceil((double)param.seq_len*((double)(i - 1)/(double)param.numchunks)) + 1;
    param.pchunks[i].end   = (int)ceil((double)param.seq_len*(double)i/(double)param.numchunks);
  }

  double msmtime = 0.0;


  std::cerr << "chunk distance:" << param.chunk_dist << std::endl;
  std::cerr << "the number of blocks:" << param.numblocks << std::endl;
  param.pos = new pstat[param.numblocks + 1];
  for (int i = 1; i <= (int) param.numchunks; i++) {
    param.chunk_len   = param.pchunks[i].end - param.pchunks[i].start + 1;
    param.start_chunk = param.pchunks[i].start; 
    param.end_chunk   = param.pchunks[i].end;
    param.cchunk      = i;
    for (int j = 1; j <= (int)param.numblocks; j++) {
      param.pos[j].start = (int)ceil((double)param.chunk_len*((double)(j - 1)/(double)param.numblocks)) + param.pchunks[i].start;
      param.pos[j].end   = (int)ceil((double)param.chunk_len*(double)j/(double)param.numblocks) + param.pchunks[i].start - 1;
    }
    std::cerr << "start enumeration chunk no " << i << std::endl;
    double msmstart = clock();
    multi_classification(sig, 1, 0, param.num_seq - 1, param);
    double msmend   = clock();
    msmtime += (msmend - msmstart)/(double)CLOCKS_PER_SEC;
  }
  std::cout << "msmtime:" << msmtime << std::endl;

  double totalend = clock();
  std::cout << "cputime:" << (totalend - totalstart)/(double)CLOCKS_PER_SEC << std::endl;

  std::cout << "numSort:" << combination(param.numblocks, param.chunk_dist) * param.numchunks << std::endl;
  std::cout << "numHamDist:" << numHamDist << std::endl;
  std::cout << "numCosDist:" << numCosDist << std::endl;
  ofs.close();
  // destructor
  delete p;
  delete[] param.counter;
  delete[] param.pchunks;
  delete[] param.pos;

  return;
}
