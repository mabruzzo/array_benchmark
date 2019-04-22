#include <cstdio>
#include <cstdlib>

#define BOOST_DISABLE_ASSERTS
#include <boost/multi_array.hpp>
#include "enzo_EnzoArray.hpp"
#include "enzo_EnzoArray_n.hpp"

// https://stackoverflow.com/questions/2349776/how-can-i-benchmark-c-code-easily
#include <sys/time.h>
#include <sys/resource.h>

inline double get_time()
{
  struct timeval t;
  struct timezone tzp;
  gettimeofday(&t, &tzp);
  return t.tv_sec + t.tv_usec*1e-6;
}

// For flushing the cache: https://stackoverflow.com/a/34461372
// Sizes of L1 & L2 +100 caches divided by 4 - I don't think this accomplishes anything
//#define BIGGER_THAN_CACHE_SIZES (32 + 256 + 100)*24
#define BIGGER_THAN_CACHE_SIZES (10*1024*1024)
#define RAND_SEED 97

static unsigned int flush_seed;

#define DEFAULT_FlUSH_SEED() flush_seed = BIGGER_THAN_CACHE_SIZES;           \

#define FLUSH_CACHE()                                                        \
  std::srand(flush_seed);                                                    \
  long *p = new long[BIGGER_THAN_CACHE_SIZES];                               \
  for(int i = 0; i < BIGGER_THAN_CACHE_SIZES; i++)                           \
  {                                                                          \
    p[i] = std::rand();				                             \
  }                                                                          \
  long seed_val = p[std::rand() % (BIGGER_THAN_CACHE_SIZES +1)];             \
  flush_seed = (unsigned int)((seed_val) % 1000);                            \ 
  delete[] p;




double* prepare_array(long N){
  std::srand(RAND_SEED);
  double* array = new double[N*N*N];
  for (int i = 0; i<N*N*N; i++){
    array[i] = (double)std::rand();
  }
  return array;
}

inline double test_array(double* array, long N){
  double val=0;
  for (int iz = 0; iz<N; iz++){
    for (int iy = 0; iy<N; iy++){
      for (int ix = 0; ix<N; ix++){
	val+=array[(iz*N+iy)*N +ix];
      }
    }
  }
  return val;
}

void collect_array_benchmarks(long n_entries, int n_reps,
			      double* results){
  DEFAULT_FlUSH_SEED()
  for (int i=0; i<n_reps; i++){
    double* array = prepare_array(n_entries);
    FLUSH_CACHE()
    double t1 = get_time();
    double val = test_array(array, n_entries);
    double t2 = get_time();
    delete[] array;
    results[i] = t2-t1;
    printf("%.e,", val);
    fflush(stdout);
  }
  printf("current flush_seed = %u\n", flush_seed);
}

////////////////////////////////////////////////////////////////////////////
// boost::multi_array using [][][] access


inline boost::multi_array<double,3> prepare_multi_array(long N){
  std::srand(RAND_SEED);
  boost::multi_array<double, 3> array{boost::extents[N][N][N]};
  for (int iz = 0; iz<N; iz++){
    for (int iy = 0; iy<N; iy++){
      for (int ix = 0; ix<N; ix++){
	array[iz][iy][ix] = (double)std::rand();
      }
    }
  }
  return array;
}


inline double test_multi_array_bracket(boost::multi_array<double, 3> &array,
				     long N){
  double val = 0;
  for (int iz = 0; iz<N; iz++){
    for (int iy = 0; iy<N; iy++){
      for (int ix = 0; ix<N; ix++){
	val+=array[iz][iy][ix];
      }
    }
  }
  return val;
}

////////////////////////////////////////////////////////////////////////////
// boost::multi_array using () access

void collect_multi_array_bracket_benchmarks(long n_entries, int n_reps,
					    double* results){
  DEFAULT_FlUSH_SEED()

  for (int i=0; i<n_reps; i++){
    boost::multi_array<double, 3> array = prepare_multi_array(n_entries);
    FLUSH_CACHE()
    double t1 = get_time();
    double val = test_multi_array_bracket(array,n_entries);
    double t2 = get_time();
    results[i] = t2-t1;
    printf("%.e,", val);
    fflush(stdout);
  }
  printf("current flush_seed = %u\n", flush_seed);
}

// Prepares boost::multi_array with container

inline double test_multi_array_paren(boost::multi_array<double, 3> &array,
				     long N){
  double val = 0;
  std::array<int,3> idx = {{0,0,0}};
  for (int iz = 0; iz<N; iz++){
    idx[0] = iz;
    for (int iy = 0; iy<N; iy++){
      idx[1] = iy;
      for (int ix = 0; ix<N; ix++){
	idx[2] = ix;
	val+=array(idx);
      }
    }
  }
  return val;
}

void collect_multi_array_bracket_paren(long n_entries, int n_reps,
				       double* results){
  DEFAULT_FlUSH_SEED()

  for (int i=0; i<n_reps; i++){
    boost::multi_array<double, 3> array = prepare_multi_array(n_entries);
    FLUSH_CACHE()
    double t1 = get_time();
    double val = test_multi_array_paren(array,n_entries);
    double t2 = get_time();
    results[i] = t2-t1;
    printf("%.e,", val);
    fflush(stdout);
  }
  printf("current flush_seed = %u\n", flush_seed);
}

/////////////////////////////////////////////////////////////////////////
// EnzoArray

inline EnzoArray<double,3> prepare_EnzoArray(long N){
  std::srand(RAND_SEED);
  EnzoArray<double, 3> array(N,N,N);
  for (int iz = 0; iz<N; iz++){
    for (int iy = 0; iy<N; iy++){
      for (int ix = 0; ix<N; ix++){
	array(iz,iy,ix) = (double)std::rand();
      }
    }
  }
  return array;
}


inline double test_EnzoArray(EnzoArray<double, 3> &array,long N){
  double val = 0;
  for (int iz = 0; iz<N; iz++){
    for (int iy = 0; iy<N; iy++){
      for (int ix = 0; ix<N; ix++){
	val += array(iz,iy,ix);
      }
    }
  }
  return val;
}


// CHECK MY ENZO CODE FOR std::abs where it should be std::absf!!!!

void collect_EnzoArray_benchmarks(long n_entries, int n_reps,
				  double* results){
  DEFAULT_FlUSH_SEED()

  for (int i=0; i<n_reps; i++){
    EnzoArray<double, 3> array = prepare_EnzoArray(n_entries);
    FLUSH_CACHE()
    double t1 = get_time();
    double val = test_EnzoArray(array,n_entries);
    double t2 = get_time();
    results[i] = t2-t1;
    printf("%.e,", (val));
    fflush(stdout);
  }
  printf("current flush_seed = %u\n", flush_seed);
}

/////////////////////////////////////////////////////////////////////////
// EnzoArray with negative indexing support

inline EnzoArrayN<double,3> prepare_EnzoArrayN(long N){
  std::srand(RAND_SEED);
  EnzoArrayN<double, 3> array(N,N,N);
  for (int iz = 0; iz<N; iz++){
    for (int iy = 0; iy<N; iy++){
      for (int ix = 0; ix<N; ix++){
	array(iz,iy,ix) = (double)std::rand();
      }
    }
  }
  return array;
}


inline double test_EnzoArrayN(EnzoArrayN<double, 3> &array,long N){
  double val = 0;
  for (int iz = 0; iz<N; iz++){
    for (int iy = 0; iy<N; iy++){
      for (int ix = 0; ix<N; ix++){
	val += array(iz,iy,ix);
      }
    }
  }
  return val;
}

void collect_EnzoArrayN_benchmarks(long n_entries, int n_reps,
				   double* results){
  DEFAULT_FlUSH_SEED()

  for (int i=0; i<n_reps; i++){
    EnzoArrayN<double, 3> array = prepare_EnzoArrayN(n_entries);
    FLUSH_CACHE()
    double t1 = get_time();
    double val = test_EnzoArrayN(array,n_entries);
    double t2 = get_time();
    results[i] = t2-t1;
    printf("%.e,", (val));
    fflush(stdout);
  }
  printf("current flush_seed = %u\n", flush_seed);
}


void save_results(FILE *file, int n_reps, double *results){
  for (int i = 0; i<n_reps-1; i++){
    fprintf(file,"%.15e, ",results[i]);
  }
  fprintf(file,"%.15e\n",results[n_reps-1]);
}


int main(int argc, char *argv[])
{
  if (argc != 4){
    printf("Need to provide of number of entries per dimension, "
	   "number of repetitions and file name.\n");
    return 1;
  }

  long n_entries = atol(argv[1]);
  int n_reps = atoi(argv[2]);

  if ( (n_reps <= 0) || (n_entries <= 0) ){
    printf("number of entries per dimension and number of repetitions must be"
	   " >0.\n");
    return 1;
  }

  FILE * outfile = fopen(argv[3],"w");
  double* results = new double[n_reps];
  // C style array
  collect_array_benchmarks(n_entries, n_reps, results);
  save_results(outfile, n_reps, results);
  // boost::multi_array
  collect_multi_array_bracket_benchmarks(n_entries, n_reps, results);
  save_results(outfile, n_reps, results);
  collect_multi_array_bracket_paren(n_entries, n_reps, results);
  save_results(outfile, n_reps, results);
  // EnzoArray
  collect_EnzoArray_benchmarks(n_entries, n_reps, results);
  save_results(outfile, n_reps, results);

  // EnzoArray negative indexing
  collect_EnzoArrayN_benchmarks(n_entries, n_reps, results);
  save_results(outfile, n_reps, results);

  fclose(outfile);

  delete[] results;
  

  return 0;
  
}
