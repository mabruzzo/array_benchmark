#include <cstdio>
#include <cstdlib>

#include "benchmark.hpp"


#define BOOST_DISABLE_ASSERTS
#include <boost/multi_array.hpp>
#include "enzo_EnzoArray.hpp"
#include "enzo_EnzoArray_n.hpp"

using namespace boost::accumulators;


// Prepares the array
class pointer_ops{
public:
  typedef double* array_type;

  array_type prepare_array(long N){
    std::srand(RAND_SEED);
    double* array = new double[N*N*N];
    for (int i = 0; i<N*N*N; i++){
      array[i] = (double)std::rand();
    }
    return array;
  }

  inline double test_array(array_type array, long N){
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

  void cleanup(array_type *array){
    delete[] *array;
  }
};

void collect_array_benchmarks(long n_entries, int n_reps,
			      double results[2]){
  DEFAULT_FlUSH_SEED()
  accumulator_set<double, features<tag::mean, tag::variance>> acc;
  pointer_ops obj;

  for (int i=0; i<n_reps; i++){
    double* array = obj.prepare_array(n_entries);
    FLUSH_CACHE()
    double t1 = get_time();
    double val = obj.test_array(array, n_entries);
    double t2 = get_time();
    delete[] array;
    acc(t2-t1);
    printf("%.e,", val);
    fflush(stdout);
  }
  printf("current flush_seed = %u", flush_seed);
  results[0] = mean(acc);
  results[1] = std::sqrt(variance(acc));
}

// Prepares boost::multi_array

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

// also want to test what happens when we use the () access


void collect_multi_array_bracket_benchmarks(long n_entries, int n_reps,
					    double results[2]){
  DEFAULT_FlUSH_SEED()
  accumulator_set<double, features<tag::mean, tag::variance>> acc;

  for (int i=0; i<n_reps; i++){
    boost::multi_array<double, 3> array = prepare_multi_array(n_entries);
    FLUSH_CACHE()
    double t1 = get_time();
    double val = test_multi_array_bracket(array,n_entries);
    double t2 = get_time();
    acc(t2-t1);
    printf("%.e,", val);
    fflush(stdout);
  }
  printf("current flush_seed = %u", flush_seed);
  results[0] = mean(acc);
  results[1] = std::sqrt(variance(acc));
}

// Prepares boost::multi_array

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
				  double results[2]){
  DEFAULT_FlUSH_SEED()
  accumulator_set<double, features<tag::mean, tag::variance>> acc;

  for (int i=0; i<n_reps; i++){
    EnzoArray<double, 3> array = prepare_EnzoArray(n_entries);
    FLUSH_CACHE()
    double t1 = get_time();
    double val = test_EnzoArray(array,n_entries);
    double t2 = get_time();
    acc(t2-t1);
    printf("%.e,", (val));
    fflush(stdout);
  }
  printf("current flush_seed = %u", flush_seed);
  results[0] = mean(acc);
  results[1] = std::sqrt(variance(acc));
}


////////////////////////////////////////////////////////////////////////




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
  double results[2];
  // C style array
  //collect_array_benchmarks(n_entries, n_reps, results);
  pointer_ops obj;
  collect_benchmark(n_entries, n_reps, results, obj);
  fprintf(outfile,"%.15e, %.15e\n",results[0],results[1]);
  // boost::multi_array
  collect_multi_array_bracket_benchmarks(n_entries, n_reps, results);
  fprintf(outfile,"%.15e, %.15e\n",results[0],results[1]); 
  // EnzoArray
  collect_EnzoArray_benchmarks(n_entries, n_reps, results);
  fprintf(outfile,"%.15e, %.15e\n",results[0],results[1]); 
    
  fclose(outfile); 
  

  return 0;
  
}
