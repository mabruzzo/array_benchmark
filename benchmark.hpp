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

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
using namespace boost::accumulators;

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


template<typename S>
void collect_benchmark(long &n_entries, int &n_reps, double *results, S obj){
  DEFAULT_FlUSH_SEED()
  accumulator_set<double, features<tag::mean, tag::variance>> acc;

  for (int i=0; i<n_reps; i++){
    typename S::array_type array;
    array = obj.prepare_array(n_entries);
    FLUSH_CACHE()
    double t1 = get_time();
    double val = obj.test_array(array, n_entries);
    double t2 = get_time();
    obj.cleanup(&array);
    acc(t2-t1);
    printf("%.e,", val);
    fflush(stdout);
  }
  printf("current flush_seed = %u", flush_seed);
  results[0] = mean(acc);
  results[1] = std::sqrt(variance(acc));
}
