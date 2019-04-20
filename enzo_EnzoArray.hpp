#ifndef ENZO_ENZO_ARRAY_HPP
#define ENZO_ENZO_ARRAY_HPP

#include <stdio.h>
#include <cstddef>
#include <type_traits>
#include <memory>

// EnzoArray
// Like arrays in Athena++, the indices are listed in order of increasing
//   access speed. Imagine a 3D array with shape {mz,my,mx} array(k,j,i) is
//   equivalent to accessing index ((k*my + j)*mx + i) of the pointer
// Dimensions are numbered with increasing indexing speed (dim0, dim1, ...)
// Should not be confused with the EnzoArray in the original Enzo. (Maybe it
//   it should be renamed CelloArray?)

// Going to define an indexing type (borrowing naming convention from numpy)
using intp = std::ptrdiff_t;

#define ASSERT(F,M,A) /* ... */
#define ASSERT1(F,M,A,A1) /* ... */
#define ASSERT2(F,M,A,A1,A2) /* ... */
#define ASSERT3(F,M,A,A1,A2,A3) /* ... */
#define ASSERT4(F,M,A,A1,A2,A3,A4) /* ... */


class ESlice
{
  // Represents a slice
public:
  // Exists only to allow for construction of arrays of slices
  ESlice() : start(0), stop(0) {}

  // Include values from start through stop-1
  ESlice(intp start, intp stop) : start(start), stop(stop)
  { ASSERT("ESlice", "start must be less than stop.", stop>start); }

  ESlice(int start, int stop) : start((intp)start), stop((intp)stop)
  { ASSERT("ESlice", "start must be less than stop.", stop>start); }

public:
  // Locations of start and stop values
  intp start;
  intp stop;
};

// Defining the macro called CHECK_BOUNDS macro means that the bounds of an
// EnzoArray, are checked everytime operator() is called
template<typename T>
bool check_bounds_(std::size_t *shape, T first) {return *shape > first;}
template<typename T, typename... Rest>
bool check_bounds_(std::size_t *shape, T first, Rest... rest){
  //(get's unrolled at compile time)
  return (*shape > first) && check_bounds_(++shape, rest...);
}

#ifdef CHECK_BOUNDS
#  define CHECK_BOUND3D(shape, k, j, i)                                       \
  ASSERT("FixedDimArray_","Invalid index", check_bounds_(shape,k,j,i));
#  define CHECK_BOUNDND(shape, ARGS)                                          \
  ASSERT("FixedDimArray_","Invalid index", check_bounds_(shape, ARGS...));
#else
#  define CHECK_BOUND3D(shape, k, j, i)   /* ... */
#  define CHECK_BOUNDND(shape, ARGS)      /* ... */
#endif

// To define EnzoArray to have arbitrary dimensions we need to accept a variable
// number of arguments to indicate shape during construction, to produce a
// subarray, and to access elements. We specify the number of dimensions of the
// array as a template argument and accept values with variadic template
// arguments to check that the appropriate the number of values are specified
// at compile time. In each case, we need to guaruntee that all arguments
// are a given type. The solution is based on:
//   - https://stackoverflow.com/a/28253503
//   - https://stackoverflow.com/a/31767710
template <bool...> struct bool_pack;
template <bool... vals> using all_true = std::is_same<bool_pack<true, vals...>,
						      bool_pack<vals..., true>>;
// May want to add more possible int types (e.g. long,short,etc.)
#define REQUIRE_TYPE(T,type1)                                                \
  class = std::enable_if<all_true<(std::is_same<T,type1>::value)...>::value>
#define REQUIRE_TYPE2(T,type1,type2)					     \
  class = std::enable_if<all_true<(std::is_same<T,type1>::value ||           \
				   std::is_same<T,type2>::value)...>::value>
#define REQUIRE_INT(T) REQUIRE_TYPE2(T,intp,int)


// Helper function to compute pointer index
// T is expected to be intp or int
template<typename T>
intp calc_index_(intp* stride, T first){return first;}

template<typename T>
intp calc_index_(intp* stride, T first, T last){
  return (*stride)*first + last;
}

template<typename T, typename... Rest>
intp calc_index_(intp* stride, T first, Rest... rest){
  // gets unrolled at compile time
  return (*stride)*first + calc_index_(++stride, rest...);
}

// The following 2 functions are helpful while dynamically iterating
template<typename T>
intp calc_index_(const std::size_t D, const intp offset,
		 const intp* stride, const T* indices){
  std::size_t out = offset + indices[D-1];
  for (std::size_t i = 0; i+1<D;i++){
    out += stride[i]*indices[i];
  }
  return out;
}

// outer means indices for "outer" loops
inline void increment_outer_indices_(std::size_t D, intp *indices, intp *shape,
				     bool &continue_outer_iter){
  std::size_t i = D-1;
  while (0 != (i--)){
    indices[i]+=1;
    if (indices[i] != shape[i]){
      return;
    } else if (i > 0){
      indices[i] = 0;
    }
  }
  continue_outer_iter = false;
}


template<typename T>
class dataWrapper
{
  // tracks the underlying EnzoArray pointer and ownership of it
public:

  dataWrapper(T* data, bool owns_ptr) : data_(data), owns_ptr_(owns_ptr) { };
  ~dataWrapper() {if (owns_ptr_) { delete[] data_; }}
  T* get() const noexcept { return data_; }

private:
  T *data_;
  bool owns_ptr_;
};



template<typename T, std::size_t D>
class EnzoArray;

template<typename T, std::size_t D>
class TempArray_;

template<typename T, std::size_t D>
class FixedDimArray_
{
public:

  // Destructor
  ~FixedDimArray_() { cleanup_helper_();}

  // operator() - used to access array Elements
  template<typename... Args, REQUIRE_INT(Args)>
  T &operator() (Args... args) {
    static_assert(D==sizeof...(args),
		  "Number of indices don't match number of dimensions");
    CHECK_BOUNDND(shape, args)
    return data_[offset_ + calc_index_(stride_,args...)];
  }
  template<typename... Args, REQUIRE_INT(Args)>
  T operator() (Args... args) const {
    static_assert(D==sizeof...(args),
		  "Number of indices don't match number of dimensions");
    CHECK_BOUNDND(shape, args)
    return data_[offset_ + calc_index_(stride_,args...)];
  }

  // Specialized implementation for 3D arrays (reduces compile time)
  T &operator() (const int k, const int j, const int i){
    static_assert(D==3, "3 indices should only be specified for 3D arrays");
    CHECK_BOUND3D(shape, k, j, i)
    return data_[offset_ + k*stride_[0] + j*stride_[1] + i];
  }
  T operator() (const int k, const int j, const int i) const{
    static_assert(D==3, "3 indices should only be specified for 3D arrays");
    CHECK_BOUND3D(shape, k, j, i)
    return data_[offset_ + k*stride_[0] + j*stride_[1] + i];
  }


  // Produce a copy of the array.
  TempArray_<T,D> deepcopy();
  
  // Returns a subarray with same D. Expects an instance of ESlice for each
  // dimension. For a 3D EnzoArray, the function declaration might look like:
  // EnzoArray<T,3> subarray(ESlice k_slice, ESlice j_slice, ESlice i_slice);
  //
  // The Special Case of no arguments returns the full array
  template<typename... Args, REQUIRE_TYPE(Args,ESlice)>
  TempArray_<T,D> subarray(Args... args);

  int shape(unsigned int dim){
    ASSERT1("FixedDimArray_", "%ui is greater than the number of dimensions",
	    dim, dim<D);
    return (int)shape_[dim];
  }

  intp size(){
    intp out = 1;
    for (std::size_t i=0; i<D; i++){
      out*=shape_[i];
    }
    return out;
  }

  // Only arrays with the same numbers of dimensions can be swapped
  friend void swap(FixedDimArray_<T,D> &first, FixedDimArray_<T,D> &second){
    std::swap(first.dataMgr_, second.dataMgr_);
    std::swap(first.data_, second.data_);
    std::swap(first.offset_, second.offset_);
    std::swap(first.shape_, second.shape_);
    std::swap(first.stride_, second.stride_);
  }

protected: // methods to be reused by subclasses

  FixedDimArray_()
    : dataMgr_(),
      data_(NULL),
      offset_(0),
      shape_(),
      stride_()
  { };

  // Construct a numeric array that allocates its own data
  // args expects an integer for each dimension
  template<typename... Args, REQUIRE_INT(Args)>
  FixedDimArray_(Args... args);

  // Construct a numeric array that wraps an existing pointer
  // args expects an integer for each dimension
  template<typename... Args, REQUIRE_INT(Args)>
  FixedDimArray_(T* array, Args... args);
  
  void init_helper_(std::shared_ptr<dataWrapper<T>> &dataMgr, intp shape[D],
		    intp offset){
    data_ = dataMgr->get();
    dataMgr_ = dataMgr;
    offset_ = offset;

    std::size_t i = D;
    while (i>0){
      --i;
      shape_[i] = shape[i];
      if (i + 1 == D){
	stride_[i] = 1;
      } else {
	stride_[i] = shape_[i+1] * stride_[i+1];
      }
    }
  }

  void cleanup_helper_(){ data_ = NULL; }

protected: // attributes
  std::shared_ptr<dataWrapper<T>> dataMgr_; // manages ownership of data_
  // pointer to data (copied from dataMgr to provide faster access to elements)
  T* data_;
  intp offset_; // offset of the first element from the start of the pointer
  intp shape_[D]; // lists dimensions with increasing indexing speed
  intp stride_[D]; // stride_[D-1] is always 1
};


// Constructor of EnzoArray by allocating new data
template<typename T, std::size_t D>
template<typename... Args, class>
FixedDimArray_<T,D>::FixedDimArray_(Args... args)
{
  static_assert(D==sizeof...(args), "Incorrect number of dimensions");
  intp shape[D] = {((intp)args)...};
  intp size = 1;
  for (std::size_t i=0; i < D; i++){
    ASSERT("FixedDimArray_", "Positive dimensions are required.", shape[i]>0);
    size *= shape[i];
  }
  T* data = new T[size](); // allocate and set entries to 0
  std::shared_ptr<dataWrapper<T>> dataMgr;
  dataMgr = std::make_shared<dataWrapper<T>>(data,true);
  init_helper_(dataMgr, shape, 0);
}


// Constructor of array that wraps an existing c-style array
template<typename T, std::size_t D>
template<typename... Args, class>
FixedDimArray_<T,D>::FixedDimArray_(T* array, Args... args)
{
  static_assert(D==sizeof...(args), "Incorrect number of dimensions");
  intp shape[D] = {((intp)args)...};
  for (std::size_t i=0; i < D; i++){
    ASSERT("FixedDimArray_", "Positive dimensions are required.", shape[i]>0);
  }
  std::shared_ptr<dataWrapper<T>> dataMgr;
  dataMgr = std::make_shared<dataWrapper<T>>(array,false);
  init_helper_(dataMgr, shape, 0);
}


// Prepares an array of slices that refer to absolute start and stop values
// along each dimension. Also checks that the slices are valid
inline void prep_slices_(const ESlice* slices, const intp shape[],
			 const std::size_t D, ESlice* out_slices)
{
   for (std::size_t i=0; i<D; i++){
     intp start, stop;
     start = (slices[i].start<0) ? slices[i].start+shape[i] : slices[i].start;
     stop  = (slices[i].stop<0 ) ? slices[i].stop+shape[i]  : slices[i].stop;
     ASSERT3("FixedDimArray_",
	     "slice.start of %ld doesn't lie in bound of dim %ld of size %ld.",
	     (long)slices[i].start, (long)i, (long)shape[i], start < shape[i]);
     ASSERT3("FixedDimArray_",
	     "slice.stop of %d doesn't lie in bound of dim %ld of size %ld.",
	     (long)slices[i].stop, (long)i, (long)shape[i], stop <= shape[i]);
     ASSERT4("FixedDimArray_", ("slice.stop (%ld) doesn't exceed slice.start "
				"(%ld) must for dim %ld of size %ld."),
	     (long)slices[i].start, (long)slices[i].stop, (long)i,
	     (long)shape[i], stop>start);
     out_slices[i] = ESlice(start,stop);
  }
}

// Returnd TempArray_ representing a view of a subarray of the current instance
template<typename T, std::size_t D>
template<typename... Args, class>
TempArray_<T,D> FixedDimArray_<T,D>::subarray(Args... args){
  static_assert(D == sizeof...(args) || 0 == sizeof...(args),
		"Number of slices don't match number of dimensions");
  TempArray_<T,D> subarray;
  if (sizeof...(args) == 0) {
    subarray.init_helper_(dataMgr_,shape_,offset_);
  } else {
    ESlice in_slices[] = {args...};
    ESlice slices[D];
    prep_slices_(in_slices, shape_, D, slices);

    intp new_shape[D];
    intp new_offset = offset_;
    for (std::size_t dim=0; dim<D; dim++){
      new_shape[dim] = slices[dim].stop - slices[dim].start;
      new_offset += slices[dim].start * stride_[dim];
    }

    subarray.init_helper_(dataMgr_,new_shape,new_offset);
    for (std::size_t dim=0; dim<D; dim++){
      subarray.stride_[dim] = stride_[dim];
    }
  }
  return subarray;
}



template<typename T, std::size_t D>
TempArray_<T,D> FixedDimArray_<T,D>::deepcopy()
{
  T* data = new T[size()]; // allocate but don't initialize values
  std::shared_ptr<dataWrapper<T>> dataMgr;
  dataMgr = std::make_shared<dataWrapper<T>>(data,true);
  TempArray_<T,D> out;
  out.init_helper_(dataMgr, shape_, 0);
  out = *this;
  return out;
}



template<typename T, std::size_t D>
class TempArray_ : public FixedDimArray_<T,D>
{
  friend class FixedDimArray_<T,D>;
  friend class EnzoArray<T,D>;

public:
  // Assigns to each element of *this the value of val
  TempArray_<T,D>& operator=(const T& val);
  
  // Assigns to each element of *this the value of the corresponding element in
  // other. Sizes must match
  TempArray_<T,D>& operator=(const TempArray_<T,D>& other){
    assign_helper_(other.offset_,other.stride_,other.shape_, other.data_);
    return *this;
  }

  TempArray_<T,D>& operator=(const EnzoArray<T,D>& other){
    assign_helper_(other.offset_,other.stride_, other.shape_, other.data_);
    return *this;
  }

private:
  TempArray_() : FixedDimArray_<T,D>() { }

  void assign_helper_(const intp o_offset, const intp *o_stride,
		      const intp *o_shape, const T* o_data){
    for (std::size_t i = 0; i<D; i++){
      ASSERT("TempArray_","shapes aren't the same.",this->shape_[i]==o_shape[i]);
    }
    bool continue_outer_iter = true;
    intp indices[D] = {}; // all elements to 0
    while (continue_outer_iter){
      intp index = calc_index_(D, this->offset_, this->stride_, indices);
      intp o_index = calc_index_(D, o_offset, o_stride, indices);
      for (intp i = 0; i<this->shape_[D-1]; i++){
	this->data_[index] = o_data[o_index];
	index++; o_index++;
      }
      increment_outer_indices_(D, indices, this->shape_, continue_outer_iter);
    }
  }

};


template<typename T, std::size_t D>
TempArray_<T,D>& TempArray_<T,D>::operator=(const T& val)
{
  bool continue_outer_iter = true;
  intp indices[D] = {}; // all elements to 0
  while (continue_outer_iter){
    intp index = calc_index_(D, this->offset_, this->stride_, indices);
    for (intp i = 0; i<this->shape_[D-1]; i++){
      this->data_[index] = val;
      index++;
    }
    increment_outer_indices_(D, indices, this->shape_,continue_outer_iter);
  }
  return *this;
}


template<typename T, std::size_t D>
class EnzoArray : public FixedDimArray_<T,D>
{
public:
  // Default constructor. Constructs an unallocated EnzoArray.
  EnzoArray() : FixedDimArray_<T,D>() { }

  // Construct a numeric array that allocates its own data
  // args expects an integer for each dimension
  template<typename... Args, REQUIRE_INT(Args)>
  EnzoArray(Args... args) : FixedDimArray_<T,D>(args...) { }

  // Construct a numeric array that wraps an existing pointer
  // args expects an integer for each dimension
  template<typename... Args, REQUIRE_INT(Args)>
  EnzoArray(T* array, Args... args) : FixedDimArray_<T,D>(array, args...) { }

  // Copy constructor. Constructs a shallow copy of other.
  EnzoArray(const EnzoArray<T,D>& other){
    this->init_helper_(other.dataMgr_, other.shape_, other.offset_);
  }

  // Move constructor. Constructs the array with the contents of other
  EnzoArray(EnzoArray<T,D>&& other) : EnzoArray() {swap(*this,other);}
  EnzoArray(TempArray_<T,D>&& other) : EnzoArray() {swap(*this,other);}

  // Copy assignment operator. Makes *this a shallow copy of other. (Contents
  // of shallow copies and subarrays of *this are unaffected)
  EnzoArray<T,D>& operator=(const EnzoArray<T,D>& other){
    this->cleanup_helper_();
    init_helper_(other.dataMgr_, other.shape_, other.offset_);
    return *this;
  }

  // Move assignment operator. Replaces the contents of *this with those of
  // other. (Contents of shallow copies and subarrays of *this are unaffected)
  EnzoArray<T,D>& operator=(EnzoArray<T,D>&& other) {
    swap(*this,other);
    return *this;
  }
  EnzoArray<T,D>& operator=(TempArray_<T,D>&& other) {
    swap(*this,other);
    return *this;
  }
};

#endif /* ENZO_ENZO_ARRAY_HPP */
