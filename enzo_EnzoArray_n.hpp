#ifndef ENZO_ENZO_ARRAY_N_HPP
#define ENZO_ENZO_ARRAY_N_HPP

// Here we define the array with negative indexing


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


// Helper function to compute pointer index
// T is expected to be intp or int
template<typename T>
intp calc_index_(intp* stride, intp*shape, T first){
  if (first < 0){
    return *shape+first;
  } else {
    return first;
  }
}

template<typename T>
intp calc_index_(intp* stride, intp*shape, T first, T last){
  if (first < 0){
    if (last < 0){
      return (*stride)*(*shape + first) + last + shape[1];
    } else {
      return (*stride)*(*shape + first) + last;
    }
  } else if (last < 0){
    return (*stride)*first + last + shape[1];
  } else {
    return (*stride)*first + last;
  }
}

template<typename T, typename... Rest>
intp calc_index_(intp* stride, intp* shape, T first, Rest... rest){
  // gets unrolled at compile time
  if (first < 0){
    return (*stride)*(*shape + first) + calc_index_(++stride, ++shape, rest...);
  } else {
    return (*stride)*first + calc_index_(++stride, ++shape, rest...);
  }
}

template<typename T, std::size_t D>
class EnzoArrayN;

template<typename T, std::size_t D>
class TempArrayN_;

template<typename T, std::size_t D>
class FixedDimArrayN_
{
public:

  // Destructor
  ~FixedDimArrayN_() { cleanup_helper_();}

  // operator() - used to access array Elements
  template<typename... Args, REQUIRE_INT(Args)>
  T &operator() (Args... args) {
    static_assert(D==sizeof...(args),
		  "Number of indices don't match number of dimensions");
    return data_[offset_ + calc_index_(stride_,shape_,args...)];
  }
  template<typename... Args, REQUIRE_INT(Args)>
  T operator() (Args... args) const {
    static_assert(D==sizeof...(args),
		  "Number of indices don't match number of dimensions");
    return data_[offset_ + calc_index_(stride_,shape_,args...)];
  }

  // Specialized implementation for 3D arrays (reduces compile time)
  T &operator() (const int k, const int j, const int i){
    static_assert(D==3, "3 indices should only be specified for 3D arrays");
    if (k < 0){
      if (j < 0) {
	if (i < 0) {
	  return data_[offset_ +
		       (shape_[0]+k)*stride_[0] +
		       (shape_[1]+j)*stride_[1] +
		       (shape_[2]+i)];
	} else {
	  return data_[offset_ +
		       (shape_[0]+k)*stride_[0] +
		       (shape_[1]+j)*stride_[1] +
		       i];
	}
      } else {
	if (i < 0) {
	  return data_[offset_ +
		       (shape_[0]+k)*stride_[0] +
		       j*stride_[1] +
		       (shape_[2]+i)];
	} else {
	  return data_[offset_ +
		       (shape_[0]+k)*stride_[0] +
		       j*stride_[1] +
		       i];
	}
      }
    } else {
      if (j < 0) {
	if (i < 0) {
	  return data_[offset_ + k*stride_[0] +
		       (shape_[1]+j)*stride_[1] +
		       (shape_[2]+i)];
	} else {
	  return data_[offset_ + k*stride_[0] +
		       (shape_[1]+j)*stride_[1] +
		       i];
	}
      } else {
	if (i < 0) {
	  return data_[offset_ + k*stride_[0] +
		       j*stride_[1] +
		       (shape_[2]+i)];
	} else {
	  return data_[offset_ + k*stride_[0] + j*stride_[1] + i];
	}
      }
    } 
  }

  T operator() (const int k, const int j, const int i) const{
    static_assert(D==3, "3 indices should only be specified for 3D arrays");
    if (k < 0){
      if (j < 0) {
	if (i < 0) {
	  return data_[offset_ +
		       (shape_[0]+k)*stride_[0] +
		       (shape_[1]+j)*stride_[1] +
		       (shape_[2]+i)];
	} else {
	  return data_[offset_ +
		       (shape_[0]+k)*stride_[0] +
		       (shape_[1]+j)*stride_[1] +
		       i];
	}
      } else {
	if (i < 0) {
	  return data_[offset_ +
		       (shape_[0]+k)*stride_[0] +
		       j*stride_[1] +
		       (shape_[2]+i)];
	} else {
	  return data_[offset_ +
		       (shape_[0]+k)*stride_[0] +
		       j*stride_[1] +
		       i];
	}
      }
    } else {
      if (j < 0) {
	if (i < 0) {
	  return data_[offset_ + k*stride_[0] +
		       (shape_[1]+j)*stride_[1] +
		       (shape_[2]+i)];
	} else {
	  return data_[offset_ + k*stride_[0] +
		       (shape_[1]+j)*stride_[1] +
		       i];
	}
      } else {
	if (i < 0) {
	  return data_[offset_ + k*stride_[0] +
		       j*stride_[1] +
		       (shape_[2]+i)];
	} else {
	  return data_[offset_ + k*stride_[0] + j*stride_[1] + i];
	}
      }
    } 
  }


  // Produce a copy of the array.
  TempArrayN_<T,D> deepcopy();
  
  // Returns a subarray with same D. Expects an instance of ESlice for each
  // dimension. For a 3D EnzoArrayN, the function declaration might look like:
  // EnzoArrayN<T,3> subarray(ESlice k_slice, ESlice j_slice, ESlice i_slice);
  //
  // The Special Case of no arguments returns the full array
  template<typename... Args, REQUIRE_TYPE(Args,ESlice)>
  TempArrayN_<T,D> subarray(Args... args);

  int shape(unsigned int dim){
    ASSERT1("FixedDimArrayN_", "%ui is greater than the number of dimensions",
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
  friend void swap(FixedDimArrayN_<T,D> &first, FixedDimArrayN_<T,D> &second){
    std::swap(first.dataMgr_, second.dataMgr_);
    std::swap(first.data_, second.data_);
    std::swap(first.offset_, second.offset_);
    std::swap(first.shape_, second.shape_);
    std::swap(first.stride_, second.stride_);
  }

protected: // methods to be reused by subclasses

  FixedDimArrayN_()
    : dataMgr_(),
      data_(NULL),
      offset_(0),
      shape_(),
      stride_()
  { };

  // Construct a numeric array that allocates its own data
  // args expects an integer for each dimension
  template<typename... Args, REQUIRE_INT(Args)>
  FixedDimArrayN_(Args... args);

  // Construct a numeric array that wraps an existing pointer
  // args expects an integer for each dimension
  template<typename... Args, REQUIRE_INT(Args)>
  FixedDimArrayN_(T* array, Args... args);
  
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


// Constructor of EnzoArrayN by allocating new data
template<typename T, std::size_t D>
template<typename... Args, class>
FixedDimArrayN_<T,D>::FixedDimArrayN_(Args... args)
{
  static_assert(D==sizeof...(args), "Incorrect number of dimensions");
  intp shape[D] = {((intp)args)...};
  intp size = 1;
  for (std::size_t i=0; i < D; i++){
    ASSERT("FixedDimArrayN_", "Positive dimensions are required.", shape[i]>0);
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
FixedDimArrayN_<T,D>::FixedDimArrayN_(T* array, Args... args)
{
  static_assert(D==sizeof...(args), "Incorrect number of dimensions");
  intp shape[D] = {((intp)args)...};
  for (std::size_t i=0; i < D; i++){
    ASSERT("FixedDimArrayN_", "Positive dimensions are required.", shape[i]>0);
  }
  std::shared_ptr<dataWrapper<T>> dataMgr;
  dataMgr = std::make_shared<dataWrapper<T>>(array,false);
  init_helper_(dataMgr, shape, 0);
}

// Returnd TempArrayN_ representing a view of a subarray of the current instance
template<typename T, std::size_t D>
template<typename... Args, class>
TempArrayN_<T,D> FixedDimArrayN_<T,D>::subarray(Args... args){
  static_assert(D == sizeof...(args) || 0 == sizeof...(args),
		"Number of slices don't match number of dimensions");
  TempArrayN_<T,D> subarray;
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
TempArrayN_<T,D> FixedDimArrayN_<T,D>::deepcopy()
{
  T* data = new T[size()]; // allocate but don't initialize values
  std::shared_ptr<dataWrapper<T>> dataMgr;
  dataMgr = std::make_shared<dataWrapper<T>>(data,true);
  TempArrayN_<T,D> out;
  out.init_helper_(dataMgr, shape_, 0);
  out = *this;
  return out;
}



template<typename T, std::size_t D>
class TempArrayN_ : public FixedDimArrayN_<T,D>
{
  friend class FixedDimArrayN_<T,D>;
  friend class EnzoArrayN<T,D>;

public:
  // Assigns to each element of *this the value of val
  TempArrayN_<T,D>& operator=(const T& val);
  
  // Assigns to each element of *this the value of the corresponding element in
  // other. Sizes must match
  TempArrayN_<T,D>& operator=(const TempArrayN_<T,D>& other){
    assign_helper_(other.offset_,other.stride_,other.shape_, other.data_);
    return *this;
  }

  TempArrayN_<T,D>& operator=(const EnzoArrayN<T,D>& other){
    assign_helper_(other.offset_,other.stride_, other.shape_, other.data_);
    return *this;
  }

private:
  TempArrayN_() : FixedDimArrayN_<T,D>() { }

  void assign_helper_(const intp o_offset, const intp *o_stride,
		      const intp *o_shape, const T* o_data){
    for (std::size_t i = 0; i<D; i++){
      ASSERT("TempArrayN_","shapes aren't the same.",this->shape_[i]==o_shape[i]);
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
TempArrayN_<T,D>& TempArrayN_<T,D>::operator=(const T& val)
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
class EnzoArrayN : public FixedDimArrayN_<T,D>
{
public:
  // Default constructor. Constructs an unallocated EnzoArrayN.
  EnzoArrayN() : FixedDimArrayN_<T,D>() { }

  // Construct a numeric array that allocates its own data
  // args expects an integer for each dimension
  template<typename... Args, REQUIRE_INT(Args)>
  EnzoArrayN(Args... args) : FixedDimArrayN_<T,D>(args...) { }

  // Construct a numeric array that wraps an existing pointer
  // args expects an integer for each dimension
  template<typename... Args, REQUIRE_INT(Args)>
  EnzoArrayN(T* array, Args... args) : FixedDimArrayN_<T,D>(array, args...) { }

  // Copy constructor. Constructs a shallow copy of other.
  EnzoArrayN(const EnzoArrayN<T,D>& other){
    this->init_helper_(other.dataMgr_, other.shape_, other.offset_);
  }

  // Move constructor. Constructs the array with the contents of other
  EnzoArrayN(EnzoArrayN<T,D>&& other) : EnzoArrayN() {swap(*this,other);}
  EnzoArrayN(TempArrayN_<T,D>&& other) : EnzoArrayN() {swap(*this,other);}

  // Copy assignment operator. Makes *this a shallow copy of other. (Contents
  // of shallow copies and subarrays of *this are unaffected)
  EnzoArrayN<T,D>& operator=(const EnzoArrayN<T,D>& other){
    this->cleanup_helper_();
    init_helper_(other.dataMgr_, other.shape_, other.offset_);
    return *this;
  }

  // Move assignment operator. Replaces the contents of *this with those of
  // other. (Contents of shallow copies and subarrays of *this are unaffected)
  EnzoArrayN<T,D>& operator=(EnzoArrayN<T,D>&& other) {
    swap(*this,other);
    return *this;
  }
  EnzoArrayN<T,D>& operator=(TempArrayN_<T,D>&& other) {
    swap(*this,other);
    return *this;
  }
};

#endif /* ENZO_ENZO_ARRAY_N_HPP */
