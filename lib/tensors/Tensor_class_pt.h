/*************************************************************************************
Grid physics library, www.github.com/paboyle/Grid
Source file: ./lib/tensors/Tensor_class_pt.h
Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Gianluca Filaci <g.filaci@ed.ac.uk>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
#ifndef GRID_MATH_TENSORS_PT_H
#define GRID_MATH_TENSORS_PT_H

namespace Grid {

///////////////////////////////////////////////////
// Perturbative series object.
///////////////////////////////////////////////////

template <class vtype, int N>
class iPert {
 public:
  vtype _internal[N];

  typedef vtype element;
  typedef typename GridTypeMapper<vtype>::scalar_type scalar_type;
  typedef typename GridTypeMapper<vtype>::vector_type vector_type;
  typedef typename GridTypeMapper<vtype>::vector_typeD vector_typeD;
  typedef typename GridTypeMapper<vtype>::tensor_reduced tensor_reduced_v;
  typedef typename GridTypeMapper<vtype>::scalar_object recurse_scalar_object;
  typedef iScalar<tensor_reduced_v> tensor_reduced;
  typedef iPert<recurse_scalar_object, N> scalar_object;

  // substitutes a real or complex version with same tensor structure
  typedef iPert<typename GridTypeMapper<vtype>::Complexified, N> Complexified;
  typedef iPert<typename GridTypeMapper<vtype>::Realified, N> Realified;

  // get double precision version
  typedef iPert<typename GridTypeMapper<vtype>::DoublePrecision, N> DoublePrecision;
  
  template <class T, typename std::enable_if<!isGridTensor<T>::value, T>::type
                         * = nullptr>
  strong_inline auto operator=(T arg) -> iPert<vtype, N> {
    zeroit(*this);
    for (int i = 0; i < N; i++) _internal[i] = arg;
    return *this;
  }

  enum { TensorLevel = GridTypeMapper<vtype>::TensorLevel + 1 };
  iPert(const Zero &z) { *this = zero; };
  iPert() = default;
  /*
  iPert(const iPert<vtype,N> &copyme)=default;
  iPert(iPert<vtype,N> &&copyme)=default;
  iPert<vtype,N> & operator= (const iPert<vtype,N> &copyme) = default;
  iPert<vtype,N> & operator= (iPert<vtype,N> &&copyme) = default;
  */

  iPert<vtype, N> &operator=(const Zero &hero) {
    zeroit(*this);
    return *this;
  }
  friend strong_inline void zeroit(iPert<vtype, N> &that) {
    for (int i = 0; i < N; i++) {
      zeroit(that._internal[i]);
    }
  }
  friend strong_inline void prefetch(iPert<vtype, N> &that) {
    for (int i = 0; i < N; i++) prefetch(that._internal[i]);
  }
  friend strong_inline void vstream(iPert<vtype, N> &out,
                                    const iPert<vtype, N> &in) {
    for (int i = 0; i < N; i++) {
      vstream(out._internal[i], in._internal[i]);
    }
  }
  friend strong_inline void vbroadcast(iPert<vtype,N> &out,const iPert<vtype,N> &in,int lane){
    for(int i=0;i<N;i++){
      vbroadcast(out._internal[i],in._internal[i],lane);
    }
  }
  friend strong_inline void permute(iPert<vtype,N> &out,const iPert<vtype,N> &in,int permutetype){
    for(int i=0;i<N;i++){
      permute(out._internal[i],in._internal[i],permutetype);
    }
  }
  friend strong_inline void rotate(iPert<vtype,N> &out,const iPert<vtype,N> &in,int rot){
    for(int i=0;i<N;i++){
      rotate(out._internal[i],in._internal[i],rot);
    }
  }
  friend strong_inline void exchange(iPert<vtype,N> &out1,iPert<vtype,N> &out2,
				     const iPert<vtype,N> &in1,const iPert<vtype,N> &in2,int type){
    for(int i=0;i<N;i++){
      exchange(out1._internal[i],out2._internal[i],
	        in1._internal[i], in2._internal[i],type);
    }
  }

  // Unary negation
  friend strong_inline iPert<vtype, N> operator-(const iPert<vtype, N> &r) {
    iPert<vtype, N> ret;
    for (int i = 0; i < N; i++) ret._internal[i] = -r._internal[i];
    return ret;
  }
  // *=,+=,-= operators inherit from corresponding "*,-,+" behaviour
  strong_inline iPert<vtype, N> &operator*=(const iScalar<vtype> &r) {
    *this = (*this) * r;
    return *this;
  }
  strong_inline iPert<vtype, N> &operator-=(const iPert<vtype, N> &r) {
    *this = (*this) - r;
    return *this;
  }
  strong_inline iPert<vtype, N> &operator+=(const iPert<vtype, N> &r) {
    *this = (*this) + r;
    return *this;
  }
  strong_inline vtype &operator()(int i) { return _internal[i]; }
  strong_inline const vtype &operator()(int i) const { return _internal[i]; }
  friend std::ostream &operator<<(std::ostream &stream,
                                  const iPert<vtype, N> &o) {
    stream << "P<" << N << ">{";
    for (int i = 0; i < N; i++) {
      stream << o._internal[i];
      if (i < N - 1) stream << ",";
    }
    stream << "}";
    return stream;
  };
  //    strong_inline vtype && operator ()(int i) {
  //      return _internal[i];
  //    }
};


template <class v, int N>
void vprefetch(const iPert<v, N> &vv) {
  for (int i = 0; i < N; i++) {
    vprefetch(vv._internal[i]);
  }
}

}
#endif
