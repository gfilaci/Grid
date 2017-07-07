    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/tensors/Tensor_reality_pt.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: neo <cossu@post.kek.jp>
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

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#ifndef GRID_MATH_REALITY_PT_H
#define GRID_MATH_REALITY_PT_H
namespace Grid {

/////////////////////////////////////////////// 
// multiply by I; make recursive.
/////////////////////////////////////////////// 

template<class vtype,int N> inline iPert<vtype,N> timesI(const iPert<vtype,N>&r)
{
  iPert<vtype,N> ret;
  for(int i=0;i<N;i++){
    timesI(ret._internal[i],r._internal[i]);
  }
  return ret;
}
template<class vtype,int N> inline void timesI(iPert<vtype,N> &ret,const iPert<vtype,N>&r)
{
  for(int i=0;i<N;i++){
    timesI(ret._internal[i],r._internal[i]);
  }
}
template<class vtype,int N> inline iPert<vtype,N> timesMinusI(const iPert<vtype,N>&r)
{
  iPert<vtype,N> ret;
  for(int i=0;i<N;i++){
    timesMinusI(ret._internal[i],r._internal[i]);
  }
  return ret;
}
template<class vtype,int N> inline void timesMinusI(iPert<vtype,N> &ret,const iPert<vtype,N>&r)
{
  for(int i=0;i<N;i++){
    timesMinusI(ret._internal[i],r._internal[i]);
  }
}


/////////////////////////////////////////////// 
// Conj function for perturbative series
/////////////////////////////////////////////// 

template<class vtype,int N> inline iPert<vtype,N> conjugate(const iPert<vtype,N>&r)
{
  iPert<vtype,N> ret;
  for(int i=0;i<N;i++){
    ret._internal[i] = conjugate(r._internal[i]);
  }
  return ret;
}


/////////////////////////////////////////////// 
// Adj function for perturbative series
/////////////////////////////////////////////// 

template<class vtype,int N> inline iPert<vtype,N> adj(const iPert<vtype,N>&r)
{
    iPert<vtype,N> ret;
    for(int i=0;i<N;i++){
        ret._internal[i] = adj(r._internal[i]);
    }
    return ret;
}


/////////////////////////////////////////////////////////////////
// Can only take the real/imag part of scalar objects, since
// lattice objects of different complex nature are non-conformable.
/////////////////////////////////////////////////////////////////

template<class itype,int N> inline auto real(const iPert<itype,N> &z) -> iPert<decltype(real(z._internal[0])),N>
{
    iPert<decltype(real(z._internal[0])),N> ret;
    for(int c1=0;c1<N;c1++){
        ret._internal[c1] = real(z._internal[c1]);
    }
    return ret;
}
template<class itype,int N> inline auto imag(const iPert<itype,N> &z) -> iPert<decltype(imag(z._internal[0])),N>
{
    iPert<decltype(imag(z._internal[0])),N> ret;
    for(int c1=0;c1<N;c1++){
        ret._internal[c1] = imag(z._internal[c1]);
    }
    return ret;
}

}
#endif
