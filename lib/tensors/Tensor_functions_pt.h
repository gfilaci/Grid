    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/tensors/Tensor_functions_pt.h

    Copyright (C) 2015

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
#ifndef GRID_FUNCTIONS_PT_H
#define GRID_FUNCTIONS_PT_H

namespace Grid {

  /////////////////////////////////////////////// 
  // Exponentiate perturbative series
  ///////////////////////////////////////////////
   
template<class vtype> inline iScalar<vtype> Exponentiate(const iScalar<vtype>&r)
    {
      iScalar<vtype> ret;
      ret._internal = Exponentiate(r._internal);
      return ret;
    }

template<class vtype, int N> inline iVector<vtype, N> Exponentiate(const iVector<vtype,N>&r)
    {
      iVector<vtype, N> ret;
      for (int i = 0; i < N; i++)
        ret._internal[i] = Exponentiate(r._internal[i]);
      return ret;
    }

template<class vtype, int N> inline iMatrix<vtype, N> Exponentiate(const iMatrix<vtype,N>&r)
    {
      iMatrix<vtype, N> ret;
      for (int i = 0; i < N; i++)
        ret._internal[i] = Exponentiate(r._internal[i]);
      return ret;
    }
    
template<class vtype, int N> inline iPert<vtype, N> Exponentiate(const iPert<vtype,N> &P)
  {
      // ASSUMING P(0) = 0
      iPert<vtype, N> ret(P), newtmp(P), tmp;
      
      typedef vtype mytype;
      mytype unit(1.0);
      
      ret._internal[0] = unit;
      
      for(int k=2; k<N; k++){
          zeroit(tmp);
          // i runs from 1 to N-1 (interval where P is defined)
          // but newtemp starts at order k-1
          // so there is a contribution only for i<N+1-k
          for(int i=1; i<N+1-k; i++){
              // j runs from k-1 to N-1 (interval where newtmp is defined)
              // but imposing i+j<N leads to j<N-i
              for(int j=k-1; j<N-i; j++){
                  // now k<=i+k<N
                  tmp._internal[i+j] += newtmp._internal[j] * P._internal[i];
              }
          }
          
          newtmp = (1./(double)k) * tmp;
          ret += newtmp;
      }
      
      return ret;
  }
  
  ///////////////////////////////////////////////
  // Logarithm of perturbative series
  ///////////////////////////////////////////////
  
  template<class vtype> inline iScalar<vtype> Logarithm(const iScalar<vtype>&r)
    {
      iScalar<vtype> ret;
      ret._internal = Logarithm(r._internal);
      return ret;
    }

template<class vtype, int N> inline iVector<vtype, N> Logarithm(const iVector<vtype,N>&r)
    {
      iVector<vtype, N> ret;
      for (int i = 0; i < N; i++)
        ret._internal[i] = Logarithm(r._internal[i]);
      return ret;
    }

template<class vtype, int N> inline iMatrix<vtype, N> Logarithm(const iMatrix<vtype,N>&r)
    {
      iMatrix<vtype, N> ret;
      for (int i = 0; i < N; i++)
        ret._internal[i] = Logarithm(r._internal[i]);
      return ret;
    }
    
template<class vtype, int N> inline iPert<vtype, N> Logarithm(const iPert<vtype,N> &P)
  {
      // ASSUMING P(0) = 1
      iPert<vtype, N> ret(P), newtmp(P), tmp;
      double factor, sign = 1.;
      
      zeroit(ret._internal[0]);
      
      for(int k=2; k<N; k++){
          zeroit(tmp);
          // i runs from 1 to N-1 (interval where P is defined)
          // but newtemp starts at order k-1
          // so there is a contribution only for i<N+1-k
          for(int i=1; i<N+1-k; i++){
              // j runs from k-1 to N-1 (interval where newtmp is defined)
              // but imposing i+j<N leads to j<N-i
              for(int j=k-1; j<N-i; j++){
                  // now k<=i+k<N
                  tmp._internal[i+j] += newtmp._internal[j] * P._internal[i];
              }
          }
          
          newtmp = tmp;
          sign = -sign;
          factor = sign/(double)k;
          ret += factor * newtmp;
      }
      
      return ret;
  }

}
#endif
