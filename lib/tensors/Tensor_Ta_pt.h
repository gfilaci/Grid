    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/tensors/Tensor_Ta_pt.h

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
#ifndef GRID_MATH_TA_PT_H
#define GRID_MATH_TA_PT_H


namespace Grid {

  /////////////////////////////////////////////// 
  // Ta function for perturbative series
  ///////////////////////////////////////////////

  template<class vtype,int N> inline iPert<vtype,N> Ta(const iPert<vtype,N>&r)
    {
      iPert<vtype,N> ret;
      for(int i=0;i<N;i++){
        ret._internal[i] = Ta(r._internal[i]);
      }
      return ret;
    }

  /////////////////////////////////////////////// 
  // ProjectOnGroup function for perturbative series
  /////////////////////////////////////////////// 

/************************************************/
/*           QUAD PRECISION FUNCTIONS           */
/************************************************/
#ifdef USE_QUADPREC
////////////
// overloading POG in quad precision
///////////
  template<int N,int M> inline iPert<iMatrix<vComplexD,M>,N> ProjectOnGroup(const iPert<iMatrix<vComplexD,M>,N>&P)
{
  typedef iMatrix<   vComplexD,M> vtype;
  typedef iMatrix<    ComplexD,M> stype;
  typedef iMatrix<__complex128,M> qtype;
  
  qtype unit(1.0);
  
  // extract simd vector
  int Nsimd = sizeof(vtype::vector_type) / sizeof(vtype::scalar_type);
  std::vector<iPert<stype,N>> dbuf(Nsimd);
  std::vector<iPert<qtype,N>> qbuf(Nsimd);
  extract(P,dbuf);
  
  // cast to quad precision
  for(int i=0; i<Nsimd; i++){
    for(int j=0; j<N; j++){
      for(int l=0; l<M; l++){
	for(int m=0; m<M; m++){
	  __real__ qbuf[i](j)(l,m) = dbuf[i](j)(l,m).real();
	  __imag__ qbuf[i](j)(l,m) = dbuf[i](j)(l,m).imag();
	}
      }
    }
  }
  
  // POG
  // force the expansion to start from the identity...
  for(int z=0; z<Nsimd; z++){
    qbuf[z]._internal[0] = unit;
    qbuf[z] = QuadExponentiate(Ta(QuadLogarithm(qbuf[z])));
  }
  
  // cast to double precision
  for(int i=0; i<Nsimd; i++){
    for(int j=0; j<N; j++){
      for(int l=0; l<M; l++){
	for(int m=0; m<M; m++){
	  dbuf[i](j)(l,m) = ComplexD(__real__ qbuf[i](j)(l,m), __imag__ qbuf[i](j)(l,m));
	}
      }
    }
  }
  
  // merge simd vector
  iPert<vtype, N> result;
  merge(result,dbuf);
  
  return result;
}
#endif

template<class vtype,int N> inline iPert<vtype,N> ProjectOnGroup(const iPert<vtype,N>&r)
    {
      iPert<vtype,N> ret(r);

      // force the expansion to start from the identity...
      typedef vtype mytype;
      mytype unit(1.0);
      ret._internal[0] = unit;

      return Exponentiate(Ta(Logarithm(ret)));
    }

}

#endif
