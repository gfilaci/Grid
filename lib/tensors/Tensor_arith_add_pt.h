    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/tensors/Tensor_arith_add.h

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
#ifndef GRID_MATH_ARITH_ADD_PT_H
#define GRID_MATH_ARITH_ADD_PT_H

namespace Grid {

    
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////// ADD         ///////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    

// ADD is simple for now; cannot mix types and straightforward template
// Pert series +/- Pert series
  template<class vtype,class ltype,class rtype,int N> strong_inline void add(iPert<vtype,N> * __restrict__ ret,
								      const iPert<ltype,N> * __restrict__ lhs,
								      const iPert<rtype,N> * __restrict__ rhs)
  {
    for(int c=0;c<N;c++){
      ret->_internal[c]=lhs->_internal[c]+rhs->_internal[c];
    }
    return;
  }
  
  // + operator for perturbative series
  template<class ltype,class rtype,int N>
    strong_inline auto operator + (const iPert<ltype,N>& lhs,const iPert<rtype,N>& rhs) ->iPert<decltype(lhs._internal[0]+rhs._internal[0]),N>
    {
      typedef iPert<decltype(lhs._internal[0]+rhs._internal[0]),N> ret_t;
    ret_t ret;
    add(&ret,&lhs,&rhs);
    return ret;
    }

}

#endif
