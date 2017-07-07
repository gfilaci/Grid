    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/tensors/Tensor_arith_mac_pt.h

    Copyright (C) 2015

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

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#ifndef GRID_MATH_ARITH_MAC_PT_H
#define GRID_MATH_ARITH_MAC_PT_H

namespace Grid {


///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////// MAC         ///////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////

    ///////////////////////////
    // Legal multiplication table
    ///////////////////////////
    // scal x pert = pert
    // pert x scal = pert
    // pert x pert = pert
    ///////////////////////////
template<class rrtype,class ltype,class rtype,int N>
strong_inline void mac(iPert<rrtype,N> * __restrict__ ret,const iScalar<ltype> * __restrict__ lhs,const iPert<rtype,N> * __restrict__ rhs)
{
    for(int c1=0;c1<N;c1++){
        mac(&ret->_internal[c1],&lhs->_internal,&rhs->_internal[c1]);
    }
    return;
}
template<class rrtype,class ltype,class rtype,int N>
strong_inline void mac(iPert<rrtype,N> * __restrict__ ret,const iPert<ltype,N> * __restrict__ lhs,const iScalar<rtype> * __restrict__ rhs)
{
    for(int c1=0;c1<N;c1++){
        mac(&ret->_internal[c1],&lhs->_internal[c1],&rhs->_internal);
    }
    return;
}
template<class rrtype,class ltype,class rtype,int N>
strong_inline void mac(iPert<rrtype,N> * __restrict__ ret,const iPert<ltype,N> * __restrict__ lhs,const iPert<rtype,N> * __restrict__ rhs)
{
    for(int c1=0;c1<N;c1++){
    for(int c2=0;c2<=c1;c2++){
        mac(&ret->_internal[c1],&lhs->_internal[c2],&rhs->_internal[c1-c2]);
    }}
    return;
}
}

#endif
