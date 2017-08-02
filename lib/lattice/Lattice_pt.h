    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/lattice/Lattice_trace.h

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
#ifndef GRID_LATTICE_PT_H
#define GRID_LATTICE_PT_H

namespace Grid {
    
template<class vobj>
    inline auto Exponentiate(const Lattice<vobj> &lhs)
    -> Lattice<decltype(Exponentiate(lhs._odata[0]))>
    {
        Lattice<decltype(Exponentiate(lhs._odata[0]))> ret(lhs._grid);
        parallel_for(int ss=0;ss<lhs._grid->oSites();ss++){
            ret._odata[ss] = Exponentiate(lhs._odata[ss]);
        }
        return ret;
    }

template<class vobj>
    inline auto Logarithm(const Lattice<vobj> &lhs)
    -> Lattice<decltype(Logarithm(lhs._odata[0]))>
    {
        Lattice<decltype(Logarithm(lhs._odata[0]))> ret(lhs._grid);
        parallel_for(int ss=0;ss<lhs._grid->oSites();ss++){
            ret._odata[ss] = Logarithm(lhs._odata[ss]);
        }
        return ret;
    }

template<class vobj1,class vobj2>
    inline auto AddToOrd(const int &ord, const Lattice<vobj1> &lhs, const Lattice<vobj2> &rhs)
    -> Lattice<decltype(AddToOrd(ord,lhs._odata[0],rhs._odata[0]))>
    {
        Lattice<decltype(AddToOrd(ord,lhs._odata[0],rhs._odata[0]))> ret(lhs._grid);
        parallel_for(int ss=0;ss<lhs._grid->oSites();ss++){
            ret._odata[ss] = AddToOrd(ord,lhs._odata[ss],rhs._odata[ss]);
        }
        return ret;
    }

template<class vobj1,class vobj2>
    inline void AddToOrdVoid(const int &ord, Lattice<vobj1> &lhs, const Lattice<vobj2> &rhs)
    {
        parallel_for(int ss=0;ss<lhs._grid->oSites();ss++){
            AddToOrdVoid(ord,lhs._odata[ss],rhs._odata[ss]);
        }
    }

template<class vobj1,class vobj2>
    inline void AddToOrdVoid(const int &ord, Lattice<vobj1> &lhs, const Lattice<vobj2> &rhs, const RealD factor)
    {
        parallel_for(int ss=0;ss<lhs._grid->oSites();ss++){
            AddToOrdVoid(ord,lhs._odata[ss],rhs._odata[ss],factor);
        }
    }
    
template<class vobj1,class vobj2>
    inline auto ShiftedSum(const int &ord, const Lattice<vobj1> &lhs, const Lattice<vobj2> &rhs)
    -> Lattice<decltype(ShiftedSum(ord,lhs._odata[0],rhs._odata[0]))>
    {
        Lattice<decltype(ShiftedSum(ord,lhs._odata[0],rhs._odata[0]))> ret(lhs._grid);
        parallel_for(int ss=0;ss<lhs._grid->oSites();ss++){
            ret._odata[ss] = ShiftedSum(ord,lhs._odata[ss],rhs._odata[ss]);
        }
        return ret;
    }
  
template<class vobj1,class vobj2>
    inline void ShiftedSumVoid(const int &ord, Lattice<vobj1> &lhs, const Lattice<vobj2> &rhs, const RealD factor)
    {
        parallel_for(int ss=0;ss<lhs._grid->oSites();ss++){
            ShiftedSumVoid(ord,lhs._odata[ss],rhs._odata[ss],factor);
        }
    }
    
template<class vobj>
    inline auto Pnorm2_internal(const Lattice<vobj> &lhs)
    -> Lattice<decltype(Pnorm2(lhs._odata[0]))>
    {
        Lattice<decltype(Pnorm2(lhs._odata[0]))> ret(lhs._grid);
        parallel_for(int ss=0;ss<lhs._grid->oSites();ss++){
            ret._odata[ss] = Pnorm2(lhs._odata[ss]);
        }
        return ret;
    }
    
template<class vobj>
    inline auto Pnorm2(const Lattice<vobj> &lhs)
    -> decltype(sum(Pnorm2_internal(lhs)))
    {
        return sum(Pnorm2_internal(lhs)) / (double)lhs._grid->gSites();
    }
    
    // specialisation needed for the critical mass
template<class sobj,class vobj, int N, int M> strong_inline
  void axpy(Lattice<iScalar<iVector<iPert<vobj,M>,N>>> &ret,iPert<sobj,M> a,const Lattice<iScalar<iVector<iPert<vobj,M>,N>>> &x,const Lattice<iScalar<iVector<iPert<vobj,M>,N>>> &y){
    ret.checkerboard = x.checkerboard;
    conformable(ret,x);
    conformable(x,y);
    parallel_for(int ss=0;ss<x._grid->oSites();ss++){
#ifdef STREAMING_STORES
        iScalar<iVector<iPert<vobj,M>,N>> tmp = zero;
        for (int alpha=0; alpha<N; alpha++) {
            for(int c1=0;c1<M;c1++){
                for(int c2=0;c2<=c1;c2++){
                    tmp._internal._internal[alpha]._internal[c1]+=a._internal[c2]*(x._odata[ss]._internal._internal[alpha]._internal[c1-c2]);
                }}
            tmp._internal._internal[alpha]+=y._odata[ss]._internal._internal[alpha];
        }
        vstream(ret._odata[ss],tmp);
#else
        ret._odata[ss] = zero;
        for (int alpha=0; alpha<N; alpha++) {
            for(int c1=0;c1<M;c1++){
                for(int c2=0;c2<=c1;c2++){
                    ret._odata[ss]._internal._internal[alpha]._internal[c1]+=a._internal[c2]*(x._odata[ss]._internal._internal[alpha]._internal[c1-c2]);
                }}
            ret._odata[ss]._internal._internal[alpha]+=y._odata[ss]._internal._internal[alpha];
        }
#endif
    }
  }
    
}
#endif

