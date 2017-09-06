/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/utils/TwistedFFT.h

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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
#ifndef QCD_UTIL_TWISTEDFFT_H
#define QCD_UTIL_TWISTEDFFT_H

#include "../action/gauge/twistmatrices_pt.h"

#define pPerpLoop(n1,n2) for(int n1=0; n1<Nc; n1++) for(int n2=0; n2<Nc; n2++)

namespace Grid {

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //  TwistMult is a multiplication.
  //  The only difference bewteen mult and TwistMult is that here a matrix is supposed to collect
  //  the fine momentum (pperp) degree of freedom.
  //  Therefore a matrix multiplication has to be done 'element by element' and not 'row by column'.
  //////////////////////////////////////////////////////////////////////////////////////////////////////

template<class rtype,class vtype,class mtype>
strong_inline void TwistMult(iScalar<rtype> * __restrict__ ret,
                       const iScalar<vtype> * __restrict__ lhs,
                       const iScalar<mtype> * __restrict__ rhs){
    TwistMult(&ret->_internal,&lhs->_internal,&rhs->_internal);
}
template<class rrtype,class ltype,class rtype,int N>
strong_inline void TwistMult(iMatrix<rrtype,N>* __restrict__ ret,
                       const iMatrix<ltype,N> * __restrict__ lhs,
                       const iMatrix<rtype,N> * __restrict__ rhs){
  for(int c1=0;c1<N;c1++){
    for(int c2=0;c2<N;c2++){
        mult(&ret->_internal[c1][c2],&lhs->_internal[c1][c2],&rhs->_internal[c1][c2]);
    }
  }
  return;
}
template<class rtype,class vtype,class mtype,int N>
strong_inline void TwistMult(iVector<rtype,N> * __restrict__ ret,
                       const iScalar<vtype>   * __restrict__ lhs,
                       const iVector<mtype,N> * __restrict__ rhs){
    for(int c1=0;c1<N;c1++){
        TwistMult(&ret->_internal[c1],&lhs->_internal,&rhs->_internal[c1]);
    }
}
template<class rtype,class vtype,class mtype,int N>
strong_inline void TwistMult(iVector<rtype,N> * __restrict__ ret,
                       const iVector<vtype,N> * __restrict__ lhs,
                       const iScalar<mtype>   * __restrict__ rhs){
    for(int c1=0;c1<N;c1++){
        TwistMult(&ret->_internal[c1],&lhs->_internal[c1],&rhs->_internal);
    }
}
template<class rtype,class vtype,class mtype,int N>
strong_inline void TwistMult(iPert<rtype,N> * __restrict__ ret,
                       const iScalar<vtype> * __restrict__ lhs,
                       const iPert<mtype,N> * __restrict__ rhs){
    for(int c1=0;c1<N;c1++){
        TwistMult(&ret->_internal[c1],&lhs->_internal,&rhs->_internal[c1]);
    }
}
template<class rtype,class vtype,class mtype,int N>
strong_inline void TwistMult(iPert<rtype,N> * __restrict__ ret,
                       const iPert<vtype,N> * __restrict__ lhs,
                       const iScalar<mtype> * __restrict__ rhs){
    for(int c1=0;c1<N;c1++){
        TwistMult(&ret->_internal[c1],&lhs->_internal[c1],&rhs->_internal);
    }
}
template<class obj1,class obj2,class obj3>
strong_inline void TwistMult(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const Lattice<obj3> &rhs){
    ret.checkerboard = lhs.checkerboard;
    conformable(ret,rhs);
    conformable(lhs,rhs);
    parallel_for(int ss=0;ss<lhs._grid->oSites();ss++){
#ifdef STREAMING_STORES
      obj1 tmp;
      TwistMult(&tmp,&lhs._odata[ss],&rhs._odata[ss]);
      vstream(ret._odata[ss],tmp);
#else
      TwistMult(&ret._odata[ss],&lhs._odata[ss],&rhs._odata[ss]);
#endif
    }
  }

namespace QCD {
namespace QCDpt {

template <class Gimpl>
class TwistedFFT {

    typedef typename std::remove_reference<decltype(Gimpl::twist.omega[0])>::type TwistBase;
    typedef typename Gimpl::ScalarField LatticeScalar;
    typedef typename Gimpl::MatrixField LatticeMatrix;
    
private:
    
    GridCartesian* grid;
    FFT theFFT;
    LatticeMatrix pPerpPhase;
    
public:
    
    iMatrix<TwistBase, Nc> Gamma, adjGamma;
    std::vector<Real> boundary_exp;
    int t1,t2;
    
    TwistedFFT(GridCartesian* grid_, std::vector<Complex> boundary_phases_):
    grid(grid_),
    theFFT(grid_),
    pPerpPhase(grid_)
    {
        // This TwistedFFT is tailored for twist on a plane
        // with generic orientation.
        int ntwisteddir = 0;
        for (int mu=0; mu<Nd; mu++) {
            if(istwisted(mu)) {
                ntwisteddir++;
                if (ntwisteddir==1) t1 = mu;
                if (ntwisteddir==2) t2 = mu;
            }
        }
        assert(ntwisteddir==2);
        
        // The FFTW library has FFTW_FORWARD = -1 by default.
        // Here the same convention is assumed.
        if(FFTW_FORWARD==+1) assert(0);
        
        // initialise boundary exponents
        for (int d=0; d<Nd; d++) {
            boundary_exp.push_back(log(boundary_phases_[d]).imag());
        }
        
        // initialise Fourier twist base
        pPerpLoop(n1,n2) {
            Gamma(n1,n2) = BuildGamma(n1,n2);
            adjGamma(n1,n2) = adj(Gamma(n1,n2));
            Gamma(n1,n2) = (1./(double)Nc)*Gamma(n1,n2);
        }
        
        // initialise phases for pPerp projection
        LatticeScalar xmu(grid), tmp(grid);
        pPerpLoop(n1,n2) {
            tmp = zero;
            for (int mu=0; mu<Nd; mu++) {
                Real pPerp = boundary_exp[mu] / (double)(grid->_fdimensions[mu]);
                if(mu==t1) pPerp +=  (double)n1 * 2. * M_PI / (double)(grid->_fdimensions[mu]) / (double)Nc;
                if(mu==t2) pPerp +=  (double)n2 * 2. * M_PI / (double)(grid->_fdimensions[mu]) / (double)Nc;
                LatticeCoordinate(xmu,mu);
                tmp = tmp + pPerp*xmu;
            }
            Complex ci(0.,1.);
            tmp = exp(ci*tmp);
            pokeColour(pPerpPhase,tmp,n1,n2);
        }
    }


TwistBase BuildGamma(int n1, int n2){
    TwistBase tmp(1.0);
    for (int i=0; i<n1; i++) tmp *= Gimpl::twist.omega[t2];
    tmp = adj(tmp);
    for (int i=0; i<n2; i++) tmp = Gimpl::twist.omega[t1]*tmp;
    return tmp;
}

void OrthonormalityTest(){
    bool testpassed = true;
    double precision = 1e-13;
    pPerpLoop(n1,n2){
        pPerpLoop(m1,m2){
            double tmp = abs(TensorRemove(trace(adjGamma(m1,m2)*Gamma(n1,n2))));
            if(m1==n1 && m2==n2) tmp -= 1.;
            if(std::abs(tmp)>precision) testpassed = false;
        }
    }
    
    if(testpassed) std::cout << GridLogMessage << "Orthonormality test SUCCEEDED" << std::endl;
    else std::cout << GridLogMessage << "Orthonormality test FAILED" << std::endl;
}

template<class vobj>
void pPerpProjectionForward(Lattice<vobj> &result, const Lattice<vobj> &source){
    
    decltype(peekColour(source,0,0)) tmp(grid);
    
    pPerpLoop(n1,n2){
        tmp = trace(adjGamma(n1,n2)*source);
        pokeColour(result,tmp,n1,n2);
    }
    TwistMult(result,conjugate(pPerpPhase),result);
}

template<class vobj>
void pPerpProjectionBackward(Lattice<vobj> &result,const Lattice<vobj> &source){
    
    Lattice<vobj> tmp(grid);
    TwistMult(tmp,pPerpPhase,source);
    
    result = zero;
    pPerpLoop(n1,n2){
        result = result + Gamma(n1,n2) * peekColour(tmp,n1,n2);
    }
}

template<class vobj>
void FFTforward(Lattice<vobj> &result,const Lattice<vobj> &source){
    conformable(result._grid,grid);
    conformable(source._grid,grid);
    
    Lattice<vobj> tmp(grid);
    pPerpProjectionForward(tmp,source);
    theFFT.FFT_all_dim(result,tmp,FFT::forward);
}

template<class vobj>
void FFTbackward(Lattice<vobj> &result,const Lattice<vobj> &source){
    conformable(result._grid,grid);
    conformable(source._grid,grid);
    
    Lattice<vobj> tmp(grid);
    theFFT.FFT_all_dim(tmp,source,FFT::backward);
    pPerpProjectionBackward(result,tmp);
}


};

}
}
}
#endif
