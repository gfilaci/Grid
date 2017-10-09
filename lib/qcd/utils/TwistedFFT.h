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

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //  TwistDivide is a division.
  //  Here a matrix is supposed to collect the fine momentum (pperp) degree of freedom.
  //  Therefore a division has to be done 'element by element'.
  //////////////////////////////////////////////////////////////////////////////////////////////////////

template<class rtype,class vtype,class mtype>
strong_inline void TwistDivide(iScalar<rtype> * __restrict__ ret,
                       const iScalar<vtype> * __restrict__ lhs,
                       const iScalar<mtype> * __restrict__ rhs){
    TwistDivide(&ret->_internal,&lhs->_internal,&rhs->_internal);
}
template<class rrtype,class ltype,class rtype,int N>
strong_inline void TwistDivide(iMatrix<rrtype,N>* __restrict__ ret,
                       const iMatrix<ltype,N> * __restrict__ lhs,
                       const iMatrix<rtype,N> * __restrict__ rhs){
  for(int c1=0;c1<N;c1++){
    for(int c2=0;c2<N;c2++){
        ret->_internal[c1][c2] = lhs->_internal[c1][c2] / rhs->_internal[c1][c2];
    }
  }
  return;
}
template<class rtype,class vtype,class mtype,int N>
strong_inline void TwistDivide(iVector<rtype,N> * __restrict__ ret,
                       const iScalar<vtype>   * __restrict__ lhs,
                       const iVector<mtype,N> * __restrict__ rhs){
    for(int c1=0;c1<N;c1++){
        TwistDivide(&ret->_internal[c1],&lhs->_internal,&rhs->_internal[c1]);
    }
}
template<class rtype,class vtype,class mtype,int N>
strong_inline void TwistDivide(iVector<rtype,N> * __restrict__ ret,
                       const iVector<vtype,N> * __restrict__ lhs,
                       const iScalar<mtype>   * __restrict__ rhs){
    for(int c1=0;c1<N;c1++){
        TwistDivide(&ret->_internal[c1],&lhs->_internal[c1],&rhs->_internal);
    }
}
template<class rtype,class vtype,class mtype,int N>
strong_inline void TwistDivide(iPert<rtype,N> * __restrict__ ret,
                       const iScalar<vtype> * __restrict__ lhs,
                       const iPert<mtype,N> * __restrict__ rhs){
    for(int c1=0;c1<N;c1++){
        TwistDivide(&ret->_internal[c1],&lhs->_internal,&rhs->_internal[c1]);
    }
}
template<class rtype,class vtype,class mtype,int N>
strong_inline void TwistDivide(iPert<rtype,N> * __restrict__ ret,
                       const iPert<vtype,N> * __restrict__ lhs,
                       const iScalar<mtype> * __restrict__ rhs){
    for(int c1=0;c1<N;c1++){
        TwistDivide(&ret->_internal[c1],&lhs->_internal[c1],&rhs->_internal);
    }
}
template<class obj1,class obj2,class obj3>
strong_inline void TwistDivide(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const Lattice<obj3> &rhs){
    ret.checkerboard = lhs.checkerboard;
    conformable(ret,rhs);
    conformable(lhs,rhs);
    parallel_for(int ss=0;ss<lhs._grid->oSites();ss++){
#ifdef STREAMING_STORES
      obj1 tmp;
      TwistDivide(&tmp,&lhs._odata[ss],&rhs._odata[ss]);
      vstream(ret._odata[ss],tmp);
#else
      TwistDivide(&ret._odata[ss],&lhs._odata[ss],&rhs._odata[ss]);
#endif
    }
}

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //  TWISTED FAST FOURIER TRANSFORM
  //////////////////////////////////////////////////////////////////////////////////////////////////////
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
    LatticeMatrix pPerpPhase;            // exp(i (p_perp + theta/L ) * x)
    LatticeMatrix halfphatsq;            // 2 * sum_mu sin^2(p_mu/2) / den
    std::vector<LatticeMatrix> pbarmu;   // sin(p_mu) / den
    // den = halfphatsq^2 + sum_mu pbarmu^2
    
public:
    
    iMatrix<TwistBase, Nc> Gamma, adjGamma;
    std::vector<Real> boundary_exp;
    std::vector<Complex> boundary_phases;
    int t1,t2;
    
    TwistedFFT(GridCartesian* grid_, std::vector<Complex> boundary_phases_):
    grid(grid_),
    theFFT(grid_),
    pPerpPhase(grid_),
    halfphatsq(grid_),
    boundary_phases(boundary_phases_)
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
        
        // INITIALISE Fourier twist base
        pPerpLoop(n1,n2) {
            Gamma(n1,n2) = BuildGamma(n1,n2);
            adjGamma(n1,n2) = adj(Gamma(n1,n2));
            Gamma(n1,n2) = (1./(double)Nc)*Gamma(n1,n2);
        }
        
        for (int mu=0; mu<Nd; mu++) {
            pbarmu.push_back(LatticeMatrix(grid_));
            boundary_exp.push_back(0.);
        }
        
        FFTinitialisation(boundary_phases_);
    }
    
    
    void FFTinitialisation(std::vector<Complex> boundary_phases_){
        
        // INITIALISE boundary exponents
        for (int d=0; d<Nd; d++) {
            boundary_exp[d] = log(boundary_phases_[d]).imag();
        }
        
        // INITIALISE momenta for Wilson propagator
        LatticeScalar xmu(grid), tmp(grid);
        // --- pbar ---
        pbarmu.reserve(Nd);
        for (int mu=0; mu<Nd; mu++) {
            LatticeCoordinate(xmu,mu);
            xmu *= 2. * M_PI / (double)(grid->_fdimensions[mu]);
            xmu += boundary_exp[mu] / (double)(grid->_fdimensions[mu]);
            pPerpLoop(n1,n2) {
                tmp = xmu;
                if(mu==t1) tmp +=  (double)n1 * 2. * M_PI / (double)(grid->_fdimensions[mu]) / (double)Nc;
                if(mu==t2) tmp +=  (double)n2 * 2. * M_PI / (double)(grid->_fdimensions[mu]) / (double)Nc;
                pokeColour(pbarmu[mu],sin(tmp),n1,n2);
            }
        }
        // --- khatsq ---
        pPerpLoop(n1,n2) {
            tmp = zero;
            for (int mu=0; mu<Nd; mu++) {
                LatticeCoordinate(xmu,mu);
                xmu *= M_PI / (double)(grid->_fdimensions[mu]);
                xmu += 0.5 * boundary_exp[mu] / (double)(grid->_fdimensions[mu]);
                if(mu==t1) xmu +=  (double)n1 * M_PI / (double)(grid->_fdimensions[mu]) / (double)Nc;
                if(mu==t2) xmu +=  (double)n2 * M_PI / (double)(grid->_fdimensions[mu]) / (double)Nc;
                tmp = tmp + 2.*sin(xmu)*sin(xmu);
            }
            pokeColour(halfphatsq,tmp,n1,n2);
        }
        
        // compute denominator, use pPerpPhase for its temporary storage
        LatticeMatrix tmpMatrix(grid);
        pPerpPhase = zero;
        for (int mu=0; mu<Nd; mu++) {
            TwistMult(tmpMatrix,pbarmu[mu],pbarmu[mu]);
            pPerpPhase += tmpMatrix;
        }
        TwistMult(tmpMatrix,halfphatsq,halfphatsq);
        pPerpPhase += tmpMatrix;
        
        TwistDivide(halfphatsq,halfphatsq,pPerpPhase);
        for (int mu=0; mu<Nd; mu++){
            TwistDivide(pbarmu[mu],pbarmu[mu],pPerpPhase);
        }
        
        
        // INITIALISE phases for pPerp projection
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
    
    Lattice<vobj> tmpresult(grid);
    decltype(peekColour(source,0,0)) tmp(grid);
    
    pPerpLoop(n1,n2){
        tmp = trace(adjGamma(n1,n2)*source);
        pokeColour(tmpresult,tmp,n1,n2);
    }
    
    TwistMult(result,conjugate(pPerpPhase),tmpresult);
}

template<class vobj>
void pPerpProjectionBackward(Lattice<vobj> &result,const Lattice<vobj> &source){
    
    Lattice<vobj> tmp(grid), tmpresult(grid);
    TwistMult(tmp,pPerpPhase,source);
    
    tmpresult = zero;
    pPerpLoop(n1,n2){
        tmpresult = tmpresult + Gamma(n1,n2) * peekColour(tmp,n1,n2);
    }
    
    result = tmpresult;
}

template<class vobj>
void FFTforward(Lattice<vobj> &result,const Lattice<vobj> &source){
    conformable(result._grid,grid);
    conformable(source._grid,grid);
    
    pPerpProjectionForward(result,source);
    theFFT.FFT_all_dim(result,result,FFT::forward);
}

template<class vobj>
void FFTbackward(Lattice<vobj> &result,const Lattice<vobj> &source){
    conformable(result._grid,grid);
    conformable(source._grid,grid);
    
    theFFT.FFT_all_dim(result,source,FFT::backward);
    pPerpProjectionBackward(result,result);
}

template<class vobj>
void FreeWilsonOperatorInverse(Lattice<vobj> &result,const Lattice<vobj> &source){
    
    FFTforward(result,source);
    TwistedWilsonPropagator(result,result);
    FFTbackward(result,result);
}

template<class vobj>
void TwistedWilsonPropagator(Lattice<vobj> &result,const Lattice<vobj> &source){
    
    Complex mci(0.,-1.);
    Lattice<vobj> tmp(grid), tmpresult(grid);
    
    Gamma::Algebra Gmu [] = {
        Gamma::Algebra::GammaX,
        Gamma::Algebra::GammaY,
        Gamma::Algebra::GammaZ,
        Gamma::Algebra::GammaT
    };
    
    tmpresult = zero;
    
    for (int mu=0; mu<Nd; mu++) {
        tmp = QCD::Gamma(Gmu[mu])*source;
        TwistMult(tmp,pbarmu[mu],tmp);
        tmpresult += tmp;
    }
    
    tmpresult *= mci;
    TwistMult(tmp,halfphatsq,source);
    
    result = tmpresult + tmp;
}

};

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //  TWISTED FAST FOURIER TRANSFORM FOR GAUGE FIELD
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  
template <class Gimpl>
class TwistedGaugeFFT {

    typedef typename std::remove_reference<decltype(Gimpl::twist.omega[0])>::type TwistBase;
    typedef typename Gimpl::ScalarField LatticeScalar;
    typedef typename Gimpl::MatrixField LatticeMatrix;
    
private:
    
    GridCartesian* grid;
    FFT theFFT;
    LatticeMatrix pPerpPhase;            // exp(i (p_perp) * x)
    LatticeMatrix phatsq;            // 4 * sum_mu sin^2(p_mu/2) / alpha
    RealD alpha;
    
public:
    
    iMatrix<TwistBase, Nc> Gamma, adjGamma;
    int t1,t2;
    
    TwistedGaugeFFT(GridCartesian* grid_, RealD alpha_):
    grid(grid_),
    alpha(alpha_),
    theFFT(grid_),
    pPerpPhase(grid_),
    phatsq(grid_)
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
        
        // INITIALISE Fourier twist base
        pPerpLoop(n1,n2) {
            Gamma(n1,n2) = BuildGamma(n1,n2);
            adjGamma(n1,n2) = adj(Gamma(n1,n2));
            Gamma(n1,n2) = (1./(double)Nc)*Gamma(n1,n2);
        }
        
        FFTinitialisation(alpha_);
    }
    
    
    void FFTinitialisation(RealD alpha){
        
        LatticeScalar xmu(grid), tmp(grid);
        
        // INITIALISE phases for pPerp projection
        pPerpLoop(n1,n2) {
            tmp = zero;
            for (int mu=0; mu<Nd; mu++) {
                Real pPerp = 0;
                if(mu==t1) pPerp +=  (double)n1 * 2. * M_PI / (double)(grid->_fdimensions[mu]) / (double)Nc;
                if(mu==t2) pPerp +=  (double)n2 * 2. * M_PI / (double)(grid->_fdimensions[mu]) / (double)Nc;
                LatticeCoordinate(xmu,mu);
                tmp = tmp + pPerp*xmu;
            }
            Complex ci(0.,1.);
            tmp = exp(ci*tmp);
            pokeColour(pPerpPhase,tmp,n1,n2);
        }
        
        // INITIALISE momenta for Fourier acceleration
        pPerpLoop(n1,n2) {
            tmp = zero;
            for (int mu=0; mu<Nd; mu++) {
                LatticeCoordinate(xmu,mu);
                xmu *= M_PI / (double)(grid->_fdimensions[mu]);
                if(mu==t1) xmu +=  (double)n1 * M_PI / (double)(grid->_fdimensions[mu]) / (double)Nc;
                if(mu==t2) xmu +=  (double)n2 * M_PI / (double)(grid->_fdimensions[mu]) / (double)Nc;
                tmp = tmp + 4.*sin(xmu)*sin(xmu);
            }
            tmp = (1./alpha)*tmp;
            pokeColour(phatsq,tmp,n1,n2);
        }
        
        // give a non-zero value to the zero mode
        TwistBase tmpmat;
        peekSite(tmpmat,phatsq,std::vector<int>({0,0,0,0}));
        tmpmat()()()(0,0) = 1.;
        pokeSite(tmpmat,phatsq,std::vector<int>({0,0,0,0}));
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
    
    Lattice<vobj> tmpresult(grid);
    decltype(peekColour(source,0,0)) tmp(grid);
    
    pPerpLoop(n1,n2){
        tmp = trace(adjGamma(n1,n2)*source);
        pokeColour(tmpresult,tmp,n1,n2);
    }
    
    TwistMult(result,conjugate(pPerpPhase),tmpresult);
}

template<class vobj>
void pPerpProjectionBackward(Lattice<vobj> &result,const Lattice<vobj> &source){
    
    Lattice<vobj> tmp(grid), tmpresult(grid);
    TwistMult(tmp,pPerpPhase,source);
    
    tmpresult = zero;
    pPerpLoop(n1,n2){
        tmpresult = tmpresult + Gamma(n1,n2) * peekColour(tmp,n1,n2);
    }
    
    result = tmpresult;
}

template<class vobj>
void FFTforward(Lattice<vobj> &result,const Lattice<vobj> &source){
    conformable(result._grid,grid);
    conformable(source._grid,grid);
    
    pPerpProjectionForward(result,source);
    theFFT.FFT_all_dim(result,result,FFT::forward);
}

template<class vobj>
void FFTbackward(Lattice<vobj> &result,const Lattice<vobj> &source){
    conformable(result._grid,grid);
    conformable(source._grid,grid);
    
    theFFT.FFT_all_dim(result,source,FFT::backward);
    pPerpProjectionBackward(result,result);
}

template<class vobj>
void mult_invphatsq(Lattice<vobj> &source){
    TwistDivide(source,source,phatsq);
}

};

}
}
}
#endif
