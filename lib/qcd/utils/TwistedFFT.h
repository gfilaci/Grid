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
namespace QCD {
namespace QCDpt {

template <class Gimpl>
class TwistedFFT {

    typedef typename std::remove_reference<decltype(Gimpl::twist.omega[0])>::type TwistBase;
    typedef typename Gimpl::ScalarField LatticeScalar;
    typedef typename Gimpl::MatrixField LatticeMatrix;
    
private:
    iMatrix<TwistBase, Nc> Gamma, adjGamma;
    int t1,t2;
    std::vector<Real> boundary_exp;
    GridCartesian* grid;
    FFT theFFT;
    LatticeMatrix pPerpPhase;
    
public:
    
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
        // The FFTW library has FFTW_FORWARD = -1 by default.
        // The following if statement is needed when
        // somehow a different convention is used (hopefully never...).
        if(FFTW_FORWARD==+1) pPerpPhase = conjugate(pPerpPhase);
        
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
    
    if(testpassed) std::cout << GridLogMessage << "Orthonormality test PASSED" << std::endl;
    else std::cout << GridLogMessage << "Orthonormality test FAILED" << std::endl;
}

template<class vobj>
void pPerpProjectionForward(Lattice<vobj> &result,const Lattice<vobj> &source){

}

template<class vobj>
void pPerpProjectionBackward(Lattice<vobj> &result,const Lattice<vobj> &source){

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
