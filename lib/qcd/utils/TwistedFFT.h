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

#define pperploop(n1,n2) for(int n1=0; n1<Nc; n1++) for(int n2=0; n2<Nc; n2++)

namespace Grid {
namespace QCD {
namespace QCDpt {

template <class Gimpl>
class TwistedFFT {

public:
    
    typedef typename std::remove_reference<decltype(Gimpl::twist.omega[0])>::type TwistBase;
    iMatrix<TwistBase, Nc> Gamma, adjGamma;
    int t1,t2;
    std::vector<Complex> boundary_phases;
    
    TwistedFFT()
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
        
        // initialise Fourier twist base
        pperploop(n1,n2) {
            Gamma(n1,n2) = BuildGamma(n1,n2);
            adjGamma(n1,n2) = adj(Gamma(n1,n2));
            Gamma(n1,n2) = (1./(double)Nc)*Gamma(n1,n2);
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
    pperploop(n1,n2){
        pperploop(m1,m2){
            double tmp = abs(TensorRemove(trace(adjGamma(m1,m2)*Gamma(n1,n2))));
            if(m1==n1 && m2==n2) tmp -= 1.;
            if(std::abs(tmp)>precision) testpassed = false;
        }
    }
    
    if(testpassed) std::cout << GridLogMessage << "Orthonormality test PASSED" << std::endl;
    else std::cout << GridLogMessage << "Orthonormality test FAILED" << std::endl;
}

};

}
}
}
#endif
