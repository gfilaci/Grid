/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: ./tests/nspt/Test_twistedFFT.cc

Copyright (C) 2015-2017

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
#include <Grid/Grid.h>

using namespace std;
using namespace Grid;
using namespace QCD;
using namespace QCDpt;

const double tolerance = 1e-15;

#define coordloop(n,mu) for (n[mu]=0; n[mu]<Grid._fdimensions[mu]; n[mu]++)
#define perploop(n) for (n=0; n<Nc; n++)

int main(int argc, char *argv[]) {
    
    Grid_init(&argc,&argv);
    
    std::vector<int> latt_size   = GridDefaultLatt();
    std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
    std::vector<int> mpi_layout  = GridDefaultMpi();
    
    GridCartesian               Grid(latt_size,simd_layout,mpi_layout);
    
    int threads = GridThread::GetThreads();
    std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;
    
    typedef TwistedGimpl_ptR   gimpl;
    gimpl::MatrixField M(&Grid), Mtilde(&Grid), Mnew(&Grid);
    
    vector<Complex> boundary_phases = {-1.,1.,1.,1.};
    
    TwistedFFT<gimpl> TheFFT(&Grid,boundary_phases);
    
    int volume = 1;
    for (int d=0; d<Nd; d++) {
        volume *= Grid._fdimensions[d];
    }
    
    cout << GridLogMessage << "Testing the twisted FFT";
    
    // build plane wave with momentum n
    vector<int> n = {1,1,1,1};
    int n1tilde = 1;
    int n2tilde = 0;
    
    perploop(n1tilde)
    perploop(n2tilde)
    coordloop(n,0)
    coordloop(n,1)
    coordloop(n,2)
    coordloop(n,3)
    {
    
    gimpl::ScalarField xmu(&Grid), tmp(&Grid);
    tmp = zero;
    for (int mu=0; mu<Nd; mu++) {
        Real qmu = TheFFT.boundary_exp[mu] / (double)(Grid._fdimensions[mu]);
        qmu +=  (double)n[mu] * 2. * M_PI / (double)(Grid._fdimensions[mu]);
        if(mu==TheFFT.t1) qmu +=  (double)n1tilde * 2. * M_PI / (double)(Grid._fdimensions[mu]) / (double)Nc;
        if(mu==TheFFT.t2) qmu +=  (double)n2tilde * 2. * M_PI / (double)(Grid._fdimensions[mu]) / (double)Nc;
        LatticeCoordinate(xmu,mu);
        tmp = tmp + qmu*xmu;
    }
    Complex ci(0.,1.);
    tmp = (double)Nc*exp(ci*tmp);
    
    M = tmp*TheFFT.Gamma(n1tilde,n2tilde);
    
    TheFFT.FFTforward(Mtilde,M);
    TheFFT.FFTbackward(Mnew,Mtilde);
    
    typename gimpl::MatrixField::vector_object::scalar_object ss;
    peekSite(ss,Mtilde,n);
    ss()()()(n1tilde,n2tilde) -= (double)volume*(double)Nc;
    pokeSite(ss,Mtilde,n);
    
    Mnew -= M;
    
    if(norm2(Mtilde)<tolerance && norm2(Mnew)<tolerance){
        cout<<".";
        fflush(stdout);
    }
    else{
        cout << endl << GridLogMessage << "Test FAILED" << endl;
        exit(EXIT_FAILURE);
    }
    
    
    }
    cout << endl << GridLogMessage << "Test SUCCEEDED" << endl;
    
    Grid_finalize();
    return EXIT_SUCCESS;
    
}
