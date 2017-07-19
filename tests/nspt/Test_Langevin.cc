/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: ./tests/Test_nspt.cc

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
using namespace QCDpt;

int main(int argc, char *argv[]) {
    
    Grid_init(&argc,&argv);
    
    std::vector<int> latt_size   = GridDefaultLatt();
    std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
    std::vector<int> mpi_layout  = GridDefaultMpi();
    
    GridCartesian               Grid(latt_size,simd_layout,mpi_layout);
    
    int threads = GridThread::GetThreads();
    std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;
    
    GridParallelRNG          pRNG(&Grid);
    pRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));
    
    
    double beta = 1.;
    double tau = 0.0025*2*2;
    double stau = std::sqrt(tau);
    
    WilsonGaugeAction<PeriodicGaugeImpl<GimplTypes_ptR>> Action(beta);
    
    QCDpt::LatticeGaugeField U(&Grid);
    QCDpt::LatticeGaugeField F(&Grid);
    QCDpt::LatticeLorentzColourMatrix noise(&Grid);
    QCDpt::LatticeColourMatrix tmp(&Grid);

    // cold start
    zeroit(U);
    
    // random start
//    random(pRNG,U);

    U = ProjectOnGroup(U);
    
    for (int i=0; i<10000; i++) {
    
        //noise
        for (int mu=0; mu<Nd; mu++) {
            QCDpt::SU<Nc>::GaussianFundamentalLieAlgebraMatrix(pRNG, tmp, M_SQRT2);
            pokeLorentz(noise,tmp,mu);
        }
        Action.deriv(U,F);
        F = - stau * F;
        F = stau * AddToOrd(1,F,noise);
        U = Exponentiate(F) * U;
//        if (i%100==0) {cout<<Action.S(U)<<endl;cout<<"\t"<<Pnorm2(U)<<endl;}
        if (i%100==0) cout<<WilsonLoops<PeriodicGaugeImpl<GimplTypes_ptR>>::avgPlaquette(U)<<endl;
    }
    
    
    Grid_finalize();
    return EXIT_SUCCESS;
    
}
