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
using namespace Grid::QCD::QCDpt;

int main(int argc, char *argv[]) {
    
    Grid_init(&argc,&argv);
    
    std::vector<int> latt_size   = GridDefaultLatt();
    std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
    std::vector<int> mpi_layout  = GridDefaultMpi();
    
    GridCartesian               Grid(latt_size,simd_layout,mpi_layout);
    GridRedBlackCartesian     RBGrid(latt_size,simd_layout,mpi_layout);
    
    int threads = GridThread::GetThreads();
    std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;
    
    std::vector<int> seeds({1,2,3,4});
    
    GridParallelRNG          pRNG(&Grid);
    pRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));
    
    QCDpt::LatticeGaugeField U(&Grid);
    QCDpt::LatticeGaugeField F(&Grid);
    gaussian(pRNG,U);
    U = ProjectOnGroup(U);
    
    double beta = 1.0;
//    WilsonGaugeActionR Action(beta);
//    Action.deriv(U,F);

    RealD factor = 0.5 * beta / RealD(Nc);
//
//    
//    for (int mu = 0; mu < Nd; mu++) {
//        WilsonLoops<Gimpl>::StapleMult(dSdU_mu, U, mu);
//        F = Ta(dSdU_mu) * factor;
//        PokeIndex<LorentzIndex>(dSdU, dSdU_mu, mu);
//    }
//    
    
    
    Grid_finalize();
    return EXIT_SUCCESS;
    
}
