/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: ./tests/Test_twistedBC.cc

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

int main(int argc, char *argv[]) {
    
    Grid_init(&argc,&argv);
    
//    std::vector<int> latt_size   = GridDefaultLatt();
    std::vector<int> latt_size({4,4,4,4});
    std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
    std::vector<int> mpi_layout  = GridDefaultMpi();
    
    GridCartesian               Grid(latt_size,simd_layout,mpi_layout);
    GridRedBlackCartesian     RBGrid(latt_size,simd_layout,mpi_layout);
    
    int threads = GridThread::GetThreads();
    std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;
    
    GridParallelRNG pRNG(&Grid);
    pRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));
    
    QCD::LatticeGaugeField U(&Grid);
    QCD::LatticeFermion psi(&Grid),psi2(&Grid);
    
    gaussian(pRNG,U);
    gaussian(pRNG,psi);
    
    RealD mass=-4.0;
//    WilsonFermionR Dw(U,Grid,RBGrid,mass);
//    WilsonFermion<WilsonImplR> Dw(U,Grid,RBGrid,mass);
    WilsonFermion<QCDpt::PWilsonSmellImpl<vComplex,FundamentalRepresentation,CoeffReal>> Dw(U,Grid,RBGrid,mass);
    
    
    Dw.M(psi,psi2);
    
    
    cout<<psi<<endl;
    cout<<psi2<<endl;
    
    
    Grid_finalize();
    return EXIT_SUCCESS;
    
}
