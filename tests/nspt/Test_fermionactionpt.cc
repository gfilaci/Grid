/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: ./tests/Test_fermionactionpt.cc

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
    
    std::vector<int> latt_size({4,4,4,4});
    std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
    std::vector<int> mpi_layout  = GridDefaultMpi();
    
    GridCartesian               Grid(latt_size,simd_layout,mpi_layout);
    GridRedBlackCartesian     RBGrid(latt_size,simd_layout,mpi_layout);
    
    int threads = GridThread::GetThreads();
    std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;
    
    GridParallelRNG pRNG(&Grid);
    pRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));
    
    PRealD mass = zero;
    int Nf = 2;
    WilsonImplParams Params;
    Params.boundary_phases = {-1.,1.,1.,1};
    
    
    QCDpt::LatticeGaugeField U(&Grid);
    StochasticFermionAction<PWilsonSmellImplR> FermionAction(pRNG,&Grid,&RBGrid,mass,Params,Nf);
    FermionAction.deriv(U,U);
    
    
//    AccessTypes<Action, Representations<FundamentalRep_pt<Nc>>>::VectorCollection actions_hirep;
//    vector<Action<PWilsonSmellImplR::GaugeField>*>& actions(std::get<0>(actions_hirep));
    
    
    ActionLevel<PWilsonSmellImplR::GaugeField,PNoHirep> Level0;
    Level0.push_back(&FermionAction);
    ActionSet<PWilsonSmellImplR::GaugeField,PNoHirep> TheActions;
//    TheActions.push_back(Level0);
//    for(int i=0; i<TheActions.size(); i++) TheActions[i].actions[0]->deriv(U,U);
    
    
    
//    FFT theFFT(&Grid);
//    
//    QCDpt::LatticeGaugeField Ucopy(&Grid);
//    random(pRNG,U);
//    Ucopy = U;
//    theFFT.FFT_all_dim(U,U,FFT::backward);
//    theFFT.FFT_all_dim(U,U,FFT::forward);
//    U-=Ucopy;
//    cout<<U<<endl;
    
    Grid_finalize();
    return EXIT_SUCCESS;
    
}
