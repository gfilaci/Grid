/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: ./tests/nspt/Test_fermionactionpt.cc

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
    
    
    typedef TwistedGimpl_ptR   gimpl;
    typedef PWilsonSmellImplR  fimpl;
    gimpl::GaugeField U(&Grid);
    PertVacuum(U);
    
    
    WilsonGaugeAction<gimpl> GaugeAction;
    ActionLevel<gimpl::GaugeField,PNoHirep> GaugeLevel;
    GaugeLevel.push_back(&GaugeAction);
    
    
    PRealD mass = zero;
    int Nf = 2;
    WilsonImplParams Params;
    Params.boundary_phases = {-1.,1.,1.,1};
    
    StochasticFermionAction<fimpl> FermionAction(&pRNG,&Grid,&RBGrid,mass,Params,Nf);
    ActionLevel<fimpl::GaugeField,PNoHirep> FermionLevel;
    FermionLevel.push_back(&FermionAction);
    
    
    double tau = 0.01;
    double alpha = -0.5*tau;
    PertLangevin<gimpl> L(&Grid,&pRNG,tau,alpha);
    L.TheActions.push_back(GaugeLevel);
    L.TheActions.push_back(FermionLevel);
    
    
    // gnuplot:
    // plot for [i=0:6] "plaquette.txt" every 7::i u 1
    ofstream plaqfile("plaquette.txt");
    plaqfile.precision(30);
    plaqfile << scientific;
    
    PRealD plaq;
    for (int i=0; i<5000; i++) {
        L.QuenchRKStep(U);
        plaq = WilsonLoops<TwistedGimpl_ptR>::avgPlaquette(U);
        for (int k=0; k<Np; k++) plaqfile << plaq(k) << endl;
    }
    
    Grid_finalize();
    return EXIT_SUCCESS;
    
}
