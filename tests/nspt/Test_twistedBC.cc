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
using namespace QCDpt;

int main(int argc, char *argv[]) {
    
    Grid_init(&argc,&argv);
    
    std::vector<int> latt_size   = GridDefaultLatt();
    std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
    std::vector<int> mpi_layout  = GridDefaultMpi();
    
    GridCartesian               Grid(latt_size,simd_layout,mpi_layout);
    
    int threads = GridThread::GetThreads();
    std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;
    
    GridParallelRNG pRNG(&Grid);
    pRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));
    
    cout<<twist.twist_tensor()<<endl;
    
    double tau = 0.01;
    double alpha = -0.5*tau;
    double gftolerance = 0.1;
    
    QCDpt::LatticeGaugeField U(&Grid);
    PertVacuum(U);
    
    PertLangevin<WilsonGaugeAction<TwistedGimpl_ptR>> L(&Grid,pRNG,tau,alpha);
    
    
    // checkpointer
    CheckpointerParameters CPparams;
    CPparams.config_prefix = "NSPTckpoint_lat";
    CPparams.rng_prefix = "NSPTckpoint_rng";
    CPparams.saveInterval = 25;
    CPparams.format = "IEEE64BIG";
    BinaryHmcCheckpointer<TwistedGimpl_ptR> CP(CPparams);
    // is it needed?
    GridSerialRNG sRNG;
    sRNG.SeedFixedIntegers(std::vector<int>({1,2,3,4}));
    
    
    // gnuplot:
    // plot for [i=0:6] "plaquette.txt" every 7::i u 1
    ofstream plaqfile("plaquette.txt");
    plaqfile.precision(30);
    plaqfile << scientific;
    
    PRealD plaq;
    for (int i=0; i<25; i++) {
        L.QuenchRKStep(U);
        plaq = WilsonLoops<TwistedGimpl_ptR>::avgPlaquette(U);
        for (int k=0; k<Np; k++) plaqfile << plaq(k) << endl;
//        if (i%25==0) cout<<WilsonLoops<TwistedGimpl_ptR>::avgPlaquette(U)<<endl;
        CP.TrajectoryComplete(i,U,sRNG,pRNG);
    }
    
    L.LandauGF(U, gftolerance);
    
    Grid_finalize();
    return EXIT_SUCCESS;
    
}
