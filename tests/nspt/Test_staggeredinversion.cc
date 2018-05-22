/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: ./tests/nspt/Test_staggeredinversion.cc

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

#define testadjoint

using namespace std;
using namespace Grid;
using namespace QCD;
using namespace QCDpt;

///////////////////////////////////////////////
// Read parameters from command line
///////////////////////////////////////////////

int main(int argc, char *argv[]) {
    
    Grid_init(&argc,&argv);
    
    
    std::vector<int> latt_size = GridDefaultLatt();
    std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
    std::vector<int> mpi_layout  = GridDefaultMpi();
    
    GridCartesian               Grid(latt_size,simd_layout,mpi_layout);
    GridRedBlackCartesian     RBGrid(&Grid);

    int threads = GridThread::GetThreads();
    std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;
    
    GridParallelRNG pRNG(&Grid);
    pRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));
    
    typedef TwistedGimpl_ptR   gimpl;
#ifndef testadjoint
    typedef PStaggeredSmellImplR  fimpl;
#else
    typedef PStaggeredAdjointImplR  fimpl;
#endif
    
    gimpl::GaugeField U(&Grid);
    PertRandom(pRNG,U);
    
    typedef typename fimpl::SOimpl::FermionField SOFermionField;
    std::vector<SOFermionField> psi2;
    psi2.reserve(Np);
    for (int i=0; i<Np; i++) psi2.push_back(SOFermionField(&Grid));
    
    RealD mass = 1.;
    int Nf = 2;
    WilsonImplParams Params;
    Params.boundary_phases = {-1.,1.,1.,1};
    
#ifndef testadjoint
    StochasticStaggeredAction<fimpl> FA(&pRNG,&Grid,&RBGrid,mass,Params,Nf);
#else
    StochasticAdjointStaggeredAction<fimpl> FA(&pRNG,&Grid,&RBGrid,mass,Params,Nf);
#endif
    
    // ********** //
    // START TEST //
    // ********** //
    
    // invM works on Npf orders
    int Npf = Np - 2;
    
    random(pRNG,FA.Xi);
#ifdef testadjoint
    FA.Xi = Ta(FA.Xi);
#endif
    cout<<endl<<"Norm of random Xi:"<<endl;
    cout<<norm2(FA.Xi)<<endl;
    
    cout<<endl<<"Norm of gauge field U:"<<endl;
    cout<<Pnorm2(U)<<endl;
    
    FA.invM(FA.psi,U,FA.Xi);
    cout<<endl<<"Norms after invM:"<<endl;
    for (int i=0; i<Npf; i++) cout<<norm2(FA.psi[i])<<endl;
    
    FA.M(psi2,U,FA.psi);
    cout<<endl<<"Norms after M:"<<endl;
    for (int i=0; i<Npf; i++) cout<<norm2(psi2[i])<<endl;
    
    psi2[0] -= FA.Xi;
    cout<<endl<<"Norms after subtraction:"<<endl;
    for (int i=0; i<Npf; i++) cout<<norm2(psi2[i])<<endl;
    
    Grid_finalize();
    return EXIT_SUCCESS;
    
}
