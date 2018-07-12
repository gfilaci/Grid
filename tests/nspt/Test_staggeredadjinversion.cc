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
    typedef PStaggeredAdjointImplR  fimpl;
    
    gimpl::GaugeField U(&Grid);
    PertRandom(pRNG,U);
	
    typedef typename fimpl::FermionField FermionField;
    FermionField psi2(&Grid);
	
    RealD mass = 0.;
    int Nf = 2;
    WilsonImplParams Params;
    Params.boundary_phases = {-1.,1.,1.,1.};
	
    StochasticAdjointStaggeredAction<fimpl> FA(&pRNG,&Grid,&RBGrid,mass,Params,Nf);
	
    // ********** //
    // START TEST //
    // ********** //
    
    // invM works on Npf orders
    int Npf = Np - 2;
	
	decltype(peekPert(FA.Xi,0)) Xitmp(&Grid);
    random(pRNG,Xitmp);
	
	pokePert(FA.Xi,Xitmp,0);
	
    cout<<endl<<"Norm of random Xi:"<<endl;
    cout<<Pnorm2(FA.Xi)<<endl;
    
    cout<<endl<<"Norm of gauge field U:"<<endl;
    cout<<Pnorm2(U)<<endl;
    
    FA.invM(FA.psi,U,FA.Xi);
    cout<<endl<<"Norms after invM:"<<endl;
    cout<<Pnorm2(FA.psi)<<endl;
	
    FA.M(psi2,U,FA.psi);
	// invM works on Npf orders, so the 2 left
	// are set to zero manually
	FA.psi_so = zero;
	for(int i=Npf; i<Np; i++) pokePert(psi2,FA.psi_so,i);
    cout<<endl<<"Norms after M:"<<endl;
    cout<<Pnorm2(psi2)<<endl;

    psi2 -= FA.Xi;
    cout<<endl<<"Norms after subtraction:"<<endl;
    cout<<Pnorm2(psi2)<<endl;
	
    Grid_finalize();
    return EXIT_SUCCESS;
    
}
