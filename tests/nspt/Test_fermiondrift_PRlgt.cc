/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: ./tests/nspt/Test_fermiondrift_PRlgt.cc

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
    
    if(Np!=7||Nc!=3||Nd!=4||istwisted(0)==true||istwisted(1)==false||istwisted(2)==false||istwisted(3)==true){
    cout<<"This test has to be compiled with Nc=3, Nd=4, Np=7 with twist on the plane 12."<<endl;
    return EXIT_FAILURE;
    }
    
    // some matrices
    iMatrix<Complex,Nc> AA,BB,CC;
    Complex im(0.,1.);
    
    AA(0,0) = 0.0 + 9.0 * im;
    AA(0,1) = 2.0 + 8.0 * im;
    AA(0,2) = 4.0 + 7.0 * im;
    AA(1,0) = 6.0 + 6.0 * im;
    AA(1,1) = 8.0 + 5.0 * im;
    AA(1,2) = 7.0 + 4.0 * im;
    AA(2,0) = 5.0 + 3.0 * im;
    AA(2,1) = 3.0 + 2.0 * im;
    AA(2,2) = 1.0 + 1.0 * im;
    
    BB(0,0) = 0.1 + 5.9 * im;
    BB(0,1) = 0.2 + 5.9 * im;
    BB(0,2) = 0.3 + 5.8 * im;
    BB(1,0) = 0.4 + 8.8 * im;
    BB(1,1) = 0.5 + 8.7 * im;
    BB(1,2) = 0.6 + 8.7 * im;
    BB(2,0) = 0.7 + 3.6 * im;
    BB(2,1) = 0.8 + 3.6 * im;
    BB(2,2) = 0.9 + 3.5 * im;
    
    CC = AA*BB;
    
    
    
    // gauge field
    QCDpt::LatticeGaugeField U(&Grid),dSdU(&Grid);
    typename PWilsonSmellImplR::SOimpl::Gimpl::Field Uso(&Grid), Uforce(&Grid);
    
    // fermion field
    QCDpt::LatticeSpinColourMatrix Xi(&Grid);
    std::vector<typename PWilsonSmellImplR::SOimpl::FermionField> psi;
    for (int i=0; i<Np-2; i++) {
        psi.push_back(typename PWilsonSmellImplR::SOimpl::FermionField(&Grid));
    }
    PRealD mass = zero;
    mass(2) = -1.;
    mass(4) = -5.;
    mass(6) = -10.;
    // antiperiodic fermions
    WilsonImplParams Params;
    Params.boundary_phases = {-1.,1.,1.,1};

    
    // initialise U
    QCDpt::LorentzPertColourMatrix ltmpmat;
    for (int x0=0; x0<4; x0++) {
        for (int x1=0; x1<4; x1++) {
            for (int x2=0; x2<4; x2++) {
                for (int x3=0; x3<4; x3++) {
                    zeroit(ltmpmat);
                    int x = x0*4*4*4 + x1*4*4 + x2*4 +x3;
                    for (int mu=0; mu<Nd; mu++) {
                        for (int n=0; n<Nc; n++) ltmpmat(mu)()(0)(n,n) = 1.;
                        for (int k=0; k<Np-1; k++) {
                            ltmpmat(mu)()(k+1) = (double)x*CC + (double)mu*AA + (double)k*BB;
                        }
                    }
                    pokeSite(ltmpmat,U,std::vector<int>({x0,x1,x2,x3}));
                }
            }
        }
    }
    
    // initialise Xi
    QCDpt::SpinColourMatrix tmpmat;
    for (int x0=0; x0<4; x0++) {
        for (int x1=0; x1<4; x1++) {
            for (int x2=0; x2<4; x2++) {
                for (int x3=0; x3<4; x3++) {
                    zeroit(tmpmat);
                    int x = x0*4*4*4 + x1*4*4 + x2*4 +x3;
                    for (int alpha=0; alpha<Ns; alpha++) {
                        tmpmat()(alpha)() = (double)x*AA + (double)alpha*BB + (double)(x+alpha)*CC;
                    }
                    pokeSite(tmpmat,Xi,std::vector<int>({x0,x1,x2,x3}));
                }
            }
        }
    }
    
    int Nf = 2;
    StochasticFermionAction<PWilsonSmellImplR> FermionAction(&pRNG,&Grid,&RBGrid,mass,Params,Nf);
    
    FermionAction.invM(psi,U,Xi);
    
    // compute force
    dSdU = zero;
    
    for (int n=0; n<Np-2; n++) {
          Uforce = zero;
          for (int j=0; j<=n; j++) {
              FermionAction.Dw[n-j].MDeriv(Uso, psi[j], Xi, DaggerYes);
              Uforce += Uso;
          }
          // the "+2" shift is due to the 1/beta factor in front of the fermion drift.
          // the first two orders are already set to zero.
          pokePert(dSdU,Uforce,n+2);
      }
    
    dSdU = Ta(dSdU);
    dSdU *= -0.005*(double)Nf/(double)Nc;
    
    ofstream file("Grid_fermiondrift.txt");
    file.precision(10);
    file << scientific;
    for (int x0=0; x0<4; x0++) {
        for (int x1=0; x1<4; x1++) {
            for (int x2=0; x2<4; x2++) {
                for (int x3=0; x3<4; x3++) {
                    peekSite(ltmpmat,dSdU,std::vector<int>({x0,x1,x2,x3}));
                    for (int alpha=0; alpha<Ns; alpha++) {
                        for (int k=1; k<Np; k++) {
                            for (int n1=0; n1<Nc; n1++) {
                                for (int n2=0; n2<Nc; n2++) {
                                    file << real(ltmpmat(alpha)()(k)(n1,n2)) << endl;
                                    file << imag(ltmpmat(alpha)()(k)(n1,n2)) << endl;
                                }
                            }
                        }
                    }
                    
                }
            }
        }
    }
    file.close();
    
    Grid_finalize();
    return EXIT_SUCCESS;
    
}


///////////////////////////////////////////////////////
    /////////////////////////////////////////////
    // TEST AGAINST PRlgt
    /////////////////////////////////////////////
    // note that gamma_peter = U gamma_franz U
    // with
    //  U = | sigma2    0   |
    //      |   0    sigma2 |
    /////////////////////////////////////////////
/*
The results in Grid_fermiondrift.txt have to be compared with tests/nspt/PRlgt_fermiondrift.txt.
That file is produced by the PRlgt code with the following options:

- gauge group SU(3)
- twist directions Y,Z
- twist matrices respectively
( -0.5, -0.866025 ) 	( 0, 0 ) 	( 0, 0 ) 	
( 0, 0 ) 	( 1, 0 ) 	( 0, 0 ) 	
( 0, 0 ) 	( 0, 0 ) 	( -0.5, 0.866025 ) 	

( 0, 0 ) 	( 1, 0 ) 	( 0, 0 ) 	
( 0, 0 ) 	( 0, 0 ) 	( 1, 0 ) 	
( 1, 0 ) 	( 0, 0 ) 	( 0, 0 )
- antiperiodic bc for fermions in direction X
- critical masses mcpt[6 + 1] = {4., 0., -1., 0., -5., 0., -10.};
- perturbative order allocORD = 6
- number of flavours NF = 2
- volume L^4 = 4^4
- tau_g = 0.005

The following code computes the drift from a given stochastic noise.
It has be compiled and linked to the library prlgt.a.
The easiest way to reproduce the output is to insert the code in the main function of Unquenched.cc,
after everything is initialised.
Then it is important to comment the Xi.gauss random generation in the function fermion_drift of nspt.cc.
It must also be checked that the zero mode subtraction of the stochastic noise (STOCTRACE_ZEROMOM_SUBTR)
is NOT activated.
*/
///////////////////////////////////////////////////////
/*

    SU3 AA,BB,CC;
    Cplx im(0.,1.);
    
    AA.whr[0] = 0.0 + 9.0 * im;
    AA.whr[1] = 2.0 + 8.0 * im;
    AA.whr[2] = 4.0 + 7.0 * im;
    AA.whr[3] = 6.0 + 6.0 * im;
    AA.whr[4] = 8.0 + 5.0 * im;
    AA.whr[5] = 7.0 + 4.0 * im;
    AA.whr[6] = 5.0 + 3.0 * im;
    AA.whr[7] = 3.0 + 2.0 * im;
    AA.whr[8] = 1.0 + 1.0 * im;
    
    BB.whr[0] = 0.1 + 5.9 * im;
    BB.whr[1] = 0.2 + 5.9 * im;
    BB.whr[2] = 0.3 + 5.8 * im;
    BB.whr[3] = 0.4 + 8.8 * im;
    BB.whr[4] = 0.5 + 8.7 * im;
    BB.whr[5] = 0.6 + 8.7 * im;
    BB.whr[6] = 0.7 + 3.6 * im;
    BB.whr[7] = 0.8 + 3.6 * im;
    BB.whr[8] = 0.9 + 3.5 * im;
    
    CC = AA*BB;
    
    // initialise Umu, Fmu
    for(int i = 0; i < Umu.Z->Size; i++)
        for(int mu = 0; mu < dim; mu++){
            Umu.W[i].U[mu].flag = 1.;
            Fmu.W[i].U[mu].zero();
            Fmu.W[i].U[mu].flag = 0.;
            for (int k=0; k<PTORD; k++){
                Umu.W[i].U[mu].ptU[k] = i*CC + mu*AA + k*BB;
            }
        }
    twist_boundary(Umu);
    
    memset((void*)Xi.psi, 0, LL.Size*sizeof(SpinColorSmell));
    
    // initialise Xi
    for (int alpha=0; alpha<dim; alpha++) {
        for (int x=0; x<LL.Size; x++) {
            Xi.psi[x].psi[alpha] = x*AA + alpha*BB + (x+alpha)*CC;
        }
    }
    
    // gamma_Franz -> gamma_Peter
    for (int x=0; x<LL.Size; x++) {
        SUX app = Xi.psi[x].psi[0];
        Xi.psi[x].psi[0] = -im*Xi.psi[x].psi[1];
        Xi.psi[x].psi[1] = im*app;
        app = Xi.psi[x].psi[2];
        Xi.psi[x].psi[2] = -im*Xi.psi[x].psi[3];
        Xi.psi[x].psi[3] = im*app;
    }
    
    fermion_drift(Umu, Fmu, Pmu, Xi, 1.);
    
    // print Pmu_app
    ofstream file("PRlgt_fermiondrift.txt");
    file.precision(10);
    file << scientific;
    for (int y0=0; y0<4; y0++) {
        for (int y1=0; y1<4; y1++) {
            for (int y2=0; y2<4; y2++) {
                for (int y3=0; y3<4; y3++) {
                    int y = y0*4*4*4 + y1*4*4 + y2*4 + y3;
                    for (int alpha=0; alpha<dim; alpha++) {
                        for (int k=0; k<PTORD; k++) {
                            for (int n1=0; n1<NC; n1++) {
                                for (int n2=0; n2<NC; n2++) {
                                    file << Fmu.W[y].U[alpha].ptU[k].whr[n1*NC+n2].re << endl;
                                    file << Fmu.W[y].U[alpha].ptU[k].whr[n1*NC+n2].im << endl;
                                }
                            }
                        }
                    }
                    
                }
            }
        }
    }
    file.close();
    
    exit(0);

*/
