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
    
    std::vector<int> latt_size({4,4,4,4});
    std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
    std::vector<int> mpi_layout  = GridDefaultMpi();
    
    GridCartesian               Grid(latt_size,simd_layout,mpi_layout);
    GridRedBlackCartesian     RBGrid(latt_size,simd_layout,mpi_layout);
    
    int threads = GridThread::GetThreads();
    std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;
    
    GridParallelRNG pRNG(&Grid);
    pRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));
    RealD massstd=0;
    PRealD mass;
    zeroit(mass);
    mass(2) = -1.;
    mass(4) = -5.;
    mass(6) = -10.;
    
    // standard gauge field, set to identity
    QCD::LatticeGaugeField Ustd(&Grid);
    parallel_for(int ss=0;ss<Grid.oSites();ss++){
        for (int alpha=0; alpha<Ns; alpha++) {
            for (int mu=0; mu<Nd; mu++) {
                for (int nn=0; nn<Nc; nn++) {
                    Ustd._odata[ss](mu)()(nn,nn) = 1.;
                }
            }
        }
    }
    // standard fermion fields, vectors in colour space
    QCD::LatticeFermion psistd(&Grid),psi2std(&Grid);
    // standard Wilson Dirac operator
    WilsonFermion<QCD::WilsonImpl<vComplex,FundamentalRepresentation,CoeffReal>> Dwstd(Ustd,Grid,RBGrid,massstd);
    
    
    // perturbative gauge field, all matrices set to the identity (-> TBC are irrelevant)
    QCDpt::LatticeGaugeField U(&Grid);
    PertVacuum(U);
    parallel_for(int ss=0;ss<Grid.oSites();ss++){
        for (int mu=0; mu<Nd; mu++) {
            for (int kk=0; kk<Np; kk++) {
                for (int nn=0; nn<Nc; nn++) {
                    (U._odata[ss])(mu)()(kk)(nn,nn) = 1.0;
                }
            }
        }
    }
    // perturbative fermion fields, matrices in colour space
    QCDpt::LatticeFermion psi(&Grid),psi2(&Grid),psicmp(&Grid);
    // antiperiodic fermions
    WilsonImplParams Params;
    Params.boundary_phases = {-1.,1.,1.,1};
    // perturbative Wilson Dirac operator
    WilsonFermion<QCDpt::PWilsonSmellImpl<vComplex,FundamentalRepresentation,CoeffReal>> Dw(U,Grid,RBGrid,mass,Params);
    
    
    // random field
    random(pRNG,psi);
    
    // apply M
    Dw.M(psi,psi2);
    
    for (int kk=0; kk<Np; kk++) {
        for (int ncol=0; ncol<Nc; ncol++) {
            
            // load ncol column of the smell matrix onto psistd
            parallel_for(int ss=0;ss<Grid.oSites();ss++){
                for (int alpha=0; alpha<Ns; alpha++) {
                    for (int nn=0; nn<Nc; nn++) {
                        (psistd._odata[ss])()(alpha)(nn) = (psi._odata[ss])()(alpha)(kk)(nn,ncol);
                    }
                }
            }
            
            // apply M
            Dwstd.M(psistd,psi2std);
            
            // load psi2std onto ncol column of the smell matrix
            parallel_for(int ss=0;ss<Grid.oSites();ss++){
                for (int alpha=0; alpha<Ns; alpha++) {
                    for (int nn=0; nn<Nc; nn++) {
                        (psicmp._odata[ss])()(alpha)(kk)(nn,ncol) = (psi2std._odata[ss])()(alpha)(nn);
                    }
                }
            }
            
        }
    }
    
    // keep into account perturbative structure of Dirac operator
    for (int kk=0; kk<Np; kk++) {
        parallel_for(int ss=0;ss<Grid.oSites();ss++){
            for (int alpha=0; alpha<Ns; alpha++) {
                for (int ord=0; ord<=kk; ord++) {
                    (psi2._odata[ss])()(alpha)(kk) -= (psi2._odata[ss])()(alpha)(ord);
                }
            }
        }
    }
    
    cout<<Pnorm2(psi2)<<endl;
    
    
    
    
    /////////////////////////////////////////////
    // TEST AGAINST PRlgt
    /////////////////////////////////////////////
    // note that gamma_peter = U gamma_franz U
    // with
    //  U = | sigma2    0   |
    //      |   0    sigma2 |
    /////////////////////////////////////////////
    
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
    
    // initialise psi
    SpinPertColourMatrix tmpmat;
    for (int x0=0; x0<4; x0++) {
        for (int x1=0; x1<4; x1++) {
            for (int x2=0; x2<4; x2++) {
                for (int x3=0; x3<4; x3++) {
                    zeroit(tmpmat);
                    int x = x0*4*4*4 + x1*4*4 + x2*4 +x3;
                    for (int alpha=0; alpha<Ns; alpha++) {
                        for (int k=0; k<Np; k++) {
                            tmpmat()(alpha)(k) = (double)x*AA + (double)alpha*BB + (double)k*CC;
                        }
                    }
                    pokeSite(tmpmat,psi,std::vector<int>({x0,x1,x2,x3}));
                }
            }
        }
    }
    
    Dw.M(psi,psi2);
//    cout<<psi2<<endl;
    cout<<Dw.mass<<endl;
    Grid_finalize();
    return EXIT_SUCCESS;
    
}
