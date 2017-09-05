/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: ./tests/nspt/Test_wilsonpt.cc

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
    
    // masses and antiperiodicity
    RealD massstd = 0.;
    PRealD mass = zero;
    WilsonImplParams Params;
    Params.boundary_phases = {-1.,1.,1.,1};
    
    
    // standard gauge field, set to identity
    QCD::LatticeGaugeField Ustd(&Grid);
    zeroit(Ustd);
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
    WilsonFermion<QCD::WilsonImpl<vComplex,FundamentalRepresentation,CoeffReal>> Dwstd(Ustd,Grid,RBGrid,massstd,Params);
    
    
    
    // perturbative gauge field, all matrices set to the identity (-> TBC are irrelevant)
    QCDpt::LatticeGaugeField U(&Grid);
    parallel_for(int ss=0;ss<Grid.oSites();ss++){
        for (int mu=0; mu<Nd; mu++) {
            for (int kk=0; kk<Np; kk++) {
                for (int nn=0; nn<Nc; nn++) {
                    (U._odata[ss])(mu)()(kk)(nn,nn) = 1.0;
                }
            }
        }
    }
    // perturbative fermion fields, matrices in colour space, all matrices set to the identity (-> TBC are irrelevant)
    QCDpt::LatticeFermion psi(&Grid),psi2(&Grid),psicmp(&Grid);
    psi=zero;
    parallel_for(int ss=0;ss<Grid.oSites();ss++){
        for (int alpha=0; alpha<Ns; alpha++) {
            for (int kk=0; kk<Np; kk++) {
                for (int nn=0; nn<Nc; nn++) {
                    (psi._odata[ss])()(alpha)(kk)(nn,nn) = 1.0;
                }
            }
        }
    }
    // perturbative Wilson Dirac operator
    WilsonFermion<QCDpt::PWilsonSmellImpl<vComplex,FundamentalRepresentation,CoeffReal>> Dw(U,Grid,RBGrid,mass,Params);
    
    
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
            
            // apply Mstd
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
                    (psi2._odata[ss])()(alpha)(kk) -= (psicmp._odata[ss])()(alpha)(ord);
                }
                (psi2._odata[ss])()(alpha)(kk) += 4.*(double)kk*(psi._odata[ss])()(alpha)(kk);
            }
        }
    }
    
    PComplex normvalue = Pnorm2(psi2);
    
    bool ok=true;
    for (int k=0; k<Np; k++)
    if (norm(normvalue(k))>1e-12)
    ok=false;
    
    cout << GridLogMessage << "perturbative Wilson operator is consistent with standard Wilson operator : ";
    if(ok) std::cout << "[OK] " << endl;
    else cout << "[fail]" << endl << normvalue <<endl;
    
    Grid_finalize();
    return EXIT_SUCCESS;
    
}
