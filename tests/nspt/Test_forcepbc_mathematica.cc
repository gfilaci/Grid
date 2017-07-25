/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: ./tests/Test_forcepbc_mathematica.cc

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
    
    if(Nc!=2||Nd!=4||Np!=3){
        cout<<"This test has to be compiled with Nc=2, Nd=4, Np=3 and executed with L=2."<<endl;
        return EXIT_FAILURE;
    }
    
    Grid_init(&argc,&argv);
    
    std::vector<int> latt_size   = GridDefaultLatt();
    std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
    std::vector<int> mpi_layout  = GridDefaultMpi();
    
    GridCartesian               Grid(latt_size,simd_layout,mpi_layout);
    
    int threads = GridThread::GetThreads();
    std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;
    
    GridParallelRNG          pRNG(&Grid);
    pRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));
    
    
    double beta = 1.;
    WilsonGaugeAction<PeriodicGaugeImpl<GimplTypes_ptR>> Action(beta);
    
    QCDpt::LatticeGaugeField U(&Grid);
    QCDpt::LatticeGaugeField F(&Grid);
    
    zeroit(U);
    U = ProjectOnGroup(U);
    
    
    /////////////////////
    //// TEST OF THE FORCE
    //// to be tested with Nc=2, Nd=4, Np=3, L=2, tau=0.0025*2*Nc
    /////////////////////
    
    double tau = 0.0025*2*2;
    
    iMatrix<Complex,Nc> A1,A2,A3,A4,A5,A6,A7;
    Complex im(0,1);
    
    A1(0,0) = 0. + 0.*im;
    A1(0,1) = 0. + 1.*im;
    A1(1,0) = -1. + 0.*im;
    A1(1,1) = 0. + 0.*im;
    
    A2(0,0) = 1. + 0.*im;
    A2(0,1) = 1. + 1.*im;
    A2(1,0) = -1. + 0.*im;
    A2(1,1) = 0. + 0.*im;
    
    A3(0,0) = 0. + 0.*im;
    A3(0,1) = 0. + 0.*im;
    A3(1,0) = -1. + 0.*im;
    A3(1,1) = 1. + 0.*im;
    
    A4(0,0) = 0. + 1.*im;
    A4(0,1) = -1. + 0.*im;
    A4(1,0) = -1. + 0.*im;
    A4(1,1) = 0. + 1.*im;
    
    A5(0,0) = 0.5 + 0.*im;
    A5(0,1) = 0. + -7.*im;
    A5(1,0) = -1. + 1.*im;
    A5(1,1) = 0. + 0.*im;
    
    A6(0,0) = 0.5 + 1.*im;
    A6(0,1) = 0. + 3.*im;
    A6(1,0) = -1. + 0.*im;
    A6(1,1) = 1. + 0.*im;
    
    A7(0,0) = 0. + -1.*im;
    A7(0,1) = 4. + -1.*im;
    A7(1,0) = -1. + 0.*im;
    A7(1,1) = 2. + 0.*im;
    
//    cout<<A1<<endl;
//    cout<<A2<<endl;
//    cout<<A3<<endl;
//    cout<<A4<<endl;
//    cout<<A5<<endl;
//    cout<<A6<<endl;
//    cout<<A7<<endl;
    
    LorentzPertColourMatrix tmpmat;
    for (int x0=0; x0<2; x0++) {
        for (int x1=0; x1<2; x1++) {
            for (int x2=0; x2<2; x2++) {
                for (int x3=0; x3<2; x3++) {
                    vector<double> dcoord({(double)x0,(double)x1,(double)x2,(double)x3});
                    zeroit(tmpmat);
                    tmpmat = ProjectOnGroup(tmpmat);
                    for (int mu=0; mu<4; mu++) {
                        for (int k=0; k<2; k++) {
                            tmpmat(mu)._internal(k+1) = dcoord[0]*A1 + dcoord[1]*A2 + dcoord[2]*A3 + dcoord[3]*A4 - (double)((k+1)*mu)*A5 + dcoord[2]*(double)(k+1)*A6 - dcoord[0]*(double)mu*A7 ;
                        }
                    }
                    pokeSite(tmpmat,U,std::vector<int>({x0,x1,x2,x3}));
                }
            }
        }
    }

    Action.deriv(U,F);
    F = -tau * F;
    U = Exponentiate(F) * U;
    
    cout<<U<<endl;
    
    Grid_finalize();
    return EXIT_SUCCESS;
    
}

/////////////////////
// Compare with results from this Mathematica example code
/////////////////////

//pTimes[{A1_, A2_}, {B1_, B2_}] := {A1 + B1, A2 + B2 + A1.B1};
//Staple[U1_, U2_, U3_, U4_] := pTimes[pTimes[pTimes[U1, U2], U3], U4];
//reH[A_] := 
//  1/2 (A - ConjugateTranspose[A] - 
//     1/2 Tr[A - ConjugateTranspose[A]] IdentityMatrix[2]);
//preH[{A1_, A2_}] := {reH[A1], reH[A2]};
//dag[{A1_, A2_}] := {ConjugateTranspose[A1], ConjugateTranspose[A2]};
//pExp[{A1_, A2_}] := {A1, A2 + 1/2 A1.A1};
//CompleteStaple[U0_, Uxf1_, Uxf2_, Uxf3_, Uxb1_, Uxb2_, Uxb3_, Uyf1_, 
//   Uyf2_, Uyf3_, Uyb1_, Uyb2_, Uyb3_, Uzf1_, Uzf2_, Uzf3_, Uzb1_, 
//   Uzb2_, Uzb3_] := 
//  Staple[U0, Uxf1, Uxf2, Uxf3] + Staple[U0, Uxb1, Uxb2, Uxb3] + 
//   Staple[U0, Uyf1, Uyf2, Uyf3] + Staple[U0, Uyb1, Uyb2, Uyb3] + 
//   Staple[U0, Uzf1, Uzf2, Uzf3] + Staple[U0, Uzb1, Uzb2, Uzb3];
//EvolveLink[taug_, U0_, Uxf1_, Uxf2_, Uxf3_, Uxb1_, Uxb2_, Uxb3_, 
//   Uyf1_, Uyf2_, Uyf3_, Uyb1_, Uyb2_, Uyb3_, Uzf1_, Uzf2_, Uzf3_, 
//   Uzb1_, Uzb2_, Uzb3_] := 
//  pTimes[pExp[
//    taug*preH[
//      CompleteStaple[U0, Uxf1, Uxf2, Uxf3, Uxb1, Uxb2, Uxb3, Uyf1, 
//       Uyf2, Uyf3, Uyb1, Uyb2, Uyb3, Uzf1, Uzf2, Uzf3, Uzb1, Uzb2, 
//       Uzb3]]], U0];
//PrintLink[{A1_, A2_}] := {A1 // MatrixForm, A2 // MatrixForm};
//A1 = {{0, I}, {-1, 0}};
//A1 // MatrixForm
//A2 = {{1, 1 + I}, {-1, 0}};
//A2 // MatrixForm
//A3 = {{0, 0}, {-1, 1}};
//A3 // MatrixForm
//A4 = {{I, -1}, {-1, I}};
//A4 // MatrixForm
//A5 = {{0.5, -7 I}, {-1 + I, 0}};
//A5 // MatrixForm
//A6 = {{0.5 + I, 3 I}, {-1, 1}};
//A6 // MatrixForm
//A7 = {{-I, 4 - I}, {-1, 2}};
//A7 // MatrixForm
//myconfig[x0_, x1_, x2_, x3_, mu_, k_] := 
//  x0*A1 + x1*A2 + x2*A3 + x3*A4 - k*mu*A5 + x2*k*A6 - x0*mu*A7;
//link[x0_, x1_, x2_, x3_, mu_] := {myconfig[x0, x1, x2, x3, mu, 1], 
//   myconfig[x0, x1, x2, x3, mu, 2]};
//   
//   PrintLink[link[1, 0, 1, 0, 3]]
//   
//   
//   (* update link (0,0,0,0) mu=0 (site=0) *)
//   
//PrintLink[EvolveLink[-0.0025, link[0, 0, 0, 0, 0],
//  link[1, 0, 0, 0, 1], dag[link[0, 1, 0, 0, 0]], 
//  dag[link[0, 0, 0, 0, 1]],
//  link[1, 0, 0, 0, 2], dag[link[0, 0, 1, 0, 0]], 
//  dag[link[0, 0, 0, 0, 2]],
//  link[1, 0, 0, 0, 3], dag[link[0, 0, 0, 1, 0]], 
//  dag[link[0, 0, 0, 0, 3]],
//  dag[link[1, 1, 0, 0, 1]], dag[link[0, 1, 0, 0, 0]], 
//  link[0, 1, 0, 0, 1],
//  dag[link[1, 0, 1, 0, 2]], dag[link[0, 0, 1, 0, 0]], 
//  link[0, 0, 1, 0, 2],
//  dag[link[1, 0, 0, 1, 3]], dag[link[0, 0, 0, 1, 0]], 
//  link[0, 0, 0, 1, 3]]]
//
// (* update link (0,1,0,0) mu=1 (site=4) *)
//
//PrintLink[EvolveLink[-0.0025, link[0, 1, 0, 0, 1],
//  link[0, 0, 0, 0, 0], dag[link[1, 1, 0, 0, 1]], 
//  dag[link[0, 1, 0, 0, 0]],
//  link[0, 0, 0, 0, 2], dag[link[0, 1, 1, 0, 1]], 
//  dag[link[0, 1, 0, 0, 2]],
//  link[0, 0, 0, 0, 3], dag[link[0, 1, 0, 1, 1]], 
//  dag[link[0, 1, 0, 0, 3]],
//  dag[link[1, 0, 0, 0, 0]], dag[link[1, 1, 0, 0, 1]], 
//  link[1, 1, 0, 0, 0],
//  dag[link[0, 0, 1, 0, 2]], dag[link[0, 1, 1, 0, 1]], 
//  link[0, 1, 1, 0, 2],
//  dag[link[0, 0, 0, 1, 3]], dag[link[0, 1, 0, 1, 1]], 
//  link[0, 1, 0, 1, 3]]]
//
// (* update link (1,1,1,1) mu=3 (site=15) *)
//
//PrintLink[EvolveLink[-0.0025, link[1, 1, 1, 1, 3],
//  link[1, 1, 1, 0, 0], dag[link[0, 1, 1, 1, 3]], 
//  dag[link[1, 1, 1, 1, 0]],
//  link[1, 1, 1, 0, 1], dag[link[1, 0, 1, 1, 3]], 
//  dag[link[1, 1, 1, 1, 1]],
//  link[1, 1, 1, 0, 2], dag[link[1, 1, 0, 1, 3]], 
//  dag[link[1, 1, 1, 1, 2]],
//  dag[link[0, 1, 1, 0, 0]], dag[link[0, 1, 1, 1, 3]], 
//  link[0, 1, 1, 1, 0],
//  dag[link[1, 0, 1, 0, 1]], dag[link[1, 0, 1, 1, 3]], 
//  link[1, 0, 1, 1, 1],
//  dag[link[1, 1, 0, 0, 2]], dag[link[1, 1, 0, 1, 3]], 
//  link[1, 1, 0, 1, 2]]]
