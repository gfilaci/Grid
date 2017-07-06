/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: ./tests/Test_nspt.cc

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

#define loopv(n) for (int n=0; n<vsz; n++)
#define loopm(n) for (int n=0; n<msz; n++)
#define loopp(n) for (int n=0; n<psz; n++)

using namespace std;

static constexpr double tolerance = 1.0e-6;

const int vsz = 2;
const int msz = 2;
const int psz = 2;


template <typename Expr>
void test(const Expr &a, const Expr &b)
{
  if (norm2(a + (-b)) < tolerance)
  {
    std::cout << "[OK] ";
  }
  else
  {
    std::cout << "[fail]" << std::endl;
    std::cout << GridLogError << "a= " << a << std::endl;
    std::cout << GridLogError << "is different (tolerance= " << tolerance << ") from " << std::endl;
    std::cout << GridLogError << "b= " << b << std::endl;
    exit(EXIT_FAILURE);
  }
}

template <typename Expr>
void print_test(const string &text, const Expr &a, const Expr &b)
{
    std::cout << GridLogMessage << text << " : ";
    test(a,b);
    std::cout << std::endl;
}

int main(int argc, char *argv[]) {
    
    Grid_init(&argc,&argv);
    
    // initialise random number generator
    GridSerialRNG sRNG;
    sRNG.SeedFixedIntegers(std::vector<int>({(int)time(0)}));
    
    
    // declare and generate random objects
    iVector<iPert<iMatrix<Complex,msz>,psz>,vsz> S,P,T,Tman;
    iScalar<iScalar<iScalar<Complex>>> sc;
    iScalar<iScalar<iVector<Complex,msz>>> vec,Tvec,Tvecman;
    iScalar<iScalar<iMatrix<Complex,msz>>> mat;
    
    random(sRNG, S);
    random(sRNG, P);
    random(sRNG, sc);
    random(sRNG, vec);
    random(sRNG, mat);
    
    
    // BEGIN TESTS
    
    std::cout << GridLogMessage << "======== Test sum" << std::endl;
    
    zeroit(Tman);
    T = S + P;
    loopv(i) loopp(j)
    Tman(i)(j) = S(i)(j) + P(i)(j);
    print_test("sum (internal)",T,Tman);
    
    zeroit(Tman);
    T(0) = S(0) + P(0);
    loopp(j)
    Tman(0)(j) = S(0)(j) + P(0)(j);
    print_test("sum (overloaded)",T(0),Tman(0));
    
    
    std::cout << GridLogMessage << "======== Test multiplication table" << std::endl;
    
    zeroit(Tvecman);
    Tvec = vec * mat;
    loopm(i) loopm(j)
    Tvecman._internal._internal(i) += vec._internal._internal(j) * mat._internal._internal(j,i);
    print_test("vec x mat (internal)",Tvec,Tvecman);
    
    Tvec._internal._internal = vec._internal._internal * mat._internal._internal;
    print_test("vec x mat (internal)",Tvec,Tvecman);
    
//    T = sc * S;
//    loopv(i) loopp(j)
//    Tman(i)(j) = sc._internal._internal * S(i)(j);
//    
//    zeroit(Tman);
//    print_test("scal x pert (internal)",T,Tman);
//    zeroit(Tman);
//    print_test("scal x pert (overloaded)",T,Tman);
//    zeroit(Tman);
//    print_test("pert x scal (internal)",T,Tman);
//    zeroit(Tman);
//    print_test("pert x scal (overloaded)",T,Tman);
//    zeroit(Tman);
//    print_test("pert x pert (internal)",T,Tman);
//    zeroit(Tman);
//    print_test("pert x pert (overloaded)",T,Tman);
    
    
    
//    iVector<iMatrix<Complex,3>,4> aa,bb;
//    iScalar<iScalar<Complex> ss;
//    bb = ss * aa;
////    T = S * sc;
//    loopv(i) loopp(j)
//    Tman(i)(j) = S(i)(j) * sc._internal._internal;
//    print_test("pert x scal",T,Tman);
    
//    T = vec * S;
//    loopv(i) loopp(j)
////    Tman(i)(j) = S(i)(j) * vec._internal._internal;
//    print_test("pert x vec",T,Tman);
    
    
    iScalar<iVector<Complex,6>> ss;
    iVector<iVector<Complex,6>,10> bb;
    iVector<iMatrix<Complex,6>,10> aa;
//    iScalar<iVector<Complex,6>> ss;
//    iVector<iMatrix<Complex,6>,5> aa,bb;
    bb = aa*ss;
    
    
    
    Grid_finalize();
    return EXIT_SUCCESS;
}








//    iPert<double,2> v,w;
//    v(0)=1.;
//    v(1)=2.;
//    
//    w(0)=3.;
//    w(1)=4.;
//    cout<<v<<endl;
//    cout<<w<<endl;
//    cout<<v*w<<endl<<endl;
//
//
//    Complex im(0,1);
//    iPert<iMatrix<Complex,2>,3> S,P;
//    
//    S(0)(0,0) = 0;
//    S(0)(0,1) = 1;
//    S(0)(1,0) = -2;
//    S(0)(1,1) = -im;
//    
//    S(1)(0,0) = -23. + im;
//    S(1)(0,1) = 0;
//    S(1)(1,0) = -8;
//    S(1)(1,1) = 1;
//
//
//    P(1) = conjugate(transpose(S(1)));
//    
//    S(2) = transpose(S(0)*S(1));
//    
//    P(0) = S(2)*S(2);
//    
//    P(2)(0,0) = S(0)(0,0) + 2.;
//    P(2)(0,1) = S(0)(0,1);
//    P(2)(1,0) = S(0)(1,0);
//    P(2)(1,1) = S(0)(1,1) + 2.;
//    
//    cout<<S*P<<endl;




/*
compare in Mathematica with
A = {{0, 1}, {-2, -I}};
B = {{-23 + I, 0}, {-8, 1}};
CC = Transpose[A.B];
S0 = A;
S1 = B;
S2 = CC;
P0 = CC.CC;
P1 = ConjugateTranspose[B];
P2 = A + 2*IdentityMatrix[2];
S0.P0 // MatrixForm
S0.P1 + S1.P0 // MatrixForm
S0.P2 + S1.P1 + S2.P0 // MatrixForm
*/
