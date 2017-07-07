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

static constexpr double tolerance = 1.0e-12;

const int vsz = 10;
const int msz = 9;
const int psz = 8;


template <typename Expr>
void test(const Expr &a, const Expr &b)
{
  if (norm2(a - b) < tolerance)
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
void test(const Expr &a)
{
  if (norm2(a) < tolerance)
  {
    std::cout << "[OK] ";
  }
  else
  {
    std::cout << "[fail]" << std::endl;
    std::cout << GridLogError << "a= " << a << std::endl;
    std::cout << GridLogError << "does not vanish" << std::endl;
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

template <typename Expr>
void print_test(const string &text, const Expr &a)
{
    std::cout << GridLogMessage << text << " : ";
    test(a);
    std::cout << std::endl;
}

int main(int argc, char *argv[]) {
    
    Grid_init(&argc,&argv);
    
    // initialise random number generator
    GridSerialRNG sRNG;
    sRNG.SeedFixedIntegers(std::vector<int>({(int)time(0)}));
    
    
    
    // declare and generate random objects
    iVector<iPert<iMatrix<Complex,msz>,psz>,vsz> S,P,T,Tman;
    iScalar<iPert<iMatrix<Complex,msz>,psz>> Q;
    iScalar<iScalar<iScalar<Complex>>> sc;
    iScalar<iScalar<iVector<Complex,msz>>> vec,Tvec,Tvecman;
    iScalar<iScalar<iMatrix<Complex,msz>>> mat;
    
    random(sRNG, S);
    random(sRNG, P);
    random(sRNG, Q);
    random(sRNG, sc);
    random(sRNG, vec);
    random(sRNG, mat);
    
    
    // BEGIN TESTS
    
    std::cout << GridLogMessage << "======== Test sum" << std::endl;
    
    zeroit(Tman);
    loopv(i) loopp(j)
    Tman(i)(j) = S(i)(j) + P(i)(j);
    T = S + P;
    print_test("sum         (internal)  ",T,Tman);
    
    loopp(j)
    T(j) = S(j) + P(j);
    print_test("sum         (overloaded)",T,Tman);
    
    
    
    std::cout << GridLogMessage << "======== Test subtraction" << std::endl;
    
    zeroit(Tman);
    loopv(i) loopp(j)
    Tman(i)(j) = S(i)(j) - P(i)(j);
    T = S - P;
    print_test("sub         (internal)  ",T,Tman);
    
    loopp(j)
    T(j) = S(j) - P(j);
    print_test("sub         (overloaded)",T,Tman);
    
    std::cout << GridLogMessage << "======== Test multiplication table" << std::endl;
    
    zeroit(Tvecman);
    loopm(i) loopm(j)
    Tvecman._internal._internal(i) += vec._internal._internal(j) * mat._internal._internal(j,i);
    Tvec = vec * mat;
    print_test("vec x mat   (internal)  ",Tvec,Tvecman);
    Tvec._internal._internal = vec._internal._internal * mat._internal._internal;
    print_test("vec x mat   (overloaded)",Tvec,Tvecman);
    
    
    zeroit(Tman);
    loopv(i) loopp(j)
    Tman(i)(j) = sc._internal._internal * S(i)(j);
    T = sc * S;
    print_test("scal x pert (internal)  ",T,Tman);
    loopv(i)
    T(i) = sc._internal * S(i);
    print_test("scal x pert (overloaded)",T,Tman);
    
    
    zeroit(Tman);
    loopv(i) loopp(j)
    Tman(i)(j) = S(i)(j) * sc._internal._internal;
    T = S * sc;
    print_test("pert x scal (internal)  ",T,Tman);
    loopv(i)
    T(i) = S(i) * sc._internal;
    print_test("pert x scal (overloaded)",T,Tman);
    
    
    zeroit(Tman);
    loopv(i)
    for (int c1=0; c1<psz; c1++) {
        for (int c2=0; c2<=c1; c2++) {
            Tman(i)(c1) += S(i)(c2) * (Q._internal)(c1-c2);
        }
    }
    T = S * Q;
    print_test("pert x pert (internal)  ",T,Tman);
    loopv(i)
    T(i) = S(i) * Q._internal;
    print_test("pert x pert (overloaded)",T,Tman);
    
    
    
    std::cout << GridLogMessage << "======== Test division by scalar" << std::endl;
    
    zeroit(Tman);
    loopv(i) loopp(j)
    Tman(i)(j) = S(i)(j) / sc._internal._internal;
    T = S / sc;
    print_test("division    (internal)  ",T,Tman);
    loopv(i)
    T(i) = S(i) / sc._internal;
    print_test("division    (overloaded)",T,Tman);
    
    
    
    std::cout << GridLogMessage << "======== Test multiplication by fundamental type" << std::endl;
    
    iPert<double,psz> Pdouble,Tdouble,Tdoubleman;
    iPert<iScalar<double>,psz> Zdoubleres;
    random(sRNG,Pdouble);
    double mydouble;
    random(sRNG,mydouble);
    
    loopp(i)
    Tdoubleman(i) = Pdouble(i) * mydouble;
    Tdouble = Pdouble * mydouble;
    loopp(i)
    Zdoubleres(i) = Tdouble(i) - Tdoubleman(i);
    print_test("pert x double           ",Zdoubleres);
    Tdouble = mydouble * Pdouble;
    loopp(i)
    Zdoubleres(i) = Tdouble(i) - Tdoubleman(i);
    print_test("double x pert           ",Zdoubleres);
    
    
    iPert<ComplexD,psz> Pcomplex,Tcomplex,Tcomplexman;
    iPert<iScalar<ComplexD>,psz> Zcomplexres;
    random(sRNG,Pcomplex);
    ComplexD mycomplex;
    random(sRNG,mycomplex);
    
    loopp(i)
    Tcomplexman(i) = Pcomplex(i) * mycomplex;
    Tcomplex = Pcomplex * mycomplex;
    loopp(i)
    Zcomplexres(i) = Tcomplex(i) - Tcomplexman(i);
    print_test("pert x Complex          ",Zcomplexres);
    Tcomplex = mycomplex * Pcomplex;
    loopp(i)
    Zcomplexres(i) = Tcomplex(i) - Tcomplexman(i);
    print_test("Complex x pert          ",Zcomplexres);
    
    
    Grid_finalize();
    return EXIT_SUCCESS;
}
