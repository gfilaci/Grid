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

const int vsz = 1;
const int msz = 2;
const int psz = 1;

const Complex im(0,1);

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

void test(const int &a, const int &b)
{
  if (a == b)
  {
    std::cout << "[OK] ";
  }
  else
  {
    std::cout << "[fail]" << std::endl;
    std::cout << GridLogError << "a= " << a << std::endl;
    std::cout << GridLogError << "is different from " << std::endl;
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
    iVector<iPert<iMatrix<double,msz>,psz>,vsz> Treim,Treimman;
    iVector<iScalar<iMatrix<Complex,msz>>,vsz> Ppeek,Ppeek2;
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
    
    loopv(i) loopp(j)
    Tman(i)(j) = S(i)(j) + P(i)(j);
    T = S + P;
    print_test("sum         (internal)              ",T,Tman);
    
    loopp(j)
    T(j) = S(j) + P(j);
    print_test("sum         (overloaded)            ",T,Tman);
    
    
    
    std::cout << GridLogMessage << "======== Test subtraction" << std::endl;
    
    loopv(i) loopp(j)
    Tman(i)(j) = S(i)(j) - P(i)(j);
    T = S - P;
    print_test("sub         (internal)              ",T,Tman);
    
    loopp(j)
    T(j) = S(j) - P(j);
    print_test("sub         (overloaded)            ",T,Tman);
    
    std::cout << GridLogMessage << "======== Test multiplication table" << std::endl;
    
    zeroit(Tvecman);
    loopm(i) loopm(j)
    Tvecman._internal._internal(i) += vec._internal._internal(j) * mat._internal._internal(j,i);
    Tvec = vec * mat;
    print_test("vec x mat   (internal)              ",Tvec,Tvecman);
    Tvec._internal._internal = vec._internal._internal * mat._internal._internal;
    print_test("vec x mat   (overloaded)            ",Tvec,Tvecman);
    
    
    loopv(i) loopp(j)
    Tman(i)(j) = sc._internal._internal * S(i)(j);
    T = sc * S;
    print_test("scal x pert (internal)              ",T,Tman);
    loopv(i)
    T(i) = sc._internal * S(i);
    print_test("scal x pert (overloaded)            ",T,Tman);
    
    
    loopv(i) loopp(j)
    Tman(i)(j) = S(i)(j) * sc._internal._internal;
    T = S * sc;
    print_test("pert x scal (internal)              ",T,Tman);
    loopv(i)
    T(i) = S(i) * sc._internal;
    print_test("pert x scal (overloaded)            ",T,Tman);
    
    
    zeroit(Tman);
    loopv(i)
    for (int c1=0; c1<psz; c1++) {
        for (int c2=0; c2<=c1; c2++) {
            Tman(i)(c1) += S(i)(c2) * (Q._internal)(c1-c2);
        }
    }
    T = S * Q;
    print_test("pert x pert (internal)              ",T,Tman);
    loopv(i)
    T(i) = S(i) * Q._internal;
    print_test("pert x pert (overloaded)            ",T,Tman);
    
    
    
    std::cout << GridLogMessage << "======== Test division by scalar" << std::endl;
    
    
    loopv(i) loopp(j)
    Tman(i)(j) = S(i)(j) / sc._internal._internal;
    T = S / sc;
    print_test("division    (internal)              ",T,Tman);
    loopv(i)
    T(i) = S(i) / sc._internal;
    print_test("division    (overloaded)            ",T,Tman);
    
    
    
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
    print_test("pert x double                       ",Zdoubleres);
    Tdouble = mydouble * Pdouble;
    loopp(i)
    Zdoubleres(i) = Tdouble(i) - Tdoubleman(i);
    print_test("double x pert                       ",Zdoubleres);
    
    
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
    print_test("pert x complex                      ",Zcomplexres);
    Tcomplex = mycomplex * Pcomplex;
    loopp(i)
    Zcomplexres(i) = Tcomplex(i) - Tcomplexman(i);
    print_test("complex x pert                      ",Zcomplexres);
    
    
    
    std::cout << GridLogMessage << "======== Test trace" << std::endl;
    iVector<iPert<iScalar<Complex>,psz>,vsz> Ttrace, Ttraceman;
    loopv(i) loopp(j)
    Ttraceman(i)(j) = trace(P(i)(j));
    Ttrace = trace(P);
    print_test("trace                               ",Ttrace,Ttraceman);
    
    
    
    std::cout << GridLogMessage << "======== Test reality" << std::endl;
    
    Tman = im * P;
    timesI(T,P);
    print_test("times I  (2 arguments)              ",T,Tman);
    T = timesI(T);
    print_test("times I  (1 argument)               ",T,-P);
    
    Tman = -im * P;
    timesMinusI(T,P);
    print_test("times -I (2 arguments)              ",T,Tman);
    T = timesMinusI(T);
    print_test("times -I (1 argument)               ",T,-P);
    
    loopv(i)loopp(j)
    Tman(i)(j) = conjugate(P(i)(j));
    T = conjugate(P);
    print_test("complex conjugation                 ",T,Tman);
    
    loopv(i)loopp(j)
    Tman(i)(j) = adj(P(i)(j));
    T = adj(P);
    print_test("adjoint                             ",T,Tman);
    
    loopv(i)loopp(j)
    Treimman(i)(j) = real(P(i)(j));
    Treim = real(P);
    print_test("real part                           ",Treim,Treimman);
    
    loopv(i)loopp(j)
    Treimman(i)(j) = imag(P(i)(j));
    Treim = imag(P);
    print_test("imaginary part                      ",Treim,Treimman);
    
    
    
    std::cout << GridLogMessage << "======== Test projections" << std::endl;
    
    loopv(i)loopp(j)
    Tman(i)(j) = Ta(P(i)(j));
    T = Ta(P);
    print_test("projection on algebra   ",T,Tman);

// check projection on group not with a random matrix...
//    loopv(i)loopp(j)
//    Tman(i)(j) = ProjectOnGroup(P(i)(j));
//    T = ProjectOnGroup(P);
//    print_test("projection on group     ",T,Tman);
    
    
    std::cout << GridLogMessage << "======== Test index" << std::endl;
    
    print_test("index rank                          ",indexRank<1,decltype(P)>(),psz);
    print_test("iMatrix is iMatrix                  ",isMatrix<2,decltype(P)>(),1);
    print_test("iPert is not iMatrix                ",isMatrix<1,decltype(P)>(),0);
    print_test("iPert is iPert                      ",isPert<1,decltype(P)>(),1);
    print_test("iScalar is not iPert                ",isPert<0,decltype(sc)>(),0);
    print_test("trace index                         ",traceIndex<2>(P),trace(P));
    // reminder:
    // iVector<iPert<iMatrix<Complex,msz>,psz>,vsz> P;
    // iVector<iScalar<iMatrix<Complex,msz>>,vsz> Ppeek;
    // fill Ppeek with the 2nd perturbative order of P
    loopv(i)
    Ppeek(i)._internal = P(i)(2);
    print_test("peek index                          ",peekIndex<1>(P,2),Ppeek);
    // take Ppeek at level 1 and put it into the 0 component of T at level 1
    pokeIndex<1>(T,Ppeek,0);
    loopv(i)
    Ppeek2(i)._internal = T(i)(0);
    print_test("poke index                          ",Ppeek,Ppeek2);
    loopv(i) loopp(j)
    T(i)(j) = transpose(P(i)(j));
    print_test("transpose index (not working?)      ",transposeIndex<2>(P),conjugate(adj((T))));
    print_test("trace index                         ",traceIndex<2>(P),Ttraceman);
    


    Grid_finalize();
    return EXIT_SUCCESS;
}
