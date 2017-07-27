/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: ./tests/Test_namespace.cc

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

static constexpr double tolerance = 1.0e-12;

void test(const double &a, const double &b)
{
  if (a - b < tolerance)
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
void print_test(const string &text, const Expr &a, const Expr &b)
{
    std::cout << GridLogMessage << text << " : ";
    test(a,b);
    std::cout << std::endl;
}

int main(int argc, char *argv[]) {
    
    Grid_init(&argc,&argv);
    
    std::vector<int> latt_size   = GridDefaultLatt();
    std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
    std::vector<int> mpi_layout  = GridDefaultMpi();
    
    GridCartesian               Grid(latt_size,simd_layout,mpi_layout);
    
    int threads = GridThread::GetThreads();
    std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;
    
    GridParallelRNG          pRNG(&Grid);
    pRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));

/////////////////////////////////////////////////////////
/////////////////// TESTS
/////////////////////////////////////////////////////////
    
    std::cout << GridLogMessage << "======== Test peek perturbative and non-perturbative indices" << std::endl;

    QCD::LatticeGaugeField gfield(&Grid);
    QCDpt::LatticeGaugeField gfieldpt(&Grid);
    int counter;
    
    gaussian(pRNG,gfield);
    gaussian(pRNG,gfieldpt);
    
    counter = 0;
    for (int i=0; i<Np; i++) if (norm2(PeekIndex<2>(gfieldpt,i))-norm2(peekPert(gfieldpt,i))>tolerance) counter++;
    print_test("peek pert index                           ",0.,(double)counter);
    
    counter = 0;
    for (int i=0; i<Nc; i++) for (int j=0; j<Nc; j++) if(norm2(PeekIndex<3>(gfieldpt,i,j))-norm2(peekColour(gfieldpt,i,j))>tolerance) counter++;
    print_test("peek colour index (perturbative)          ",0.,(double)counter);
    
    counter = 0;
    for (int i=0; i<Nc; i++) for (int j=0; j<Nc; j++) if(norm2(PeekIndex<2>(gfield,i,j))-norm2(peekColour(gfield,i,j))>tolerance) counter++;
    print_test("peek colour index (non-perturbative)      ",0.,(double)counter);
    
    
    std::cout << GridLogMessage << "======== Test perturbative trait" << std::endl;
    
    GimplTypes_ptR traitpt;
    GimplTypesR traitnonpt;
    
    print_test("GimplTypes_ptR is perturbative            ",1.,(double)isPerturbative<GimplTypes_ptR::Field>::value);
    print_test("GimplTypesR is not perturbative           ",0.,(double)isPerturbative<GimplTypesR::Field>::value);
    
    Grid_finalize();
    return EXIT_SUCCESS;
    
}
