#include <Grid/Grid.h>

using namespace std;
using namespace Grid;
using namespace QCD;
using namespace QCDpt;

int main(int argc, char *argv[]) {
    
    Grid_init(&argc,&argv);
    
    std::vector<int> latt_size   = GridDefaultLatt();
    std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
    std::vector<int> mpi_layout  = GridDefaultMpi();
    
    GridCartesian Grid(latt_size,simd_layout,mpi_layout);
    
    GridParallelRNG pRNG(&Grid);
    pRNG.SeedFixedIntegers({5,12,3,5});
    
    Lattice<iScalar<vRealD>> latD(&Grid);
    //    Lattice<iScalar<__float128>> lat(&Grid);
    iPert<iMatrix<__complex128,3>,10> test;
    /*
    random(pRNG,latD);
    random128(pRNG,lat);

    cout << latD << endl;
    cout << lat << endl;
    */

    Grid_finalize();
    return EXIT_SUCCESS;
}
