#include <Grid/Grid.h>

using namespace std;
using namespace Grid;
using namespace QCD;
using namespace QCDpt;

int main(int argc, char *argv[]) {
    
    Grid_init(&argc,&argv);
    
    std::vector<int> latt_size   = {4,4,4,4};
    std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
    std::vector<int> mpi_layout  = GridDefaultMpi();
    
    GridCartesian Grid(latt_size,simd_layout,mpi_layout);
    
    GridParallelRNG pRNG(&Grid);
    pRNG.SeedFixedIntegers({5,12,3,5});
    
    GridSerialRNG sRNG;
    sRNG.SeedFixedIntegers({6,17,7,65});
    
    typedef iMatrix<vComplexD,Nc> vtype;
    typedef typename vtype::scalar_object stype;
    //    typedef typename GridTypeMapper<vtype>::scalar_type stype;
      //typedef iMatrix<ComplexD,Nc>  stype;
    
    Lattice<iPert<vtype,Np>> Dlat(&Grid);
    iPert<vtype,Np> Dsite;
    PertRandom(pRNG,Dlat);
    cout << Dlat._odata[10] << endl;
    
    int Nsimd = sizeof(vtype::vector_type) / sizeof(vtype::scalar_type);
    std::vector<iPert<stype,Np>> buf(Nsimd);
    
    extract(Dlat._odata[10],buf);
    for(int i=0; i<Nsimd; i++){
      cout << "--- i=" << i << " ---" << endl;
      cout << buf[i] << endl;
    }
    merge(Dlat._odata[10],buf);
    
    iPert<iMatrix<__float128,Nc>,71> osite;

    Grid_finalize();
    return EXIT_SUCCESS;
}
