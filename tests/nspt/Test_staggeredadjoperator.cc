    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/nspt/Test_staggeredadjoperator.cc

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>
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
using namespace Grid::QCD;
using namespace Grid::QCD::QCDpt;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
  GridCartesian               Grid(latt_size,simd_layout,mpi_layout);
  GridRedBlackCartesian     RBGrid(&Grid);

  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;
  std::cout<<GridLogMessage << "Grid floating point word size is REALF"<< sizeof(RealF)<<std::endl;
  std::cout<<GridLogMessage << "Grid floating point word size is REALD"<< sizeof(RealD)<<std::endl;
  std::cout<<GridLogMessage << "Grid floating point word size is REAL"<< sizeof(Real)<<std::endl;

  GridParallelRNG          pRNG(&Grid);
  pRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));

  typedef typename StaggeredFermion<QCDpt::PStaggeredAdjointImplD>::FermionField FermionField;
	typedef typename StaggeredFermion<QCDpt::PStaggeredAdjointImplD>::SOimpl::ComplexField ComplexField;
  typename StaggeredFermion<QCDpt::PStaggeredAdjointImplD>::ImplParams params;
  params.boundary_phases = {-1.,1.,1.,1.};
	
	FermionField src   (&Grid); gaussian(pRNG,src);
  FermionField result(&Grid); result=zero;
  FermionField    ref(&Grid);    ref=zero;
  FermionField    tmp(&Grid);    tmp=zero;
  FermionField    err(&Grid);    tmp=zero;
  FermionField phi   (&Grid); random(pRNG,phi);
  FermionField chi   (&Grid); random(pRNG,chi);
	QCDpt::LatticeGaugeField Umu(&Grid); QCDpt::PertRandom(pRNG,Umu);
	std::vector<QCDpt::LatticePertColourMatrix> U(4,&Grid);

	QCDpt::TwistedGaugeImpl<QCDpt::GimplTypes_ptR> TwistedBC;
	
  double volume=1;
  for(int mu=0;mu<Nd;mu++){
    volume=volume*latt_size[mu];
  }  

  // Only one non-zero (y)
  for(int mu=0;mu<Nd;mu++){
	  U[mu] = PeekIndex<LorentzIndex>(Umu,mu);
  // Debug force unit
//    U[mu] = 1.0;
//    PokeIndex<LorentzIndex>(Umu,U[mu],mu);
  }

  ref = zero;

  RealD mass=1.;
  
  { // Simple staggered implementation
    ref = zero;

    Lattice<iScalar<vInteger> > coor(&Grid);

    Lattice<iScalar<vInteger> > x(&Grid); LatticeCoordinate(x,0);
    Lattice<iScalar<vInteger> > y(&Grid); LatticeCoordinate(y,1);
    Lattice<iScalar<vInteger> > z(&Grid); LatticeCoordinate(z,2);
    Lattice<iScalar<vInteger> > t(&Grid); LatticeCoordinate(t,3);

    Lattice<iScalar<vInteger> > lin_z(&Grid); lin_z=x+y;
    Lattice<iScalar<vInteger> > lin_t(&Grid); lin_t=x+y+z;

    for(int mu=0;mu<Nd;mu++){

      // Staggered Phase.
      ComplexField phases(&Grid);	phases=1.0;
  
      if ( mu == 1 ) phases = where( mod(x    ,2)==(Integer)0, phases,-phases);
      if ( mu == 2 ) phases = where( mod(lin_z,2)==(Integer)0, phases,-phases);
      if ( mu == 3 ) phases = where( mod(lin_t,2)==(Integer)0, phases,-phases);
		
	  int Lmu = Grid.GlobalDimensions()[mu] - 1;
	  LatticeCoordinate(coor, mu);
		
	  auto pha = params.boundary_phases[mu];
	  ComplexD phase( real(pha),imag(pha) );
		
	  tmp = Cshift(src,mu,1);
	  if(istwisted(mu)) tmp = where(coor==Lmu, phase*TwistedBC.twist.forward(tmp,mu), tmp);
	  else tmp = where(coor==Lmu, phase*tmp, tmp);
      ref = ref +0.5*phases*U[mu]*tmp*adj(U[mu]);
		
	  tmp = Cshift(src,mu,-1);
	  if(istwisted(mu)) tmp = where(coor==0, conjugate(phase)*TwistedBC.twist.backward(tmp,mu), tmp);
	  else tmp = where(coor==0, conjugate(phase)*tmp, tmp);
	  U[mu] = Cshift(U[mu],mu,-1);
	  if(istwisted(mu)) U[mu] = where(coor==0, TwistedBC.twist.backward(U[mu],mu), U[mu]);
	  ref = ref -0.5*phases*adj(U[mu])*tmp*U[mu];
		
    }
  }

  StaggeredFermion<QCDpt::PStaggeredAdjointImplD> Ds(Umu,Grid,RBGrid,mass,params);
	
  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;
  std::cout<<GridLogMessage<<"= Testing Dhop against cshift implementation         "<<std::endl;
  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;

  std::cout<<GridLogMessage << "Calling Ds"<<std::endl;
  int ncall=1000;
  double t0=usecond();
  for(int i=0;i<ncall;i++){
    Ds.Dhop(src,result,0);
  }
  double t1=usecond();
  double t2;
  double flops=(16*(3*(6+8+8)) + 15*3*2)*volume*ncall; // == 66*16 +  == 1146

  std::cout<<GridLogMessage << "Called Ds"<<std::endl;
  std::cout<<GridLogMessage << "norm result "<< Pnorm2(result)<<std::endl;
  std::cout<<GridLogMessage << "norm ref    "<< Pnorm2(ref)<<std::endl;
  std::cout<<GridLogMessage << "mflop/s =   "<< flops/(t1-t0)<<std::endl;

  err = ref-result;
  std::cout<<GridLogMessage << "norm diff   "<< Pnorm2(err)<<std::endl;
	
	// Only one non-zero (y)
	for(int mu=0;mu<Nd;mu++){
		U[mu] = PeekIndex<LorentzIndex>(Umu,mu);
		// Debug force unit
		//    U[mu] = 1.0;
		//    PokeIndex<LorentzIndex>(Umu,U[mu],mu);
	}
	
for(int muu=0;muu<Nd;muu++){
std::cout<<GridLogMessage<<"=========================================================="<<std::endl;
  std::cout<<GridLogMessage<<"= Testing Mdir"<<muu<<" against cshift implementation         "<<std::endl;
  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;

  { // Simple staggered implementation
    ref = zero;

    Lattice<iScalar<vInteger> > coor(&Grid);

    Lattice<iScalar<vInteger> > x(&Grid); LatticeCoordinate(x,0);
    Lattice<iScalar<vInteger> > y(&Grid); LatticeCoordinate(y,1);
    Lattice<iScalar<vInteger> > z(&Grid); LatticeCoordinate(z,2);
    Lattice<iScalar<vInteger> > t(&Grid); LatticeCoordinate(t,3);

    Lattice<iScalar<vInteger> > lin_z(&Grid); lin_z=x+y;
    Lattice<iScalar<vInteger> > lin_t(&Grid); lin_t=x+y+z;

      // Staggered Phase.
      ComplexField phases(&Grid);	phases=1.0;

      if ( muu == 1 ) phases = where( mod(x    ,2)==(Integer)0, phases,-phases);
      if ( muu == 2 ) phases = where( mod(lin_z,2)==(Integer)0, phases,-phases);
      if ( muu == 3 ) phases = where( mod(lin_t,2)==(Integer)0, phases,-phases);

	  int Lmu = Grid.GlobalDimensions()[muu] - 1;
	  LatticeCoordinate(coor, muu);
	  
	  auto pha = params.boundary_phases[muu];
	  ComplexD phase( real(pha),imag(pha) );
	  
	  tmp = Cshift(src,muu,1);
	  if(istwisted(muu)) tmp = where(coor==Lmu, phase*TwistedBC.twist.forward(tmp,muu), tmp);
	  else tmp = where(coor==Lmu, phase*tmp, tmp);
	  ref = ref +0.5*phases*U[muu]*tmp*adj(U[muu]);

    }

  std::cout<<GridLogMessage << "Calling Mdir"<<std::endl;
  double ncall=1000;
  double t0=usecond();
  for(int i=0;i<ncall;i++){
    Ds.Mdir(src,result,muu,DaggerNo);
  }
  double t1=usecond();
  double t2;
  double flops=(16*(3*(6+8+8)) + 15*3*2)*volume*ncall; // == 66*16 +  == 1146

  std::cout<<GridLogMessage << "Called Mdir"<<std::endl;
  std::cout<<GridLogMessage << "norm result "<< Pnorm2(result)<<std::endl;
  std::cout<<GridLogMessage << "norm ref    "<< Pnorm2(ref)<<std::endl;
  std::cout<<GridLogMessage << "mflop/s =   "<< flops/(t1-t0)<<std::endl;

  err = ref-result;
  std::cout<<GridLogMessage << "norm diff   "<< Pnorm2(err)<<std::endl;
}
                                      std::cout<<GridLogMessage<<"=========================================================="<<std::endl;
  std::cout<<GridLogMessage<<"= Testing that Deo + Doe = Dunprec "<<std::endl;
  std::cout<<GridLogMessage<<"=========================================================="<<std::endl;

  FermionField src_e   (&RBGrid);
  FermionField src_o   (&RBGrid);
  FermionField r_e   (&RBGrid);
  FermionField r_o   (&RBGrid);
  FermionField r_eo  (&Grid);
  pickCheckerboard(Even,src_e,src);
  pickCheckerboard(Odd,src_o,src);

  Ds.Meooe(src_e,r_o);  std::cout<<GridLogMessage<<"Applied Meo"<<std::endl;
  Ds.Meooe(src_o,r_e);  std::cout<<GridLogMessage<<"Applied Moe"<<std::endl;
  Ds.Dhop (src,ref,DaggerNo);

  setCheckerboard(r_eo,r_o);
  setCheckerboard(r_eo,r_e);

  err= ref - r_eo;
  std::cout<<GridLogMessage << "EO norm diff   "<< Pnorm2(err)<< " "<<Pnorm2(ref)<< " " << Pnorm2(r_eo) <<std::endl;
	
  Grid_finalize();
}
