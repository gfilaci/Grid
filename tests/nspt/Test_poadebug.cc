/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: ./tests/nspt/Test_BinaryToScidac.cc

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

/*
 Use --configname, --rngname and --number to choose the file.
 Filenames must have the structure
 
 configname.number
 rngname.number
 
 The lattice size must be given with --grid
 and Grid must be compiled with the correct Np
 */
#include <iostream>
#include <quadmath.h>

using namespace std;

//typedef float      mytype;                                                                                                                                                 
//typedef double     mytype;                                                                                                                                                 
//typedef long double mytype;                                                                                                                                                
typedef __float128 mytype;

const int size = 71;
/*
ostream &operator<<(ostream& stream, const __float128 &x){
  stream << (double)x;
  return stream;
}
*/
template<class T>
T myabs(const T x) {return x>0 ? x : -x;}

void zero2(mytype *x)             {for(int i=0; i<size; i++) x[i] = (mytype)0;}
void setto(mytype *x, mytype *y) {for(int i=0; i<size; i++) x[i] = y[i];}
template<class T1, class T2> void settocast(T1 *x, T2 *y) {for(int i=0; i<size; i++) x[i] = (T1)y[i];}
void print(mytype *x){
  cout << "{";
  for(int i=0; i<size; i++){
    cout << x[i];
    if(i!=size-1) cout << ",";
  }
  cout << "}" << endl;
}
void compare(mytype *x, mytype *y){
  cout << "{";
  for(int i=0; i<size; i++){
    cout << (mytype)100*myabs((x[i]-y[i])/y[i]);
    if(i!=size-1) cout << ",";
  }
  cout << "}" << endl;
}

typedef void (*sumtype)(mytype&, mytype&, mytype&, mytype&, mytype&);

void compensatedLog(mytype *P, mytype *ret, sumtype sumfunc, bool finalcomp = false){
  mytype xin, sign = 1;
  mytype newtmp[size], tmp[size];
  setto(ret,P); setto(newtmp,P);
  ret[0] = (mytype)0;
  mytype e[size], eret[size];
  zero2(eret);
  for(int k=2; k<size; k++){
    zero2(tmp); zero2(e);
    for(int i=1; i<size+1-k; i++){
      for(int j=k-1; j<size-i; j++){
	mytype xin = newtmp[j]*P[i];
	sumfunc(tmp[i+j], e[i+j], tmp[i+j], xin, e[i+j]);
      }
    }
    if(finalcomp) for(int i=0; i<size; i++) tmp[i] += e[i];
    setto(newtmp,tmp);
    sign = -sign;
    for(int i=0; i<size; i++){
      xin = sign / (mytype)k * newtmp[i];
      sumfunc(ret[i], eret[i], ret[i], xin, eret[i]);
    }
  }
  if(finalcomp) for(int i=0; i<size; i++) ret[i] += eret[i];
}

void compensatedExp(mytype *P, mytype *ret, sumtype sumfunc, bool finalcomp = false){
  mytype xin;
  mytype newtmp[size], tmp[size];
  setto(ret,P); setto(newtmp,P);
  ret[0] = (mytype)1;
  mytype e[size], eret[size];
  zero2(eret);
  for(int k=2; k<size; k++){
    zero2(tmp); zero2(e);
    for(int i=1; i<size+1-k; i++){
      for(int j=k-1; j<size-i; j++){
	mytype xin = newtmp[j]*P[i];
	sumfunc(tmp[i+j], e[i+j], tmp[i+j], xin, e[i+j]);
      }
    }
    if(finalcomp) for(int i=0; i<size; i++) tmp[i] += e[i];
    for(int i=0; i<size; i++){
      xin = tmp[i] / (mytype)k;
      newtmp[i] = xin;
      sumfunc(ret[i], eret[i], ret[i], xin, eret[i]);
    }
  }
  if(finalcomp) for(int i=0; i<size; i++) ret[i] += eret[i];
}

void RecursiveSum(mytype &sumout, mytype &eout, mytype &sumin, mytype &xin, mytype &ein){
  sumout = sumin + xin;
}

void myLog      (mytype *P, mytype *ret) {compensatedLog(P,ret,RecursiveSum);}
void myExp      (mytype *P, mytype *ret) {compensatedExp(P,ret,RecursiveSum);}

int main(int argc, char *argv[]) {
    
    Grid_init(&argc,&argv);
    
    std::vector<int> latt_size = GridDefaultLatt();
    std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
    std::vector<int> mpi_layout  = GridDefaultMpi();
    
    GridCartesian Grid(latt_size,simd_layout,mpi_layout);
    
    GridParallelRNG pRNG(&Grid);
    GridSerialRNG   sRNG;
    
    typedef TwistedGimpl_ptR gimpl;
    typename gimpl::GaugeField U(&Grid), V(&Grid);
    string configname, rngname;
    int number;
    
    if( GridCmdOptionExists(argv,argv+argc,"--configname") ){
        configname = GridCmdOptionPayload(argv,argv+argc,"--configname");
    } else{
        std::cout << "ERROR: choose the name of the configuration with --configname" << std::endl;
        return EXIT_FAILURE;
    }
    if( GridCmdOptionExists(argv,argv+argc,"--rngname") ){
        rngname = GridCmdOptionPayload(argv,argv+argc,"--rngname");
    } else{
        std::cout << "ERROR: choose the name of the rng checkpoint with --rngname" << std::endl;
        return EXIT_FAILURE;
    }
    if( GridCmdOptionExists(argv,argv+argc,"--number") ){
        string arg = GridCmdOptionPayload(argv,argv+argc,"--number");
        std::stringstream ss(arg);
        ss>>number;
    } else{
        std::cout << "ERROR: choose the number of the configuration with --number" << std::endl;
        return EXIT_FAILURE;
    }
    
    CheckpointerParameters CPparams;
    
    CPparams.config_prefix = configname;
    CPparams.rng_prefix = rngname;
    CPparams.saveInterval = 1;
    CPparams.format = "IEEE64BIG";
        
    ScidacHmcCheckpointer<gimpl,ActionParameters> CP(CPparams);
    
    CP.CheckpointRestore(number,U,sRNG,pRNG);
    
    cout << std::setprecision(std::numeric_limits<Real>::digits10) << std::scientific;
    
    //iVector<iScalar<iPert<iMatrix<ComplexD,Nc>,Np>>,Nd> link;
    //iPert<iMatrix<Complex,Nc>,Np> mulink1, mulink2;
    //peekSite(link,U,std::vector<int>({12,3,8,6}));//this is the one that does not work
    //mulink1 = link(0)();
    /*
    mulink2 = Logarithm(mulink1);
    mulink2 = Exponentiate(mulink2);
    mulink2 -= mulink1;
    cout << "link that does not work after exp(log)" << endl << Pnorm2(mulink2).real()(Np-1) << endl;
    */
    //iPert<ComplexD,Np> norm1, norm2;
    //norm2 = Pnorm2(mulink1);
    //for(int k=0; k<Np; k++) norm1(k) = sqrt(norm2(k));
    //norm2 = Logarithm(norm1);
    //norm2 = Exponentiate(norm2);
    //norm2 -= norm1;
    //cout << "norm at the start" << endl << norm1 << endl;
    //cout << "norm that does not work after exp(log)" << endl << norm2 << endl;
    

    /*
    mytype sshort[10] = {1., 2.1238019542083837, 2.3355067272143577, 1.3374115080807851, 0.8167527513061994, 1.6390592537339208, 1.2146275469931898, 1.2797037856416371, 1.6626373260625067, 1.756478899439216};
    double slong[71] ={1., 2.1238019542083837, 2.3355067272143577, 1.3374115080807851, 0.8167527513061994, 1.6390592537339208, 1.2146275469931898, 1.2797037856416371, 1.6626373260625067, 1.756478899439216, 2.208241428931018, 3.070072392123372, 4.726777137894486, 6.891013435737664, 12.649584472045289, 26.471408453239214, 27.16856192583126,37.04331843602089, 82.39995786458687, 51.815917383481995, 87.11017642427103, 198.7735694522524, 234.64325454821883, 81.70071029911736, 764.6216818018614, 517.3758250636724, 1190.5106042401387, 2423.015391184666, 2891.3178116662853, 6748.547968884998, 11124.99341049452, 13043.835320180819, 29773.391044578875, 37346.72294335105, 59491.00458903986, 137429.66499669565, 193864.42573678892, 362366.3141099088, 521343.63383334083, 724154.6078615171, 1.3534599097372077*1.e6, 1.7599936603547235*1.e6, 3.9452511947024637*1.e6, 7.219380929373059*1.e6, 9.145858676793337*1.e6, 2.033867824161415*1.e7, 2.418259040333183*1.e7, 5.041092888666406*1.e7, 6.571864308286387*1.e7, 1.2405076935275543*1.e8, 1.7499257895408687*1.e8, 3.1891543321220237*1.e8, 4.3616573652778256*1.e8, 7.444590658960634*1.e8, 9.481150284681091*1.e8, 1.2430234753418758*1.e9, 2.4616499180782127*1.e9, 2.8744466329296713*1.e9, 5.1677726480783415*1.e9, 8.68511594779088*1.e9, 3.658632224631485*1.e10, 9.307413763756827*1.e10, 2.604672507971325*1.e11, 5.635722941028926*1.e11, 1.9664382463425808*1.e12, 3.711942389637628*1.e12, 1.4447869744449941*1.e13, 3.0709952323769516*1.e13, 6.563193380837487*1.e13, 2.2510151879142084*1.e14, 1.1284763461062465*1.e15};

    mytype a[size], b[size], c[size];
    if(size==10) settocast(a,sshort);
    else if(size==71) settocast(a,slong);
    else{
      cout << "ERROR" << endl;
      exit(EXIT_FAILURE);
    }

    cout << endl << "RECURSIVE SUMMATION, relative error (%)" << endl;
    myLog(a,b);
    cout << "starting a" << endl;
    print(a);
    cout << "after log" << endl;
    print(b);
    myExp(b,c);
    compare(c,a);







    Grid_finalize();
    return EXIT_SUCCESS;
    */
    U = ProjectOnGroup(U);
    V = ProjectOnGroup(U);
    auto plaq1 = WilsonLoops<gimpl>::avgPlaquette(U);
    auto plaq2 = WilsonLoops<gimpl>::avgPlaquette(V);
    V -= U;
    plaq2 -= plaq1;
    cout << "difference between plaquette after one POG and two POG" << endl << Pnorm2(plaq2) << endl;
    cout << "difference between fields after one POG and two POG" << endl << Pnorm2(V) << endl;

    cout << GridLogMessage << "start logarithm" << endl;    
    V = Logarithm(U);
    cout << GridLogMessage << "end logarithm"<< endl;
    cout << GridLogMessage << "start exponential"<< endl;
    V = Exponentiate(V);
    cout << GridLogMessage << "end exponential"<< endl;

    V -= U;
    cout << "checking exp(log) = id" << endl << Pnorm2(plaq2) << endl;
    
    Grid_finalize();
    return EXIT_SUCCESS;

    
    /************************************************/
    /*
    // test for fermion inversion //
    GridRedBlackCartesian     RBGrid(&Grid);
    //GridParallelRNG pRNG(&Grid);
    pRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));
    
    typedef PStaggeredAdjointImplR  fimpl;

    typedef typename fimpl::FermionField FermionField;
    FermionField psi2(&Grid);

    RealD mass = 0.;
    int Nf = 2;
    WilsonImplParams Params;
    Params.boundary_phases = {-1.,1.,1.,1.};

    StochasticAdjointStaggeredAction<fimpl> FA(&pRNG,&Grid,&RBGrid,mass,Params,Nf);

    // invM works on Npf orders                                                                                                  
    int Npf = Np - 2;
    FA.Xi = zero;
    decltype(peekPert(FA.Xi,0)) Xitmp(&Grid);
    Xitmp = zero;    

    //random(pRNG,Xitmp);
    QCDpt::ColourMatrix ta;
    QCDpt::LatticeComplex ca(&Grid);
    Complex shift(M_SQRT1_2,M_SQRT1_2);
    for (int a = 0; a < QCDpt::SU<Nc>::AdjointDimension; a++) {
      bernoulli(pRNG, ca);
      ca = M_SQRT2*ca - shift;
      QCDpt::SU<Nc>::generator(a, ta);
      Xitmp += ca * ta;
    }

    pokePert(FA.Xi,Xitmp,0);

    cout<<endl<<"Norm of random Xi:"<<endl;
    cout<<Pnorm2(FA.Xi)<<endl;

    cout<<endl<<"Norm of gauge field U:"<<endl;
    cout<<Pnorm2(U)<<endl;

    FA.invM(FA.psi,U,FA.Xi);
    cout<<endl<<"Norms after invM:"<<endl;
    cout<<Pnorm2(FA.psi)<<endl;

    FA.M(psi2,U,FA.psi);
    // invM works on Npf orders, so the 2 left                                                                                
    // are set to zero manually                                                                                               
    FA.psi_so = zero;
    for(int i=Npf; i<Np; i++) pokePert(psi2,FA.psi_so,i);
    cout<<endl<<"Norms after M:"<<endl;
    cout<<Pnorm2(psi2)<<endl;

    psi2 -= FA.Xi;
    cout<<endl<<"Norms after subtraction:"<<endl;
    cout<<Pnorm2(psi2)<<endl;


    cout<<endl<<"Norm of three forces:"<<endl;
    FA.deriv(U,F);
    cout<<Pnorm2(F)<<endl;
    FA.deriv(U,F);
    cout<<Pnorm2(F)<<endl;
    FA.deriv(U,F);
    cout<<Pnorm2(F)<<endl;
    
    Grid_finalize();
    return EXIT_SUCCESS;*/
    /************************************************/    
    
    
    
    
    
    //V = Logarithm(U);
    //V = Exponentiate(V);
    /*
    iVector<iScalar<iPert<iMatrix<Complex,Nc>,Np>>,Nd> link1, link2;
    PComplexD sanity1=zero, sanity2=zero;
    for (int x0=0; x0<Grid._fdimensions[0]; x0++) {
      cout << GridLogMessage << "starting x0=" << x0 << endl;
      for (int x1=0; x1<Grid._fdimensions[1]; x1++) {
        for (int x2=0; x2<Grid._fdimensions[2]; x2++) {
          for (int x3=0; x3<Grid._fdimensions[3]; x3++) {
            peekSite(link1,U,std::vector<int>({x0,x1,x2,x3}));
	    peekSite(link2,V,std::vector<int>({x0,x1,x2,x3}));
	    sanity1 += Pnorm2(link1);
	    sanity2 += Pnorm2(link2);
	    link2 -= link1;
	    PRealD diffnorm = Pnorm2(link2).real();
	    if(diffnorm(Np-1)>1.e22) cout << "bigdiff " << x0 << " " << x1 << " " << x2 << " " << x3 << endl << diffnorm << endl;
	  }
	}
      }
    }
    cout << Pnorm2(U) << endl << Pnorm2(V) << endl;
    cout << "sanity1=" << sanity1.real()/pow(16.,4) << endl;
    cout << "sanity2=" << sanity2.real()/pow(16.,4) << endl;
    */
    
    
    /*
    auto plaq1 = WilsonLoops<gimpl>::avgPlaquette(U);
    auto plaq2 = WilsonLoops<gimpl>::avgPlaquette(V);
    
    cout << GridLogMessage << "norm of fields" << endl;
    cout << Pnorm2(U) << endl;
    cout << Pnorm2(V) << endl;
    V -= U;    
    cout << Pnorm2(V) << endl;
    
    cout << GridLogMessage << "norm of plaqs" << endl;
    cout << Pnorm2(plaq1) << endl;
    cout << Pnorm2(plaq2) << endl;
    plaq2 -= plaq1;
    cout << Pnorm2(plaq2) << endl;
    */
    /*
    iVector<iScalar<iPert<iMatrix<Complex,Nc>,Np>>,Nd> link;
    iPert<iMatrix<Complex,Nc>,Np> mulink1, mulink2;
    //peekSite(link,U,std::vector<int>({12,3,8,6}));//this is the one that does not work
    peekSite(link,U,std::vector<int>({1,2,3,4}));*/
    /*mulink1 = link(0)();
    mulink2 = Logarithm(mulink1);
    mulink2 = Exponentiate(mulink2);
    mulink2 -= mulink1;
    cout << "mu=" << 0 << "\t" << Pnorm2(mulink2).real()(Np-1) << endl;*/
    
    
    /*  mulink1 = link(0)();
    cout << GridLogMessage << "MULINK1" << endl;
    
    if(Grid._processor==0){
      cout << "{";
      for(int n=0;n<Np;n++){
	cout << "{";
	for(int i=0;i<Nc;i++){
	  cout << "{";
	  for(int j=0;j<Nc;j++){
	    cout << mulink1(n)(i,j).real() << "+I*" << mulink1(n)(i,j).imag();
	    if(j==Nc-1) cout << "}"; else cout << ",";
	  }
	  if(i==Nc-1) cout << "}"; else cout <<",";
	}
	if(n==Np-1) cout << "}"; else cout <<",";
      }
      cout << endl;
    }
    
    cout << GridLogMessage << "norm before log" << endl << Pnorm2(mulink1) << endl;
    mulink2 = Logarithm(mulink1);
    cout << GridLogMessage << "norm after log" << endl << Pnorm2(mulink2) << endl;
    cout << GridLogMessage << "MULINK2" << endl;
    
    if(Grid._processor==0){
      cout << "{";
      for(int n=0;n<Np;n++){
        cout << "{";
        for(int i=0;i<Nc;i++){
	  cout << "{";
	  for(int j=0;j<Nc;j++){
            cout << mulink2(n)(i,j).real() << "+I*" << mulink2(n)(i,j).imag();
            if(j==Nc-1) cout << "}"; else cout << ",";
          }
          if(i==Nc-1) cout << "}"; else cout <<",";
        }
        if(n==Np-1) cout << "}"; else cout <<",";
      }
      cout << endl;
    }
    
    mulink2 = Exponentiate(mulink2);
    cout << GridLogMessage << "norm after exp" << endl << Pnorm2(mulink2) << endl;
    
    mulink2 -= mulink1;
    cout << GridLogMessage << "norm of the diff" << Pnorm2(mulink2) << endl;
*/    
    
    /*    
    iVector<iScalar<iPert<iMatrix<Complex,Nc>,Np>>,Nd> link1, link2;
    for (int x0=0; x0<Grid._fdimensions[0]; x0++) {
      for (int x1=0; x1<Grid._fdimensions[1]; x1++) {
	for (int x2=0; x2<Grid._fdimensions[2]; x2++) {
	  for (int x3=0; x3<Grid._fdimensions[3]; x3++) {
	    peekSite(link1,U,std::vector<int>({x0,x1,x2,x3}));
	    link2 = Logarithm(link1);
	    link2 = Exponentiate(link2);
	    link2 -= link1;
	    cout << GridLogMessage << x0 << " " << x1 << " " << x2 << " " << x3 << "\t" << Pnorm2(link2) << endl;
	    pokeSite(link2,U,std::vector<int>({x0,x1,x2,x3}));	    
	  }
	}
      }
    }
    */
    //plaq = WilsonLoops<gimpl>::avgPlaquette(U);
    //cout << plaq << endl;
    /*
    U = ProjectOnGroup(U);
    plaq = WilsonLoops<gimpl>::avgPlaquette(U);
    cout << plaq << endl;
    
    U = ProjectOnGroup(U);
    plaq = WilsonLoops<gimpl>::avgPlaquette(U);
    cout << plaq << endl;
    */
    Grid_finalize();
    return EXIT_SUCCESS;
}
