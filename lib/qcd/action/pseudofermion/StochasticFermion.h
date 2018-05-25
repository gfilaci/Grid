/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/pseudofermion/StochasticFermion.h

Copyright (C) 2015

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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
#ifndef QCD_STOCHASTICFERMION_H
#define QCD_STOCHASTICFERMION_H

//#define OSX_TO_STD_RNG//$$//

namespace Grid {
namespace QCD {
namespace QCDpt {

////////////////////////////////////////////////////////////////////////
// Stochastic fermion action for Wilson fermions
////////////////////////////////////////////////////////////////////////
template <class Impl>
class StochasticFermionAction : public Action<typename Impl::GaugeField> {
 public:
  INHERIT_IMPL_TYPES(Impl);
  typedef typename Impl::SOimpl SOImpl;
  typedef typename Impl::SOimpl::FermionField SOFermionField;
  typedef typename Impl::SOimpl::Gimpl::Field SOGaugeField;
  
  // vector of Wilson operators, one for each perturbative order
  std::vector<WilsonFermion<SOImpl>> Dw;
  std::vector<SOFermionField> psi;
  SOFermionField Xi, psitmp;
  SOGaugeField Uso, Uforce;
  TwistedFFT<Impl> TheFFT;
  
 private:
  GridParallelRNG *pRNG;
  GridCartesian* grid;
  
  
  // we can work at lower perturabtive order, because
  // there is a factor 1/beta in front of the fermion drift.
  // it arises from rescaling tau = beta * epsilon
  int Npf = Np - 2;
  
  int Nf; // number of flavours
  double Nf_over_Nc;  // 1/Nc is the factor due to the smell

 public:
  
  /////////////////////////////////////////////////
  // Pass in required objects.
  /////////////////////////////////////////////////
  StochasticFermionAction(GridParallelRNG *pRNG_, GridCartesian* grid_, GridRedBlackCartesian* rbgrid_, PRealD mass_, WilsonImplParams Params_, int Nf_)
      : pRNG(pRNG_),
        Uso(grid_),
        Xi(grid_),
        psitmp(grid_),
        Uforce(grid_),
        grid(grid_),
        Nf(Nf_),
        TheFFT(grid_,Params_.boundary_phases)
        {
            Nf_over_Nc = (double)Nf / (double)Nc;
            Dw.reserve(Np); // reserving Np, needed in measuring propagator
            psi.reserve(Np); // reserving Np, needed in measuring propagator
            for (int i=0; i<Np; i++) { // initialising Np, needed in measuring propagator
                if (i==0) {
                    Dw.push_back(WilsonFermion<SOImpl>(Uso,*grid_,*rbgrid_,mass_(i),Params_));
                } else{ // else there is no diagonal term in the Wilson operator
                    Dw.push_back(WilsonFermion<SOImpl>(Uso,*grid_,*rbgrid_,mass_(i)-4.,Params_));
                }
                psi.push_back(SOFermionField(grid));
            }
        };
  
  virtual std::string action_name(){return "StochasticFermionAction";}

  virtual std::string LogParameters(){
    std::stringstream sstream;
    sstream << GridLogMessage << "["<<action_name()<<"] has no parameters" << std::endl;
    return sstream.str();
  }  
  
  virtual void refresh(const GaugeField &U, GridParallelRNG &pRNG) {
      assert(0);
  };
  
  virtual PRealD S(const GaugeField &U) {
      assert(0);
  };

  //////////////////////////////////////////////////////
  // return dSdU = i nabla S = -i Xi^dag  dM (M)^-1 Xi
  //////////////////////////////////////////////////////
  virtual void deriv(const GaugeField &U, GaugeField &dSdU) {
      
      // Real part of Xi is gaussian with sigma = 1/sqrt(2),
      // same for imaginary part.
      // In this way < Xi^dag_a Xi_b > = delta_ab
      
      // GAUSSIAN NOISE
      gaussian(*pRNG,Xi);
      Xi *= M_SQRT1_2;
      
      // Z2 NOISE
//      Complex shift(M_SQRT1_2,M_SQRT1_2);
//      typename Impl::SOimpl::SiteSpinor shiftmat;
//      for(int a=0; a<Nc; a++) for(int b=0; b<Nc; b++) pokeColour(shiftmat,shift,a,b);
//      bernoulli(*pRNG,Xi);
//      Xi = M_SQRT2*Xi - shiftmat;
      
#ifdef OSX_TO_STD_RNG
      Complex im(0.,1.);
      Xi = imag(Xi) + im * real(Xi);
#endif
      
      // compute psi = (M)^-1 Xi
      invM(psi,U,Xi);
      
      // compute force
      dSdU = zero;
      for (int n=0; n<Npf; n++) {
          Uforce = zero;
          for (int j=0; j<=n; j++) {
              Dw[n-j].MDeriv(Uso, psi[j], Xi, DaggerYes);
              Uforce += Uso;
              // the following three lines enable the non optimised fermion drift
//              Dw[n-j].MDeriv(Uso, Xi, psi[j], DaggerNo);
//              Uforce += Uso;
//              Uforce *= 0.5;
          }
          // the "+2" shift is due to the 1/beta factor in front of the fermion drift.
          // the first two orders are already set to zero.
          pokePert(dSdU,Uforce,n+2);
      }
      
      dSdU = Ta(dSdU);
      dSdU *= Nf_over_Nc;
      
  }
  
  
    // COMPUTE PERTURBATIVELY psi = (M)^-1 Xi
    void invM(std::vector<SOFermionField>& psi, const GaugeField &U, const SOFermionField &Xi){
        for (int k=0; k<Npf; k++) {
            Uso = peekPert(U,k);
            Dw[k].ImportGauge(Uso);
        }
        
        // apply M0^-1 to Xi
        TheFFT.FreeWilsonOperatorInverse(psi[0],Xi);
        
        for (int n=1; n<Npf; n++){
            psi[n] = zero;
            for (int j=0; j<n; j++) {
                Dw[n-j].M(psi[j],psitmp);
                psi[n] -= psitmp;
            }
            // apply M0^-1 to psi[n]
            TheFFT.FreeWilsonOperatorInverse(psi[n],psi[n]);
        }
    }
    
    
    // COMPUTE PERTURBATIVELY psi2 = M psi
    void M(std::vector<SOFermionField>& psi2, const GaugeField &U, const std::vector<SOFermionField>& psi){
        for (int k=0; k<Npf; k++) {
            Uso = peekPert(U,k);
            Dw[k].ImportGauge(Uso);
        }
        
        for (int n=0; n<Npf; n++){
            psi2[n] = zero;
            for (int j=0; j<=n; j++){
                Dw[j].M(psi[n-j],psitmp);
                psi2[n] += psitmp;
            }
        }
    }
};


////////////////////////////////////////////////////////////////////////
// Stochastic fermion action for staggered fermions
////////////////////////////////////////////////////////////////////////
template <class Impl>
class StochasticStaggeredAction : public Action<typename Impl::GaugeField> {
 public:
  INHERIT_IMPL_TYPES(Impl);
  typedef typename Impl::SOimpl SOImpl;
  typedef typename Impl::SOimpl::FermionField SOFermionField;
  typedef typename Impl::SOimpl::Gimpl::Field SOGaugeField;
  
  // vector of staggered operators, one for each perturbative order
  std::vector<StaggeredFermion<SOImpl>> Dw;
  std::vector<SOFermionField> psi;
  SOFermionField Xi, psitmp;
  SOGaugeField Uso, Uforce;
  TwistedStaggeredFFT<Impl> TheFFT;
  
 private:
  GridParallelRNG *pRNG;
  GridCartesian* grid;
  
  
  // we can work at lower perturabtive order, because
  // there is a factor 1/beta in front of the fermion drift.
  // it arises from rescaling tau = beta * epsilon
  int Npf = Np - 2;
  
  int Nf; // number of flavours
  double Nf_over_fourNc;  // 1/Nc is the factor due to the smell, 1/4 is for taste

 public:
  
  /////////////////////////////////////////////////
  // Pass in required objects.
  /////////////////////////////////////////////////
  StochasticStaggeredAction(GridParallelRNG *pRNG_, GridCartesian* grid_, GridRedBlackCartesian* rbgrid_, RealD mass_, WilsonImplParams Params_, int Nf_)
      : pRNG(pRNG_),
        Uso(grid_),
        Xi(grid_),
        psitmp(grid_),
        Uforce(grid_),
        grid(grid_),
        Nf(Nf_),
        TheFFT(grid_,mass_,Params_.boundary_phases)
        {
            Nf_over_fourNc = (double)Nf / (4.*(double)Nc);
            Dw.reserve(Np); // reserving Np, needed in measuring propagator
            psi.reserve(Np); // reserving Np, needed in measuring propagator
            for (int i=0; i<Np; i++) { // initialising Np, needed in measuring propagator
                if (i==0) {
                    Dw.push_back(StaggeredFermion<SOImpl>(Uso,*grid_,*rbgrid_,mass_,Params_));
                } else{ // else there is no diagonal term in the Wilson operator
                    Dw.push_back(StaggeredFermion<SOImpl>(Uso,*grid_,*rbgrid_,0.,Params_));
                }
                psi.push_back(SOFermionField(grid));
            }
        };
  
  virtual std::string action_name(){return "StochasticStaggeredAction";}

  virtual std::string LogParameters(){
    std::stringstream sstream;
    sstream << GridLogMessage << "["<<action_name()<<"] has no parameters" << std::endl;
    return sstream.str();
  }  
  
  virtual void refresh(const GaugeField &U, GridParallelRNG &pRNG) {
      assert(0);
  };
  
  virtual PRealD S(const GaugeField &U) {
      assert(0);
  };

  //////////////////////////////////////////////////////
  // return dSdU = i nabla S = -i Xi^dag  dM (M)^-1 Xi
  //////////////////////////////////////////////////////
  virtual void deriv(const GaugeField &U, GaugeField &dSdU) {
      
      // Real part of Xi is gaussian with sigma = 1/sqrt(2),
      // same for imaginary part.
      // In this way < Xi^dag_a Xi_b > = delta_ab
      
      // GAUSSIAN NOISE
      gaussian(*pRNG,Xi);
      Xi *= M_SQRT1_2;
      
      // Z2 NOISE
//      Complex shift(M_SQRT1_2,M_SQRT1_2);
//      typename Impl::SOimpl::SiteSpinor shiftmat;
//      for(int a=0; a<Nc; a++) for(int b=0; b<Nc; b++) pokeColour(shiftmat,shift,a,b);
//      bernoulli(*pRNG,Xi);
//      Xi = M_SQRT2*Xi - shiftmat;
      
#ifdef OSX_TO_STD_RNG
      Complex im(0.,1.);
      Xi = imag(Xi) + im * real(Xi);
#endif
      
      // compute psi = (M)^-1 Xi
      invM(psi,U,Xi);
      
      // compute force
      dSdU = zero;
      for (int n=0; n<Npf; n++) {
          Uforce = zero;
          for (int j=0; j<=n; j++) {
              Dw[n-j].MDeriv(Uso, psi[j], Xi, DaggerYes);
              Uforce += Uso;
              // the following three lines enable the non optimised fermion drift
//              Dw[n-j].MDeriv(Uso, Xi, psi[j], DaggerNo);
//              Uforce += Uso;
//              Uforce *= 0.5;
          }
          // the "+2" shift is due to the 1/beta factor in front of the fermion drift.
          // the first two orders are already set to zero.
          pokePert(dSdU,Uforce,n+2);
      }
      
      dSdU = Ta(dSdU);
      dSdU *= Nf_over_fourNc;
      
  }
  
  
    // COMPUTE PERTURBATIVELY psi = (M)^-1 Xi
    void invM(std::vector<SOFermionField>& psi, const GaugeField &U, const SOFermionField &Xi){
        for (int k=0; k<Npf; k++) {
            Uso = peekPert(U,k);
            Dw[k].ImportGauge(Uso);
        }

        // apply M0^-1 to Xi
        TheFFT.FreeStaggeredOperatorInverse(psi[0],Xi);

        for (int n=1; n<Npf; n++){
            psi[n] = zero;
            for (int j=0; j<n; j++) {
                Dw[n-j].M(psi[j],psitmp);
                psi[n] -= psitmp;
            }
            // apply M0^-1 to psi[n]
            TheFFT.FreeStaggeredOperatorInverse(psi[n],psi[n]);
        }
    }
    
    // COMPUTE PERTURBATIVELY psi2 = M psi
    void M(std::vector<SOFermionField>& psi2, const GaugeField &U, const std::vector<SOFermionField>& psi){
        for (int k=0; k<Npf; k++) {
            Uso = peekPert(U,k);
            Dw[k].ImportGauge(Uso);
        }
        
        for (int n=0; n<Npf; n++){
            psi2[n] = zero;
            for (int j=0; j<=n; j++){
                Dw[j].M(psi[n-j],psitmp);
                psi2[n] += psitmp;
            }
        }
    }
};

///////////////////////////////////////////////////////////////////////////////////
// Stochastic fermion action for staggered fermions in the adjoint representation
///////////////////////////////////////////////////////////////////////////////////
template <class Impl>
class StochasticAdjointStaggeredAction : public Action<typename Impl::GaugeField> {
public:
    INHERIT_IMPL_TYPES(Impl);
    typedef typename Impl::SOimpl SOImpl;
    typedef typename Impl::SOimpl::FermionField SOFermionField;
    typedef typename Impl::SOimpl::Gimpl::Field SOGaugeField;
    
    // vector of staggered operators, one for each perturbative order
    std::vector<StaggeredFermion<SOImpl>> Dw;
    std::vector<SOFermionField> psi;
    SOFermionField Xi, psitmp;
    SOGaugeField Uso, Uforce;
    TwistedStaggeredFFT<Impl> TheFFT;
    
private:
    GridParallelRNG *pRNG;
    GridCartesian* grid;
    
    
    // we can work at lower perturabtive order, because
    // there is a factor 1/beta in front of the fermion drift.
    // it arises from rescaling tau = beta * epsilon
    int Npf = Np - 2;
    
    int Nf; // number of flavours
    double twoNf_over_four;  // 1/4 is for taste, 2 is for the trace normalisation since the field is in the matrix representation
    
public:
    
    /////////////////////////////////////////////////
    // Pass in required objects.
    /////////////////////////////////////////////////
    StochasticAdjointStaggeredAction(GridParallelRNG *pRNG_, GridCartesian* grid_, GridRedBlackCartesian* rbgrid_, RealD mass_, WilsonImplParams Params_, int Nf_)
    : pRNG(pRNG_),
    Uso(grid_),
    Xi(grid_),
    psitmp(grid_),
    Uforce(grid_),
    grid(grid_),
    Nf(Nf_),
    TheFFT(grid_,mass_,Params_.boundary_phases)
    {
        twoNf_over_four = 2. * (double)Nf / 4.;
        Dw.reserve(Np); // reserving Np, needed in measuring propagator
        psi.reserve(Np); // reserving Np, needed in measuring propagator
        for (int i=0; i<Np; i++) { // initialising Np, needed in measuring propagator
            if (i==0) {
                Dw.push_back(StaggeredFermion<SOImpl>(Uso,*grid_,*rbgrid_,mass_,Params_));
            } else{ // else there is no diagonal term in the Wilson operator
                Dw.push_back(StaggeredFermion<SOImpl>(Uso,*grid_,*rbgrid_,0.,Params_));
            }
            psi.push_back(SOFermionField(grid));
        }
    };
    
    virtual std::string action_name(){return "StochasticAdjointStaggeredAction";}
    
    virtual std::string LogParameters(){
        std::stringstream sstream;
        sstream << GridLogMessage << "["<<action_name()<<"] has no parameters" << std::endl;
        return sstream.str();
    }
    
    virtual void refresh(const GaugeField &U, GridParallelRNG &pRNG) {
        assert(0);
    };
    
    virtual PRealD S(const GaugeField &U) {
        assert(0);
    };
    
    //////////////////////////////////////////////////////
    // return dSdU = i nabla S = -i Xi^dag  dM (M)^-1 Xi
    //////////////////////////////////////////////////////
    virtual void deriv(const GaugeField &U, GaugeField &dSdU) {
        
//        // complex gaussian noise in the algebra
//        ColourMatrix ta;
//        Xi = zero;
//        LatticeComplex ca(grid);
//        for (int a = 0; a < SU<Nc>::AdjointDimension; a++) {
//            gaussian(*pRNG, ca);
//            SU<Nc>::generator(a, ta);
//            Xi += ca * ta;
//        }
//        Xi *= M_SQRT1_2;//$//
        
        // real gaussian noise in the algebra
        ColourMatrix ta;
        Xi = zero;
        LatticeReal ca(grid);
        for (int a = 0; a < SU<Nc>::AdjointDimension; a++) {
            gaussian(*pRNG, ca);
            SU<Nc>::generator(a, ta);
            Xi += toComplex(ca) * ta;
        }
        //$//
        
        // Real part of Xi is gaussian with sigma = 1/sqrt(2),
        // same for imaginary part.
        // In this way < Xi^dag_a Xi_b > = delta_ab
        
        // GAUSSIAN NOISE
//        gaussian(*pRNG,Xi);
//        Xi *= M_SQRT1_2;
        
        // Z2 NOISE
        //      Complex shift(M_SQRT1_2,M_SQRT1_2);
        //      typename Impl::SOimpl::SiteSpinor shiftmat;
        //      for(int a=0; a<Nc; a++) for(int b=0; b<Nc; b++) pokeColour(shiftmat,shift,a,b);
        //      bernoulli(*pRNG,Xi);
        //      Xi = M_SQRT2*Xi - shiftmat;
        
#ifdef OSX_TO_STD_RNG
        Complex im(0.,1.);
        Xi = imag(Xi) + im * real(Xi);
#endif
        
        // compute psi = (M)^-1 Xi
        invM(psi,U,Xi);
        
        // compute force
        dSdU = zero;
        for (int n=0; n<Npf; n++) {
            Uforce = zero;
            for (int j=0; j<=n; j++) {
                Dw[n-j].MDeriv(Uso, psi[j], Xi, DaggerYes);
                Uforce += Uso;
                // the following three lines enable the non optimised fermion drift
                Dw[n-j].MDeriv(Uso, Xi, psi[j], DaggerNo);
                Uforce += Uso;
                Uforce *= 0.5;
            }
            // the "+2" shift is due to the 1/beta factor in front of the fermion drift.
            // the first two orders are already set to zero.
            pokePert(dSdU,Uforce,n+2);
        }
        dSdU = Ta(dSdU);
        dSdU *= twoNf_over_four;
        
    }
    
    
    // COMPUTE PERTURBATIVELY psi = (M)^-1 Xi
    void invM(std::vector<SOFermionField>& psi, const GaugeField &U, const SOFermionField &Xi){
        for (int k=0; k<Npf; k++) {
            Uso = peekPert(U,k);
            Dw[k].ImportGauge(Uso);
        }

        // apply M0^-1 to Xi
        TheFFT.FreeStaggeredOperatorInverse(psi[0],Xi);

        for (int n=1; n<Npf; n++){
            psi[n] = zero;
            for (int j=0; j<n; j++) {
                Dw[n-j].M(psi[j],psitmp);
                psi[n] -= psitmp;
            }
            // apply M0^-1 to psi[n]
            TheFFT.FreeStaggeredOperatorInverse(psi[n],psi[n]);
        }
    }
    
    // COMPUTE PERTURBATIVELY psi2 = M psi
    void M(std::vector<SOFermionField>& psi2, const GaugeField &U, const std::vector<SOFermionField>& psi){
        for (int k=0; k<Npf; k++) {
            Uso = peekPert(U,k);
            Dw[k].ImportGauge(Uso);
        }
        
        for (int n=0; n<Npf; n++){
            psi2[n] = zero;
            for (int j=0; j<=n; j++){
                Dw[j].M(psi[n-j],psitmp);
                psi2[n] += psitmp;
            }
        }
    }
};

template <class T>
class circular_buffer{

private:
    std::vector<T> buf;
    int index = -1;
    int bufsize;
    bool full = false;
    double div;
    
public:
    circular_buffer(int bufsize_):
    bufsize(bufsize_)
    {
        buf.reserve(bufsize_);
        div = 1./(double)bufsize_;
    }
    
    void push(T elm){
        if((index+1)==bufsize) full = true;
        index = (index+1) % bufsize;
        buf[index] = elm;
    }
    
    T pull(int i){
        return buf[(index+i)%bufsize];
    }
    
    T average(){
        T tmp = zero;
        if(full){
            for (int i=0; i<bufsize; i++) {
                tmp += buf[i];
            }
            return div*tmp;
        } else{
            for (int i=0; i<=index; i++) {
                tmp += buf[i];
            }
            return (1./(double)(index+1))*tmp;
        }
    }
};


template <class Impl>
class TwistValencePropagator{
    
    typedef iScalar<iMatrix<iPert<iScalar<Complex>,Np>,Ns>> PropType;
    
private:
    StochasticFermionAction<Impl> *FA;
    TwistedFFT<Impl> TheFFT1,TheFFT2,TheFFT3;
    QCDpt::SpinColourMatrix delta;
    GridCartesian* grid;
    PropType prop;
    
public:
    
    TwistValencePropagator(GridCartesian* grid_, StochasticFermionAction<Impl> *FA_,std::vector<Complex> phases1, std::vector<Complex> phases2, std::vector<Complex> phases3):
    grid(grid_),
    FA(FA_),
    TheFFT1(grid_,phases1),
    TheFFT2(grid_,phases2),
    TheFFT3(grid_,phases3)
    {
        // odd powers will be discarded,
        // the series must be in "1/beta"
        assert(Np%2==1);
    };
    
    
    // N.B.: the propagator file must be read with
    // loop #propagators
    // loop three momenta
    // loop beta
    // loop orders (in "1/beta")
    // loop alpha
    // read Complex
    template <class gaugefield>
    void measurePropagator(TwistedFFT<Impl> *myFFT, std::ofstream &propfile, const gaugefield &U){
        
        for (int k=0; k<Np; k++) {
            FA->Dw[k].Params.boundary_phases = myFFT->boundary_phases;
            FA->Uso = peekPert(U,k);
            FA->Dw[k].ImportGauge(FA->Uso);
        }
        
        for (int beta=0; beta<Ns; beta++){
            FA->Xi = zero;
            delta = zero;
            delta()(beta)()(0,0) = 1.0;
            pokeSite(delta,FA->Xi,std::vector<int>({0,0,0,0}));
            
            // apply M0^-1 to Xi
            // Xi is already in momentum space,
            // no need for an FFT forward
            myFFT->TwistedWilsonPropagator(FA->psi[0],FA->Xi);
            myFFT->FFTbackward(FA->psi[0],FA->psi[0]);
            
            for (int n=1; n<Np; n++){
                FA->psi[n] = zero;
                for (int j=0; j<n; j++) {
                    FA->Dw[n-j].M(FA->psi[j],FA->psitmp);
                    FA->psi[n] -= FA->psitmp;
                }
                // apply M0^-1 to psi[n]
                myFFT->FreeWilsonOperatorInverse(FA->psi[n],FA->psi[n]);
            }
            
            // back to momentum space
            // save only "1/beta orders"
            for (int n=0; n<Np; n+=2){
                myFFT->FFTforward(FA->psi[n],FA->psi[n]);
                // put measure into proapgator
                // use delta as tmp storage
                peekSite(delta,FA->psi[n],std::vector<int>({0,0,0,0}));
                for (int alpha=0; alpha<Ns; alpha++){
                    propfile.write( reinterpret_cast<char*>( &(delta()(alpha)()(0,0)) ), sizeof(Complex) );
                }
            }
            
        } // source Dirac index loop
        propfile.flush();
    }
    
    template <class gaugefield>
    void measure(std::ofstream &propfile, const gaugefield &U){
        
        measurePropagator(&TheFFT1,propfile,U);
        measurePropagator(&TheFFT2,propfile,U);
        measurePropagator(&TheFFT3,propfile,U);
        
        // restore original phases in Dirac operators
        for (int k=0; k<Np; k++) {
            FA->Dw[k].Params.boundary_phases = FA->TheFFT.boundary_phases;
        }
    }
    
};


}
}
}

#endif
