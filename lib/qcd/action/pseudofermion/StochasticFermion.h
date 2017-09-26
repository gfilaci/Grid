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
  
 private:
  SOGaugeField Uso, Uforce;
  GridParallelRNG *pRNG;
  GridCartesian* grid;
  TwistedFFT<Impl> TheFFT;
  
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
            Dw.reserve(Np); //$// reserving Np, needed in measuring propagator
            psi.reserve(Np); //$// reserving Np, needed in measuring propagator
            for (int i=0; i<Np; i++) { //$// reserving Np, needed in measuring propagator
                if (i==0) {
                    Dw.push_back(WilsonFermion<SOImpl>(Uso,*grid_,*rbgrid_,mass_(i),Params_));
                } else{ // there is no diagonal term in the Wilson operator
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
      gaussian(*pRNG,Xi);
      Xi *= M_SQRT1_2;
      
      // compute psi = (M)^-1 Xi
      invM(psi,U,Xi);
      
      // compute force
      dSdU = zero;
      for (int n=0; n<Npf; n++) {
          Uforce = zero;
          for (int j=0; j<=n; j++) {
              Dw[n-j].MDeriv(Uso, psi[j], Xi, DaggerYes);
              Uforce += Uso;
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

};



template <class T>
class circular_buffer{

private:
    T *buf;
    int index = 0;
    int bufsize;
    bool full = false;
    double div;
    
public:
    circular_buffer(int bufsize_):
    bufsize(bufsize_)
    {
        buf = new T[bufsize_];
        div = 1./(double)bufsize;
    }
    
    void push(T elm){
        if(index==bufsize) full = true;
        index = index % bufsize;
        buf[index++] = elm;
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
            for (int i=0; i<index; i++) {
                tmp += buf[i];
            }
            return (1./(double)index)*tmp;
        }
    }
};


template <class Impl>
class TwistValencePropagator{
    
    typedef iScalar<iMatrix<iPert<iScalar<Complex>,Np>,Ns>> PropType;
    
private:
    StochasticFermionAction<Impl> *FA;
    TwistedFFT<Impl> TheFFT1,TheFFT2;
    QCDpt::SpinColourMatrix delta;
    GridCartesian* grid;
    
public:
    circular_buffer<PropType> bufprop1;
    circular_buffer<PropType> bufprop2;
    
    TwistValencePropagator(int bufsize, GridCartesian* grid_, StochasticFermionAction<Impl> *FA_,std::vector<Complex> phases1, std::vector<Complex> phases2):
    grid(grid_),
    bufprop1(bufsize),
    bufprop2(bufsize),
    FA(FA_),
    TheFFT1(grid_,phases1),
    TheFFT2(grid_,phases2)
    {};
    
    void measurePropagator(TwistedFFT<Impl> *myFFT, PropType &propagator){
        
        for (int k=0; k<Np; k++) {
            FA->Dw[k].Params.boundary_phases = myFFT->boundary_phases;
            FA->Uso = peekPert(FA->U,k);
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
            for (int n=0; n<Np; n++){
                myFFT->FFTforward(FA->psi[n],FA->psi[n]);
                
                // put measure into proapgator
                peekSite(delta,FA->psi[n],std::vector<int>({0,0,0,0}));
                for (int alpha=0; alpha<Ns; alpha++){
                    propagator()(alpha,beta)(n)() = delta()(alpha)()(0,0);
                }
            }
            
        } // source Dirac index loop
    }
    
    void measure(){
        
        PropType prop1, prop2;
        measurePropagator(TheFFT1,prop1);
        bufprop1.push(prop1);
        measurePropagator(TheFFT2,prop2);
        bufprop2.push(prop2);
        
        // restore original phases in Dirac operators
        for (int k=0; k<Np; k++) {
            FA->Dw[k].Params.boundary_phases = FA->TheFFT->boundary_phases;
        }
        
    }
    
};


}
}
}

#endif
