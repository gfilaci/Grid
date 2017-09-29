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

//#define OSX_TO_STD_RNG//$//

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
    typedef iScalar<iMatrix<iScalar<iScalar<Complex>>,Ns>> SOPropType;
    
private:
    StochasticFermionAction<Impl> *FA;
    TwistedFFT<Impl> TheFFT1,TheFFT2;
    QCDpt::SpinColourMatrix delta;
    GridCartesian* grid;
    PropType prop1, prop2;
    PropType invprop1, invprop2;
    PReal gamma1, gamma2, mc;
    Real p1sq, p2sq;
    
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
    {
        // odd powers will be discarded,
        // the series must be in "1/beta"
        assert(Np%2==1);
        
        initialiseInverseProp(invprop1,phases1);
        initialiseInverseProp(invprop2,phases2);
        
        p1sq = 0.;
        p2sq = 0.;
        for (int d=0; d<Nd; d++) {
            p1sq += (log(phases1[d]).imag()/(double)grid->_fdimensions[d])*(log(phases1[d]).imag()/(double)grid->_fdimensions[d]);
            p2sq += (log(phases2[d]).imag()/(double)grid->_fdimensions[d])*(log(phases2[d]).imag()/(double)grid->_fdimensions[d]);
        }
    };
    
    template <class gaugefield>
    void measurePropagator(TwistedFFT<Impl> *myFFT, PropType &propagator, const gaugefield &U){
        propagator = zero;
        
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
                peekSite(delta,FA->psi[n],std::vector<int>({0,0,0,0}));
                for (int alpha=0; alpha<Ns; alpha++){
                    propagator()(alpha,beta)(n)() = delta()(alpha)()(0,0);
                }
            }
            
        } // source Dirac index loop
    }
    
    template <class gaugefield>
    void measure(const gaugefield &U){
        
        measurePropagator(&TheFFT1,prop1,U);
        bufprop1.push(prop1);
        measurePropagator(&TheFFT2,prop2,U);
        bufprop2.push(prop2);
        
        // restore original phases in Dirac operators
        for (int k=0; k<Np; k++) {
            FA->Dw[k].Params.boundary_phases = FA->TheFFT.boundary_phases;
        }
        
    }
    
    void feedback_cm(std::ofstream &massfile){
        prop1 = bufprop1.average();
        prop2 = bufprop2.average();
        invertpropagator(invprop1,prop1);
        invertpropagator(invprop2,prop2);
        
        gamma1 = (1./(double)Nd) * TensorRemove(traceIndex<SpinIndex>(invprop1)).real();
        gamma2 = (1./(double)Nd) * TensorRemove(traceIndex<SpinIndex>(invprop2)).real();
        
        mc = (1./(p2sq-p1sq))*(p2sq*gamma1-p1sq*gamma2);
        
        // higher orders had mass-4 to kill the Wilson diagonal term,
        // feedback only on "1/beta" components
        for (int i=2; i<Np; i+=2) {
            FA->Dw[i].mass -= mc(i);
            massfile << FA->Dw[i].mass+4. << std::endl;
        }
    }
    
    void invertpropagator(PropType &invprop, PropType &prop){
        std::vector<SOPropType> inp(Np), outp(Np);
        for (int k=0; k<Np; k++) {
            inp[k] = peekIndex<PertIndex>(prop,k);
        }
        outp[0] = peekIndex<PertIndex>(invprop,0);
        for (int k=1; k<Np; k++) {
            outp[k] = zero;
            for (int j=0; j<k; j++) {
                outp[k] += inp[k-j]*outp[j];
            }
            outp[k] = -outp[0]*outp[k];
        }
        for (int k=1; k<Np; k++) {
            pokePert(invprop,outp[k],k);
        }
        /******************/
//        //test...
//        SOPropType test;
//        for (int i=0; i<Np; i++) {
//            test = zero;
//            for (int j=0; j<=i; j++) {
//                test += inp[j]*outp[i-j];
//            }
//            if(i==0) std::cout << GridLogMessage << "inversion order " << i << " good within " << norm2(test)-4. << std::endl;
//            else std::cout << GridLogMessage << "inversion order " << i << " good within " << norm2(test) << std::endl;
//        }
        /******************/
    }
    
    void initialiseInverseProp(PropType &invprop, const std::vector<Complex> &phases){
        
        Complex im(0.,1.);
        std::vector<Real> mom({0,0,0,0});
        for (int d=0; d<Nd; d++) {
            mom[d] = log(phases[d]).imag()/(double)grid->_fdimensions[d];
        }
        // sum_mu sin^2(p_mu/2)
        Real mompr = 0.;
        for (int mu=0; mu<Nd; mu++) {
            mompr += std::sin(0.5*mom[mu])*std::sin(0.5*mom[mu]);
        }
        
        
        Gamma::Algebra Gmu [] = {
            Gamma::Algebra::GammaX,
            Gamma::Algebra::GammaY,
            Gamma::Algebra::GammaZ,
            Gamma::Algebra::GammaT
        };
        
        iScalar<iMatrix<iScalar<Complex>,Nd>> mygamma[Nd];
        for (int mu=0; mu<Nd; mu++) {
            mygamma[mu] = zero;
            for (int i=0; i<Nd; i++) {
                mygamma[mu]()(i,i)() = 1.0;
            }
            mygamma[mu] = QCD::Gamma(Gmu[mu])*mygamma[mu];
        }
        
        SOPropType Pmygamma[Nd];
        for (int mu=0; mu<Nd; mu++) {
            for (int i=0; i<Nd; i++) {
                for (int j=0; j<Nd; j++) {
                    Pmygamma[mu]()(i,j)()() = mygamma[mu]()(i,j)();
                }
            }
        }
        
        SOPropType SOprop = zero;
        
        // Fill order 0 with
        // i ( gamma_mu sin(p_mu) ) + ...
        for (int mu=0; mu<Nd; mu++) {
            SOprop += im * std::sin(mom[mu]) * Pmygamma[mu];
        }
        
        // ... + 2 * sum_mu sin^2(p_mu/2)
        for (int i=0; i<Nd; i++) {
            SOprop()(i,i)()() += 2*mompr;
        }
        
        
        invprop = zero;
        for (int i=0; i<Nd; i++) {
            for (int j=0; j<Nd; j++) {
                invprop()(i,j)(0)() = SOprop()(i,j)()();
            }
        }
    }
};


}
}
}

#endif
