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
  
 private:
  SOGaugeField Uso, Uforce;
  std::vector<SOFermionField> psi;
  SOFermionField Xi, psitmp;
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
            Dw.reserve(Npf);
            psi.reserve(Npf);
            for (int i=0; i<Npf; i++) {
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
        
        for (int n=1; n<Npf; n++) {
            psi[n] = zero;
            for (int j=0; j<n; j++) {
                Dw[n-j].M(psi[j],psitmp);
                psi[n] -= psitmp;
            }
            // apply M0^-1 to psi[n]
            TheFFT.FreeWilsonOperatorInverse(psi[n],psi[n]);
        }
    }

void changeBoundaryPhases(std::vector<Complex> newphases){
        for (int k=0; k<Npf; k++) {
            Dw[k].Params.boundary_phases = newphases;
        }
        TheFFT.FFTinitialisation(newphases);
}

};


}
}
}

#endif
