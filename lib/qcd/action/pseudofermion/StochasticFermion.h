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
  
 private:
  // vector of Wilson operators, one for each perturbative order
  std::vector<WilsonFermion<SOImpl>> Dw;
  SOGaugeField Uso, Uforce;
  std::vector<SOFermionField> psi;
  SOFermionField psitmp;
  GridParallelRNG pRNG;
  GridCartesian* grid;
  // check this number... //$//
  int Npf = Np - 2;
  

 public:
  
  /////////////////////////////////////////////////
  // Pass in required objects.
  /////////////////////////////////////////////////
  StochasticFermionAction(GridParallelRNG pRNG_, GridCartesian* grid_, GridRedBlackCartesian* rbgrid_, PRealD mass_, WilsonImplParams Params_)
      : pRNG(pRNG_),
        Uso(grid_),
        psitmp(grid_),
        Uforce(grid_),
        grid(grid_)
        {
            Dw.reserve(Npf);
            psi.reserve(Npf);
            for (int i=0; i<Npf; i++) {
                Dw.push_back(WilsonFermion<SOImpl>(Uso,*grid_,*rbgrid_,mass_(i),Params_));
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
  // dS/du = - Xi^dag  dM (M)^-1 Xi
  //////////////////////////////////////////////////////
  virtual void deriv(const GaugeField &U, GaugeField &dSdU) {
      for (int k=0; k<Npf; k++) {
          Uso = peekPert(U,k);
          Dw[k].ImportGauge(Uso);
      }
      
      SOFermionField Xi(grid);
      SOFermionField invMXi(grid);
      
      // Real part of Xi is gaussian with sigma = 1/sqrt(2),
      // same for imaginary part.
      // In this way < Xi^dag_a Xi_b > = delta_ab
      gaussian(pRNG,Xi);
      Xi *= M_SQRT1_2;
      // (are we sure that 'gaussian' generates real and imaginary parts independently?) //$//
      
      psi[0] = Xi; // here I have to apply M0^-1...//$//
      
      for (int n=1; n<Npf; n++) {
          psi[n] = zero;
          for (int j=0; j<n; j++) {
              Dw[n-j].M(psi[j],psitmp);
              psi[n] -= psitmp;
          }
          // apply M0^-1 to psi[n] //$//
      }
      
      for (int n=0; n<Npf; n++) {
          Uforce = zero;
          for (int j=0; j<=n; j++) {
              Dw[n-j].MDeriv(Uso, Xi, psi[j], DaggerNo);
              Uforce += Uso;
          }
          //$// remember algebra projection
          //$// remember to set dSdU to zero and shift the orders
          pokePert(dSdU,Uforce,n);
      }
      
      // ********************************** //
      // CHECK CONVENTIONS AND OVERALL SIGN //
      // ********************************** //
      // Our conventions really make this UdSdU; We do not differentiate wrt Udag here.
      // So must take dSdU - adj(dSdU) and left multiply by mom to get dS/dt.
      // not taking here the traceless antihermitian component*/
  };
};
}
}
}

#endif
