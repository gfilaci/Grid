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
// Stochastic fermion action for any dop
////////////////////////////////////////////////////////////////////////
template <class Impl>
class StochasticFermionAction : public Action<typename Impl::GaugeField> {
 public:
  INHERIT_IMPL_TYPES(Impl);

 private:
  FermionOperator<Impl> &FermOp;  // the basic operator
  GridParallelRNG pRNG;
  
//  OperatorFunction<FermionField> &DerivativeSolver;//$//

  FermionField Phi;

 public:
  /////////////////////////////////////////////////
  // Pass in required objects.
  /////////////////////////////////////////////////
  StochasticFermionAction(FermionOperator<Impl> &Op, GridParallelRNG pRNG_/*,
                                OperatorFunction<FermionField> &DS*/)
      : FermOp(Op),
        pRNG(pRNG_),
        //DerivativeSolver(DS),
        Phi(Op.FermionGrid()){};


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
      FermOp.ImportGauge(U);
      
      //$// the noise should not be perturbative...
      FermionField Xi(FermOp.FermionGrid());
      FermionField invMXi(FermOp.FermionGrid());
      decltype(peekPert(Xi,0)) tmp(FermOp.FermionGrid());
      
      // Real part of Xi is gaussian with sigma = 1/sqrt(2),
      // same for imaginary part.
      // In this way < Xi^dag_a Xi_b > = delta_ab
      gaussian(pRNG,tmp);
      tmp *= M_SQRT1_2;
      pokePert(Xi,tmp,0);
      // (is it sure that gaussian generates real and imaginary parts independently?) //$//
      
      // compute invMXi = (M)^-1 Xi
//      invMXi = zero;
//      DerivativeSolver(M, Xi, invMXi);
      
      invMXi = Xi; //$// just to have something, waiting to be able to compute the inverse...
      
      // compute force
      FermOp.MDeriv(dSdU, Xi, invMXi, DaggerNo);
      //$// remember algebra projection
      
      // ********************************** //
      // CHECK CONVENTIONS AND OVERALL SIGN //
      // ********************************** //
      // Our conventions really make this UdSdU; We do not differentiate wrt Udag here.
      // So must take dSdU - adj(dSdU) and left multiply by mom to get dS/dt.
      // not taking here the traceless antihermitian component
  };
};
}
}
}

#endif
