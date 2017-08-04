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
template <class Impl, int Nf>
class StochasticFermionAction : public Action<typename Impl::GaugeField> {
 public:
  INHERIT_IMPL_TYPES(Impl);

 private:
  FermionOperator<Impl> &FermOp;  // the basic operator
  const double invNf = 1. / (double)Nf;

//  OperatorFunction<FermionField> &DerivativeSolver;//$//

  FermionField Phi;

 public:
  /////////////////////////////////////////////////
  // Pass in required objects.
  /////////////////////////////////////////////////
  StochasticFermionAction(FermionOperator<Impl> &Op/*,
                                OperatorFunction<FermionField> &DS*/)
      : FermOp(Op),
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
      
      FermionField Xi(FermOp.FermionGrid());
      FermionField invMXi(FermOp.FermionGrid());
      
      // generate random Xi
      // ...
      
      // compute invMXi = (M)^-1 Xi
//      invMXi = zero;
//      DerivativeSolver(M, Xi, invMXi);
      
      // compute force
//      FermOp.MDeriv(dSdU, Xi, invMXi, DaggerNo);
//      dSdU *= invNf;
      
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
