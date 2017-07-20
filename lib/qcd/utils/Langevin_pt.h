/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/utils/Langevin_pt.h

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
#ifndef QCD_LANGEVIN_PT_H
#define QCD_LANGEVIN_PT_H

namespace Grid {
namespace QCD {
namespace QCDpt {

  ///////////////////////////////////////////////
  // PertLangevin class
  ///////////////////////////////////////////////
  
template <class Action>
class PertLangevin {

protected:
    GridParallelRNG pRNG;
    Action ac;
    
    RealD tau;
    RealD alpha;
    double stau;
    double mtau;
    double halftau;
    double RKtau;
    
    typedef typename Action::GaugeField FieldType;
    FieldType F, F0, U0;
    typedef decltype(peekPert(F,0)) NoiseType;
    NoiseType noise;
    typedef decltype(peekLorentz(noise,0)) ColourMatrixType;
    ColourMatrixType tmp;
    typedef decltype(peekLorentz(F,0)) GTType;
    GTType gt, div;
    
public:
    explicit PertLangevin(GridBase* grid, RealD tau_, RealD alpha_):
    pRNG(grid),
    ac(1.),
    tau(tau_),
    alpha(alpha_),
    F(grid),
    F0(grid),
    U0(grid),
    noise(grid),
    tmp(grid),
    gt(grid),
    div(grid)
    {
        pRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));
        stau = std::sqrt(tau);
        mtau = - tau;
        halftau = 0.5 * tau;
        RKtau = - (double)Nc * tau * tau / 6.;
    };

    void GenerateNoise() {
        for (int mu=0; mu<Nd; mu++) {
            SU<Nc>::GaussianFundamentalLieAlgebraMatrix(pRNG, tmp, M_SQRT2);
            pokeLorentz(noise,tmp,mu);
        }
        noise *= stau;
    }
    
    void StochasticGF(FieldType &U) {
        zeroit(gt);
        for (int mu=0; mu<Nd; mu++) {
            div = peekLorentz(U,mu);
            gt += div - Cshift(div,mu,-1);
        }
        gt = alpha * Ta(gt);
        gt = Exponentiate(gt);
        SU<Nc>::GaugeTransform(U,gt);
    }
    
    void QuenchEulerStep(FieldType &U) {
        GenerateNoise();
        ac.deriv(U,F);
        F = mtau * F;
        AddToOrdVoid(1,F,noise);
        U = Exponentiate(F) * U;
        StochasticGF(U);
    }
    
    void QuenchRKStep(FieldType &U) {
        GenerateNoise();
        ac.deriv(U,F);
        F = mtau * F;
        AddToOrdVoid(1,F,noise);
        U0 = Exponentiate(F) * U;
        // Go from F = -tau*nabla_S[U] + stau*noise ...
        AddToOrdVoid(1,F,noise,0.5);
        // ... to F = -tau/2*nabla_S[U] + stau*noise .
        ac.deriv(U0,F0);
        F -= halftau * F0;
        ShiftedSumVoid(2,F,F0,RKtau); // multiplication by 1/beta
        U = Exponentiate(F) * U;
        StochasticGF(U);
    }
    
};

    
  ///////////////////////////////////////////////
  // Set to perturbative vacuum
  ///////////////////////////////////////////////
template<class vtype, int N> inline void PertVacuum(iPert<vtype,N> &P)
  {
      zeroit(P);
      typedef vtype mytype;
      mytype unit(1.0);
      P._internal[0] = unit;
  }
template<class vtype> inline void PertVacuum(iScalar<vtype>&r)
    {
      PertVacuum(r._internal);
    }
template<class vtype, int N> inline void PertVacuum(iVector<vtype,N>&r)
    {
      for (int i = 0; i < N; i++)
        PertVacuum(r._internal[i]);
    }
template<class vtype, int N> inline void PertVacuum(iMatrix<vtype,N>&r)
    {
      for (int i = 0; i < N; i++)
          for (int j = 0; j < N; j++)
              PertVacuum(r._internal[i][j]);
    }
template<class vobj>
    inline void PertVacuum(Lattice<vobj> &arg)
    {
        parallel_for(int ss=0;ss<arg._grid->oSites();ss++){
            PertVacuum(arg._odata[ss]);
        }
    }

    
  ///////////////////////////////////////////////
  // Set to random perturbative field
  ///////////////////////////////////////////////
  
template<class vobj>
    inline void PertRandom(GridParallelRNG &pRNG, Lattice<vobj> &arg)
    {
        random(pRNG,arg);
        arg = ProjectOnGroup(arg);
    }
    

}
}
}
#endif
