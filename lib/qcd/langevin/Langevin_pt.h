/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/langevin/Langevin_pt.h

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
  
template <class Gimpl>
class PertLangevin {

protected:
    
    GridBase* grid;
    GridParallelRNG *pRNG;
    
    RealD tau;
    RealD alpha;
    double stau;
    double mtau;
    double halftau;
    double RKtau;
    Complex ci;
    ColourMatrix ta;
    
    typedef typename Gimpl::GaugeField FieldType;
    FieldType F, F0, Ftmp, U0;
    typedef decltype(peekPert(F,0)) NoiseType;
    NoiseType noise;
    typedef Lattice<iVector<iScalar<iScalar<iScalar<vReal>>>,Nd>> rndType;
    rndType ca;
    typedef decltype(peekLorentz(F,0)) GTType;
    GTType gt, agt, div;
    
public:
    ActionSet<typename Gimpl::GaugeField,PNoHirep> TheActions;
    
    explicit PertLangevin(GridBase* grid_, GridParallelRNG *pRNG_, RealD tau_, RealD alpha_):
    pRNG(pRNG_),
    tau(tau_),
    alpha(alpha_),
    F(grid_),
    F0(grid_),
    Ftmp(grid_),
    U0(grid_),
    noise(grid_),
    ca(grid_),
    gt(grid_),
    agt(grid_),
    div(grid_),
    grid(grid_),
    ci(0.,1.)
    {
        stau = std::sqrt(tau);
        mtau = - tau;
        halftau = 0.5 * tau;
        // - tau^2 * C_A/6
        RKtau = - (double)Nc * tau * tau / 6.;
        // the noise has to be normalised such that
        // < eta^a eta^b > = 2 delta^{ab}
        // therefore I need devstd = sqrt(2)
        ci *= stau * M_SQRT2;
    };

    void GenerateNoise() {
        zeroit(noise);
        for (int a = 0; a < SU<Nc>::AdjointDimension; a++) {
            gaussian(*pRNG, ca);
            SU<Nc>::generator(a, ta);
            noise += toComplex(ca) * ta;
        }
        noise *= ci;
        
    // this is slower, it's better to make GaussianFundamentalLieAlgebraMatrix explicit
//        typedef decltype(peekLorentz(noise,0)) ColourMatrixType;
//        ColourMatrixType tmp(grid);
//        for (int mu=0; mu<Nd; mu++) {
//            SU<Nc>::GaussianFundamentalLieAlgebraMatrix(pRNG, tmp, M_SQRT2);
//            pokeLorentz(noise,tmp,mu);
//        }
//        noise *= stau;
    }
    
    void StochasticGF(FieldType &U) {
        zeroit(gt);
        for (int mu=0; mu<Nd; mu++) {
            div = peekLorentz(U,mu);
            gt += div - Gimpl::MoveBackward(div,mu);
        }
        gt = alpha * Ta(gt);
        gt = Exponentiate(gt);
        agt = adj(gt);
        for(int mu=0;mu<Nd;mu++){
            // use div as a tmp field
            div = PeekIndex<LorentzIndex>(U,mu);
            div = gt * div * Gimpl::MoveForward(agt,mu);
            PokeIndex<LorentzIndex>(U,div,mu);
        }
    }

    void QuenchEulerStep(FieldType &U) {
        GenerateNoise();
        // compute drift for all actions and put it into F
        F = zero;
        for(int i=0; i<TheActions.size(); i++) {
            TheActions[i].actions[0]->deriv(U,Ftmp);
            F += Ftmp;
        }
        F = mtau * F;
        AddToOrdVoid(1,F,noise);
        U = Exponentiate(F) * U;
        StochasticGF(U);
    }
    
    void QuenchRKStep(FieldType &U) {
        GenerateNoise();
        // compute drift for all actions and put it into F
        F = zero;
        for(int i=0; i<TheActions.size(); i++) {
            TheActions[i].actions[0]->deriv(U,Ftmp);
            F += Ftmp;
        }
        F = mtau * F;
        AddToOrdVoid(1,F,noise);
        
        U0 = Exponentiate(F) * U;
        // Go from F = -tau*nabla_S[U] + stau*noise ...
        AddToOrdVoid(1,F,noise,0.5);
        // ... to F = -tau/2*nabla_S[U] + stau*noise .
        
        // starting from U0, compute drift for all actions and put it into F0
        F0 = zero;
        for(int i=0; i<TheActions.size(); i++) {
            TheActions[i].actions[0]->deriv(U0,Ftmp);
            F0 += Ftmp;
        }
        F -= halftau * F0;
        ShiftedSumVoid(2,F,F0,RKtau); // multiplication by 1/beta
        U = Exponentiate(F) * U;
        StochasticGF(U);
    }
    
    void LandauGF(FieldType &U, const double gftolerance) {
        PComplexD residualdiv;
        do{
            zeroit(gt);
            // use F as a tmp field
            F = Logarithm(U);
            for (int mu=0; mu<Nd; mu++) {
                div = peekLorentz(F,mu);
                gt += div - Gimpl::MoveBackward(div,mu);
            }
            residualdiv = Pnorm2(gt);
            // Ta(gt) is not needed, gt is already in the algebra
            gt *= alpha;
            gt = Exponentiate(gt);
            agt = adj(gt);
            for(int mu=0;mu<Nd;mu++){
                // use div as a tmp field
                div = PeekIndex<LorentzIndex>(U,mu);
                div = gt * div * Gimpl::MoveForward(agt,mu);
                PokeIndex<LorentzIndex>(U,div,mu);
            }
        } while (residualdiv(Np-1).real() > gftolerance);
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
