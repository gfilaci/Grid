/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/fermion/StaggeredFermion.h

Copyright (C) 2015

Author: Azusa Yamaguchi, Peter Boyle, Gianluca Filaci

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
#ifndef GRID_QCD_STAG_FERMION_H
#define GRID_QCD_STAG_FERMION_H

namespace Grid {

namespace QCD {

class StaggeredFermionStatic {
 public:
  static const std::vector<int> directions;
  static const std::vector<int> displacements;
  static const int npoint = 8;
};

template <class Impl>
class StaggeredFermion : public NaiveStaggeredKernels<Impl>, public StaggeredFermionStatic {
 public:
  INHERIT_IMPL_TYPES(Impl);
  typedef NaiveStaggeredKernels<Impl> Kernels;

  FermionField _tmp;
  FermionField &tmp(void) { return _tmp; }

  ///////////////////////////////////////////////////////////////
  // Implement the abstract base
  ///////////////////////////////////////////////////////////////
  GridBase *GaugeGrid(void) { return _grid; }
  GridBase *GaugeRedBlackGrid(void) { return _cbgrid; }
  GridBase *FermionGrid(void) { return _grid; }
  GridBase *FermionRedBlackGrid(void) { return _cbgrid; }

  //////////////////////////////////////////////////////////////////
  // override multiply; cut number routines if pass dagger argument
  // and also make interface more uniformly consistent
  //////////////////////////////////////////////////////////////////
  RealD M(const FermionField &in, FermionField &out);
  RealD Mdag(const FermionField &in, FermionField &out);

  /////////////////////////////////////////////////////////
  // half checkerboard operations
  /////////////////////////////////////////////////////////
  void Meooe(const FermionField &in, FermionField &out);
  void MeooeDag(const FermionField &in, FermionField &out);
  void Mooee(const FermionField &in, FermionField &out);
  void MooeeDag(const FermionField &in, FermionField &out);
  void MooeeInv(const FermionField &in, FermionField &out);
  void MooeeInvDag(const FermionField &in, FermionField &out);

  ////////////////////////
  // Derivative interface
  ////////////////////////
  // Interface calls an internal routine
  void DhopDeriv  (GaugeField &mat, const FermionField &U, const FermionField &V, int dag);
  void DhopDerivOE(GaugeField &mat, const FermionField &U, const FermionField &V, int dag);
  void DhopDerivEO(GaugeField &mat, const FermionField &U, const FermionField &V, int dag);

  ///////////////////////////////////////////////////////////////
  // non-hermitian hopping term; half cb or both
  ///////////////////////////////////////////////////////////////
  void Dhop  (const FermionField &in, FermionField &out, int dag);
  void DhopOE(const FermionField &in, FermionField &out, int dag);
  void DhopEO(const FermionField &in, FermionField &out, int dag);

  ///////////////////////////////////////////////////////////////
  // Multigrid assistance; force term uses too
  ///////////////////////////////////////////////////////////////
  void Mdir(const FermionField &in, FermionField &out, int dir, int disp);
  void DhopDir(const FermionField &in, FermionField &out, int dir, int disp);

  ///////////////////////////////////////////////////////////////
  // Extra methods added by derived
  ///////////////////////////////////////////////////////////////
  void DerivInternal(StencilImpl &st, DoubledGaugeField &U, GaugeField &mat, const FermionField &A, const FermionField &B, int dag);

  void DhopInternal(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, const FermionField &in, FermionField &out, int dag);

  // Constructor
  StaggeredFermion(GaugeField &_U, GridCartesian &Fgrid,
			   GridRedBlackCartesian &Hgrid, RealD _mass,
			   const ImplParams &p = ImplParams());

  // DoubleStore impl dependent
  void ImportGauge(const GaugeField &_Uthin);

  ///////////////////////////////////////////////////////////////
  // Data members require to support the functionality
  ///////////////////////////////////////////////////////////////

  //    protected:
 public:
  // any other parameters of action ???

  RealD mass;

  GridBase *_grid;
  GridBase *_cbgrid;

  // Defines the stencils for even and odd
  StencilImpl Stencil;
  StencilImpl StencilEven;
  StencilImpl StencilOdd;

  // Copy of the gauge field , with even and odd subsets
  DoubledGaugeField Umu;
  DoubledGaugeField UmuEven;
  DoubledGaugeField UmuOdd;

  LebesgueOrder Lebesgue;
  LebesgueOrder LebesgueEvenOdd;
};

typedef StaggeredFermion<StaggeredImplF> StaggeredFermionF;
typedef StaggeredFermion<StaggeredImplD> StaggeredFermionD;


}
}
#endif
