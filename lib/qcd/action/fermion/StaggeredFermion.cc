/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/fermion/StaggeredFermion.cc

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
#include <Grid.h>

namespace Grid {
namespace QCD {

const std::vector<int> 
StaggeredFermionStatic::directions({0, 1, 2, 3, 0, 1, 2, 3});
const std::vector<int> 
StaggeredFermionStatic::displacements({1, 1, 1, 1, -1, -1, -1, -1});

/////////////////////////////////
// Constructor and gauge import
/////////////////////////////////


template <class Impl>
StaggeredFermion<Impl>::StaggeredFermion(GaugeField &_Umu, GridCartesian &Fgrid, GridRedBlackCartesian &Hgrid,
							 RealD _mass,
							 const ImplParams &p)
    : Kernels(p),
      _grid(&Fgrid),
      _cbgrid(&Hgrid),
      Stencil(&Fgrid, npoint, Even, directions, displacements),
      StencilEven(&Hgrid, npoint, Even, directions, displacements),  // source is Even
      StencilOdd(&Hgrid, npoint, Odd, directions, displacements),  // source is Odd
      mass(_mass),
      Lebesgue(_grid),
      LebesgueEvenOdd(_cbgrid),
      Umu(&Fgrid),
      UmuEven(&Hgrid),
      UmuOdd(&Hgrid),
      _tmp(&Hgrid)
{
    ImportGauge(_Umu);
}

template <class Impl>
void StaggeredFermion<Impl>::ImportGauge(const GaugeField &_Umu)
{
  GaugeLinkField U(GaugeGrid());

  Impl::DoubleStore(GaugeGrid(), Umu, _Umu );

  ////////////////////////////////////////////////////////
  // 0.5 ( U p(x+mu) - Udag(x-mu) p(x-mu) ) 
  ////////////////////////////////////////////////////////
  for (int mu = 0; mu < Nd; mu++) {

    U = PeekIndex<LorentzIndex>(Umu, mu);
    PokeIndex<LorentzIndex>(Umu, U*( 0.5), mu );
    
    U = PeekIndex<LorentzIndex>(Umu, mu+4);
    PokeIndex<LorentzIndex>(Umu, U*(-0.5), mu+4);
  }

  pickCheckerboard(Even, UmuEven, Umu);
  pickCheckerboard(Odd,  UmuOdd , Umu);
}

/////////////////////////////
// Implement the interface
/////////////////////////////

template <class Impl>
RealD StaggeredFermion<Impl>::M(const FermionField &in, FermionField &out) {
  out.checkerboard = in.checkerboard;
  Dhop(in, out, DaggerNo);
  return axpy_norm(out, mass, in, out);
}

template <class Impl>
RealD StaggeredFermion<Impl>::Mdag(const FermionField &in, FermionField &out) {
  out.checkerboard = in.checkerboard;
  Dhop(in, out, DaggerYes);
  return axpy_norm(out, mass, in, out);
}

template <class Impl>
void StaggeredFermion<Impl>::Meooe(const FermionField &in, FermionField &out) {
  if (in.checkerboard == Odd) {
    DhopEO(in, out, DaggerNo);
  } else {
    DhopOE(in, out, DaggerNo);
  }
}
template <class Impl>
void StaggeredFermion<Impl>::MeooeDag(const FermionField &in, FermionField &out) {
  if (in.checkerboard == Odd) {
    DhopEO(in, out, DaggerYes);
  } else {
    DhopOE(in, out, DaggerYes);
  }
}

template <class Impl>
void StaggeredFermion<Impl>::Mooee(const FermionField &in, FermionField &out) {
  out.checkerboard = in.checkerboard;
  typename FermionField::scalar_type scal(mass);
  out = scal * in;
}

template <class Impl>
void StaggeredFermion<Impl>::MooeeDag(const FermionField &in, FermionField &out) {
  out.checkerboard = in.checkerboard;
  Mooee(in, out);
}

// Mooee is disabled in the perturbative case
template <>
void WilsonFermion<QCDpt::PStaggeredSmellImplF>::Mooee(const FermionField &in, FermionField &out) {
  assert(0);
}
template <>
void WilsonFermion<QCDpt::PStaggeredSmellImplD>::Mooee(const FermionField &in, FermionField &out) {
  assert(0);
}
template <>
void WilsonFermion<QCDpt::StaggeredSmellImplF>::Mooee(const FermionField &in, FermionField &out) {
  assert(0);
}
template <>
void WilsonFermion<QCDpt::StaggeredSmellImplD>::Mooee(const FermionField &in, FermionField &out) {
  assert(0);
}

template <class Impl>
void StaggeredFermion<Impl>::MooeeInv(const FermionField &in, FermionField &out) {
  out.checkerboard = in.checkerboard;
  out = (1.0 / (mass)) * in;
}

// MooeeInv is disabled in the perturbative case
template <>
void WilsonFermion<QCDpt::PStaggeredSmellImplF>::MooeeInv(const FermionField &in, FermionField &out) {
  assert(0);
}
template <>
void WilsonFermion<QCDpt::PStaggeredSmellImplD>::MooeeInv(const FermionField &in, FermionField &out) {
  assert(0);
}
template <>
void WilsonFermion<QCDpt::StaggeredSmellImplF>::MooeeInv(const FermionField &in, FermionField &out) {
  assert(0);
}
template <>
void WilsonFermion<QCDpt::StaggeredSmellImplD>::MooeeInv(const FermionField &in, FermionField &out) {
  assert(0);
}

template <class Impl>
void StaggeredFermion<Impl>::MooeeInvDag(const FermionField &in,
                                      FermionField &out) {
  out.checkerboard = in.checkerboard;
  MooeeInv(in, out);
}

///////////////////////////////////
// Internal
///////////////////////////////////

template <class Impl>
void StaggeredFermion<Impl>::DerivInternal(StencilImpl &st, DoubledGaugeField &U,
						   GaugeField & mat,
						   const FermionField &A, const FermionField &B, int dag) {
  assert((dag == DaggerNo) || (dag == DaggerYes));

  Compressor compressor;

  FermionField Btilde(B._grid);
  FermionField Atilde(B._grid);

  if(dag) Atilde = -A;
  else Atilde = A;

  st.HaloExchange(B, compressor);

  for (int mu = 0; mu < Nd; mu++) {

    ////////////////////////
    // Call the single hop
    ////////////////////////
    
    parallel_for (int sss = 0; sss < B._grid->oSites(); sss++) {
      Kernels::DhopDir(st, U, st.CommBuf(), sss, sss, B, Btilde, mu);
    }

      Impl::InsertForce4D(mat, Btilde, Atilde, mu);
  }
}

template <class Impl>
void StaggeredFermion<Impl>::DhopDeriv(GaugeField &mat, const FermionField &U, const FermionField &V, int dag) {

  conformable(U._grid, _grid);
  conformable(U._grid, V._grid);
  conformable(U._grid, mat._grid);

  mat.checkerboard = U.checkerboard;

  DerivInternal(Stencil, Umu, mat, U, V, dag);
}

template <class Impl>
void StaggeredFermion<Impl>::DhopDerivOE(GaugeField &mat, const FermionField &U, const FermionField &V, int dag) {

  conformable(U._grid, _cbgrid);
  conformable(U._grid, V._grid);
  conformable(U._grid, mat._grid);

  assert(V.checkerboard == Even);
  assert(U.checkerboard == Odd);
  mat.checkerboard = Odd;

  DerivInternal(StencilEven, UmuOdd, mat, U, V, dag);
}

template <class Impl>
void StaggeredFermion<Impl>::DhopDerivEO(GaugeField &mat, const FermionField &U, const FermionField &V, int dag) {

  conformable(U._grid, _cbgrid);
  conformable(U._grid, V._grid);
  conformable(U._grid, mat._grid);

  assert(V.checkerboard == Odd);
  assert(U.checkerboard == Even);
  mat.checkerboard = Even;

  DerivInternal(StencilOdd, UmuEven, mat, U, V, dag);
}

template <class Impl>
void StaggeredFermion<Impl>::Dhop(const FermionField &in, FermionField &out, int dag) {
  conformable(in._grid, _grid);  // verifies full grid
  conformable(in._grid, out._grid);

  out.checkerboard = in.checkerboard;
  DhopInternal(Stencil, Lebesgue, Umu, in, out, dag);
}

template <class Impl>
void StaggeredFermion<Impl>::DhopOE(const FermionField &in, FermionField &out, int dag) {
  conformable(in._grid, _cbgrid);    // verifies half grid
  conformable(in._grid, out._grid);  // drops the cb check

  assert(in.checkerboard == Even);
  out.checkerboard = Odd;

  DhopInternal(StencilEven, LebesgueEvenOdd, UmuOdd, in, out, dag);
}

template <class Impl>
void StaggeredFermion<Impl>::DhopEO(const FermionField &in, FermionField &out, int dag) {
  conformable(in._grid, _cbgrid);    // verifies half grid
  conformable(in._grid, out._grid);  // drops the cb check

  assert(in.checkerboard == Odd);
  out.checkerboard = Even;

  DhopInternal(StencilOdd, LebesgueEvenOdd, UmuEven, in, out, dag);
}

template <class Impl>
void StaggeredFermion<Impl>::Mdir(const FermionField &in, FermionField &out, int dir, int disp) {
  DhopDir(in, out, dir, disp);
}

template <class Impl>
void StaggeredFermion<Impl>::DhopDir(const FermionField &in, FermionField &out, int dir, int disp) {

  Compressor compressor;
  Stencil.HaloExchange(in, compressor);

  PARALLEL_FOR_LOOP
  for (int sss = 0; sss < in._grid->oSites(); sss++) {
    Kernels::DhopDir(Stencil, Umu, Stencil.CommBuf(), sss, sss, in, out, dir);
  }
};

template <class Impl>
void StaggeredFermion<Impl>::DhopInternal(StencilImpl &st, LebesgueOrder &lo,
						  DoubledGaugeField &U,
						  const FermionField &in,
						  FermionField &out, int dag) {
  assert((dag == DaggerNo) || (dag == DaggerYes));

  Compressor compressor;
  st.HaloExchange(in, compressor);

  if (dag == DaggerYes) {
    PARALLEL_FOR_LOOP
    for (int sss = 0; sss < in._grid->oSites(); sss++) {
      Kernels::DhopSiteDag(st, lo, U, st.CommBuf(), 1, sss, in, out);
    }
  } else {
    PARALLEL_FOR_LOOP
    for (int sss = 0; sss < in._grid->oSites(); sss++) {
      Kernels::DhopSite(st, lo, U, st.CommBuf(), 1, sss, in, out);
    }
  }
};

template class StaggeredFermion<QCDpt::PStaggeredSmellImplF>;
template class StaggeredFermion<QCDpt::PStaggeredSmellImplD>;
template class StaggeredFermion<QCDpt::StaggeredSmellImplF>;
template class StaggeredFermion<QCDpt::StaggeredSmellImplD>;
template class StaggeredFermion<QCDpt::NaiveStaggeredImplF>;
template class StaggeredFermion<QCDpt::NaiveStaggeredImplD>;

}}
