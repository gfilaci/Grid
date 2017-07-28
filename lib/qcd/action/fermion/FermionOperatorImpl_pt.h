/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/fermion/FermionOperatorImpl_pt.h

Copyright (C) 2015

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>
Author: paboyle <paboyle@ph.ed.ac.uk>
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
#ifndef GRID_QCD_FERMION_OPERATOR_IMPL_PT_H
#define GRID_QCD_FERMION_OPERATOR_IMPL_PT_H

namespace Grid {
namespace QCD {
namespace QCDpt {
  
#define INHERIT_FIMPL_TYPES(Impl)\
  typedef typename Impl::FermionField           FermionField;		\
  typedef typename Impl::PropagatorField     PropagatorField;		\
  typedef typename Impl::DoubledGaugeField DoubledGaugeField;		\
  typedef typename Impl::SiteSpinor               SiteSpinor;		\
  typedef typename Impl::SitePropagator       SitePropagator;		\
  typedef typename Impl::SiteHalfSpinor       SiteHalfSpinor;		\
  typedef typename Impl::Compressor               Compressor;		\
  typedef typename Impl::StencilImpl             StencilImpl;		\
  typedef typename Impl::ImplParams               ImplParams;	        \
  typedef typename Impl::Coeff_t                     Coeff_t;           \
  
#define INHERIT_IMPL_TYPES(Base) \
  INHERIT_GIMPL_TYPES(Base)      \
  INHERIT_FIMPL_TYPES(Base)

  /////////////////////////////////////////////////////////////////////////////
  // Perturbative fermion with smell
  /////////////////////////////////////////////////////////////////////////////
  template <class S, class Representation = FundamentalRepresentation,class Options = CoeffReal >
  class PWilsonSmellImpl : public TwistedGaugeImpl<GaugeImplTypes_pt<S, Representation::Dimension > > {
    public:

    static const int Dimension = Representation::Dimension;
    static const bool LsVectorised=false;
    static const int Nhcs = Options::Nhcs;

    typedef TwistedGaugeImpl<GaugeImplTypes_pt<S, Dimension > > Gimpl;
    INHERIT_GIMPL_TYPES(Gimpl);
      
    //Necessary?
    constexpr bool is_fundamental() const{return Dimension == Nc ? 1 : 0;}//$//
    
    typedef typename Options::_Coeff_t Coeff_t;
    typedef typename Options::template PrecisionMapper<Simd>::LowerPrecVector SimdL;
      
      
    template <typename vtype> using iImplSpinor            = iScalar<iVector<iPert<iMatrix<vtype, Dimension>, Np>, Ns> >;
    template <typename vtype> using iImplPropagator        = iScalar<iMatrix<iPert<iMatrix<vtype, Dimension>, Np>, Ns> >;
    template <typename vtype> using iImplHalfSpinor        = iScalar<iVector<iPert<iMatrix<vtype, Dimension>, Np>, Nhs> >;
    template <typename vtype> using iImplHalfCommSpinor    = iScalar<iVector<iPert<iMatrix<vtype, Dimension>, Np>, Nhcs> >;
    template <typename vtype> using iImplDoubledGaugeField = iVector<iScalar<iPert<iMatrix<vtype, Dimension>, Np> >, Nds>;
      
//    template <typename vtype> using iImplSpinor            = iScalar<iVector<iVector<vtype, Dimension>, Ns> >;
//    template <typename vtype> using iImplPropagator        = iScalar<iMatrix<iMatrix<vtype, Dimension>, Ns> >;
//    template <typename vtype> using iImplHalfSpinor        = iScalar<iVector<iVector<vtype, Dimension>, Nhs> >;
//    template <typename vtype> using iImplHalfCommSpinor    = iScalar<iVector<iVector<vtype, Dimension>, Nhcs> >;
//    template <typename vtype> using iImplDoubledGaugeField = iVector<iScalar<iMatrix<vtype, Dimension> >, Nds>;
    
    typedef iImplSpinor<Simd>            SiteSpinor;
    typedef iImplPropagator<Simd>        SitePropagator;
    typedef iImplHalfSpinor<Simd>        SiteHalfSpinor;
    typedef iImplHalfCommSpinor<SimdL>   SiteHalfCommSpinor;
    typedef iImplDoubledGaugeField<Simd> SiteDoubledGaugeField;
    
    typedef Lattice<SiteSpinor>            FermionField;
    typedef Lattice<SitePropagator>        PropagatorField;
    typedef Lattice<SiteDoubledGaugeField> DoubledGaugeField;
    
    typedef WilsonCompressor<SiteHalfCommSpinor,SiteHalfSpinor, SiteSpinor> Compressor;
    typedef WilsonImplParams ImplParams;
    typedef WilsonStencil<SiteSpinor, SiteHalfSpinor> StencilImpl;
    
    ImplParams Params;
    
    PWilsonSmellImpl(const ImplParams &p = ImplParams()) : Params(p){
      assert(Params.boundary_phases.size() == Nd);
    };
      
    bool overlapCommsCompute(void) { return Params.overlapCommsCompute; };
      
    inline void multLink(SiteHalfSpinor &phi,
                         const SiteDoubledGaugeField &U,
                         const SiteHalfSpinor &chi,
                         int mu,
                         StencilEntry *SE,
                         StencilImpl &St) {
      mult(&phi(), &U(mu), &chi());
    }
      
    template <class ref>
    inline void loadLinkElement(Simd &reg, ref &memory) {
      reg = memory;
    }
      
    inline void DoubleStore(GridBase *GaugeGrid,
                            DoubledGaugeField &Uds,
                            const GaugeField &Umu) 
    {
      typedef typename Simd::scalar_type scalar_type;

      conformable(Uds._grid, GaugeGrid);
      conformable(Umu._grid, GaugeGrid);

      GaugeLinkField U(GaugeGrid);
      GaugeLinkField tmp(GaugeGrid);

      Lattice<iScalar<vInteger> > coor(GaugeGrid);
      for (int mu = 0; mu < Nd; mu++) {

	      auto pha = Params.boundary_phases[mu];
	      scalar_type phase( real(pha),imag(pha) );

        int Lmu = GaugeGrid->GlobalDimensions()[mu] - 1;

        LatticeCoordinate(coor, mu);

        U = PeekIndex<LorentzIndex>(Umu, mu);
        tmp = where(coor == Lmu, phase * U, U);
        PokeIndex<LorentzIndex>(Uds, tmp, mu);

        U = adj(Cshift(U, mu, -1));
        U = where(coor == 0, conjugate(phase) * U, U); 
        PokeIndex<LorentzIndex>(Uds, U, mu + 4);
      }
    }

    inline void InsertForce4D(GaugeField &mat, FermionField &Btilde, FermionField &A,int mu){
      GaugeLinkField link(mat._grid);
      link = TraceIndex<SpinIndex>(outerProduct(Btilde,A)); 
      PokeIndex<LorentzIndex>(mat,link,mu);
    }   
      
    inline void InsertForce5D(GaugeField &mat, FermionField &Btilde, FermionField &Atilde,int mu){
      
      int Ls=Btilde._grid->_fdimensions[0];
      GaugeLinkField tmp(mat._grid);
      tmp = zero;
      
      parallel_for(int sss=0;sss<tmp._grid->oSites();sss++){
	    int sU=sss;
	    for(int s=0;s<Ls;s++){
	        int sF = s+Ls*sU;
	        tmp[sU] = tmp[sU]+ traceIndex<SpinIndex>(outerProduct(Btilde[sF],Atilde[sF])); // ordering here
	    }
      }
      PokeIndex<LorentzIndex>(mat,tmp,mu);
      
    }
  };

typedef PWilsonSmellImpl<vComplex,  FundamentalRepresentation, CoeffReal > PWilsonSmellImplR;  // Real.. whichever prec
typedef PWilsonSmellImpl<vComplexF, FundamentalRepresentation, CoeffReal > PWilsonSmellImplF;  // Float
typedef PWilsonSmellImpl<vComplexD, FundamentalRepresentation, CoeffReal > PWilsonSmellImplD;  // Double

}}}

#endif
