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

  /////////////////////////////////////////////////////////////////////////////
  // Fermion with smell
  /////////////////////////////////////////////////////////////////////////////
  template <class S, class Representation = FundamentalRepresentation,class Options = CoeffReal >
  class WilsonSmellImpl : public TwistedGaugeImpl<GaugeImplTypes_ptscalar<S, Representation::Dimension > > {
    public:

    static const int Dimension = Representation::Dimension;
    static const bool isFundamental = Representation::isFundamental;
    static const bool LsVectorised=false;
    static const int Nhcs = Options::Nhcs;

    typedef TwistedGaugeImpl<GaugeImplTypes_ptscalar<S, Dimension > > Gimpl;
    INHERIT_GIMPL_TYPES(Gimpl);
      
    //Necessary?
    constexpr bool is_fundamental() const{return Dimension == Nc ? 1 : 0;}
    
    typedef typename Options::_Coeff_t Coeff_t;
    typedef typename Options::template PrecisionMapper<Simd>::LowerPrecVector SimdL;
      
      
    template <typename vtype> using iImplSpinor            = iScalar<iVector<iScalar<iMatrix<vtype, Dimension> >, Ns> >;
    template <typename vtype> using iImplPropagator        = iScalar<iMatrix<iScalar<iMatrix<vtype, Dimension> >, Ns> >;
    template <typename vtype> using iImplHalfSpinor        = iScalar<iVector<iScalar<iMatrix<vtype, Dimension> >, Nhs> >;
    template <typename vtype> using iImplHalfCommSpinor    = iScalar<iVector<iScalar<iMatrix<vtype, Dimension> >, Nhcs> >;
    template <typename vtype> using iImplDoubledGaugeField = iVector<iScalar<iScalar<iMatrix<vtype, Dimension> > >, Nds>;
    
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
    
    typedef Lattice<QCDpt::iSinglet<typename FermionField::vector_type> > LatticeSinglet;
    
    ImplParams Params;
    
    WilsonSmellImpl(const ImplParams &p = ImplParams()) : Params(p){
      assert(Params.boundary_phases.size() == Nd);
    };
      
    bool overlapCommsCompute(void) { return Params.overlapCommsCompute; };
      
    inline void multLink(SiteHalfSpinor &phi,
                         const SiteDoubledGaugeField &U,
                         const SiteHalfSpinor &chi,
                         int mu,
                         StencilEntry *SE,
                         StencilImpl &St) {
        
        typedef SiteHalfSpinor vobj;
        typedef typename SiteHalfSpinor::scalar_object sobj;
        
        vobj vtmp;
        
        GridBase *grid = St._grid;
        
        const int Nsimd = grid->Nsimd();
        
        int direction = St._directions[mu];
        int distance = St._distances[mu];
        int ptype = St._permute_type[mu];
        int sl = St._grid->_simd_layout[direction];
        
        // Fixme X.Y.Z.T hardcode in stencil
        int mmu = mu % Nd;
        
        // assert our assumptions
        assert((distance == 1) || (distance == -1));  // nearest neighbour stencil hard code
        assert((sl == 1) || (sl == 2));
        
        std::vector<int> icoor;
        
        if ( SE->_around_the_world && istwisted(mmu) ) {
            
            if ( sl == 2 ) {
                
                std::vector<sobj> vals(Nsimd);
                
                extract(chi,vals);
                for(int s=0;s<Nsimd;s++){
                    
                    grid->iCoorFromIindex(icoor,s);
                    
                    assert((icoor[direction]==0)||(icoor[direction]==1));
                    
                    int permute_lane;
                    if ( distance == 1) {
                        permute_lane = icoor[direction]?1:0;
                    } else {
                        permute_lane = icoor[direction]?0:1;
                    }
                    
                    if ( permute_lane ) {
                        // distance = +1  -->  (U*omega) (psi*omegadag)
                        // distance = -1  -->  (omegadag*U) (psi*omega)
                        if(distance == 1) vals[s] = vals[s]*Gimpl::twist.adjomega[mmu];
                        else vals[s] = vals[s]*Gimpl::twist.omega[mmu];
                    }
                }
                merge(vtmp,vals);
                
            } else {
                // distance = +1  -->  (U*omega) (psi*omegadag)
                // distance = -1  -->  (omegadag*U) (psi*omega)
                if(distance == 1) vtmp = chi*Gimpl::twist.adjomega[mmu];
                else vtmp = chi*Gimpl::twist.omega[mmu];
            }
            mult(&phi(), &U(mu), &vtmp());
            
        } else { 
            mult(&phi(), &U(mu), &chi());
        }
        
    }
    
    inline void multLinkProp(SitePropagator &phi, const SiteDoubledGaugeField &U,
                             const SitePropagator &chi, int mu)
    {
        assert(0);
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
        if(istwisted(mu)) tmp = where(coor == Lmu, phase * U * Gimpl::twist.omega[mu], U);
        else tmp = where(coor == Lmu, phase * U, U);
        PokeIndex<LorentzIndex>(Uds, tmp, mu);

        U = adj(Cshift(U, mu, -1));
        if(istwisted(mu)) U = where(coor == 0, conjugate(phase) * Gimpl::twist.adjomega[mu] * U, U);
        else U = where(coor == 0, conjugate(phase) * U, U);
        PokeIndex<LorentzIndex>(Uds, U, mu + 4);
      }
    }

    inline void InsertForce4D(GaugeField &mat, FermionField &Btilde, FermionField &A,int mu){
        GaugeLinkField link(mat._grid);
        // matrix product automatically performs the sum over smells
        // and gives a matrix in colour space.
        // Need to resolve nesting because mult(iVector,iVector) is not defined...
        FermionField tmp = adj(A);
        link = zero;
        parallel_for(int ss=0;ss<mat._grid->oSites();ss++){
            for (int alpha=0; alpha<Ns; alpha++)
                link._odata[ss]._internal._internal += Btilde._odata[ss]._internal._internal[alpha] * tmp._odata[ss]._internal._internal[alpha];
        }
        PokeIndex<LorentzIndex>(mat,link,mu);
    }
      
    inline void InsertForce5D(GaugeField &mat, FermionField &Btilde, FermionField &Atilde,int mu){
        assert(0);
    }
  };
  
  
  /////////////////////////////////////////////////////////////////////////////
  // Perturbative fermion with smell
  /////////////////////////////////////////////////////////////////////////////
  template <class S, class Representation = FundamentalRepresentation,class Options = CoeffReal >
  class PWilsonSmellImpl : public TwistedGaugeImpl<GaugeImplTypes_pt<S, Representation::Dimension > > {
    public:
    // single order implementation
    typedef typename QCDpt::WilsonSmellImpl<S,Representation,CoeffReal> SOimpl;
    
    static const int Dimension = Representation::Dimension;
    static const bool isFundamental = Representation::isFundamental;
    static const bool LsVectorised=false;
    static const int Nhcs = Options::Nhcs;

    typedef TwistedGaugeImpl<GaugeImplTypes_pt<S, Dimension > > Gimpl;
    INHERIT_GIMPL_TYPES(Gimpl);
      
    //Necessary?
    constexpr bool is_fundamental() const{return Dimension == Nc ? 1 : 0;}
    
    typedef typename Options::_Coeff_t Coeff_t;
    typedef typename Options::template PrecisionMapper<Simd>::LowerPrecVector SimdL;
      
      
    template <typename vtype> using iImplSpinor            = iScalar<iVector<iPert<iMatrix<vtype, Dimension>, Np>, Ns> >;
    template <typename vtype> using iImplPropagator        = iScalar<iMatrix<iPert<iMatrix<vtype, Dimension>, Np>, Ns> >;
    template <typename vtype> using iImplHalfSpinor        = iScalar<iVector<iPert<iMatrix<vtype, Dimension>, Np>, Nhs> >;
    template <typename vtype> using iImplHalfCommSpinor    = iScalar<iVector<iPert<iMatrix<vtype, Dimension>, Np>, Nhcs> >;
    template <typename vtype> using iImplDoubledGaugeField = iVector<iScalar<iPert<iMatrix<vtype, Dimension>, Np> >, Nds>;
    
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
    
    typedef Lattice<QCDpt::iSinglet<typename FermionField::vector_type> > LatticeSinglet;
    
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
        
        typedef SiteHalfSpinor vobj;
        typedef typename SiteHalfSpinor::scalar_object sobj;
        
        vobj vtmp;
        
        GridBase *grid = St._grid;
        
        const int Nsimd = grid->Nsimd();
        
        int direction = St._directions[mu];
        int distance = St._distances[mu];
        int ptype = St._permute_type[mu];
        int sl = St._grid->_simd_layout[direction];
        
        // Fixme X.Y.Z.T hardcode in stencil
        int mmu = mu % Nd;
        
        // assert our assumptions
        assert((distance == 1) || (distance == -1));  // nearest neighbour stencil hard code
        assert((sl == 1) || (sl == 2));
        
        std::vector<int> icoor;
        
        if ( SE->_around_the_world && istwisted(mmu) ) {
            
            if ( sl == 2 ) {
                
                std::vector<sobj> vals(Nsimd);
                
                extract(chi,vals);
                for(int s=0;s<Nsimd;s++){
                    
                    grid->iCoorFromIindex(icoor,s);
                    
                    assert((icoor[direction]==0)||(icoor[direction]==1));
                    
                    int permute_lane;
                    if ( distance == 1) {
                        permute_lane = icoor[direction]?1:0;
                    } else {
                        permute_lane = icoor[direction]?0:1;
                    }
                    
                    if ( permute_lane ) {
                        // distance = +1  -->  (U*omega) (psi*omegadag)
                        // distance = -1  -->  (omegadag*U) (psi*omega)
                        if(distance == 1) vals[s] = vals[s]*Gimpl::twist.adjomega[mmu];
                        else vals[s] = vals[s]*Gimpl::twist.omega[mmu];
                    }
                }
                merge(vtmp,vals);
                
            } else {
                // distance = +1  -->  (U*omega) (psi*omegadag)
                // distance = -1  -->  (omegadag*U) (psi*omega)
                if(distance == 1) vtmp = chi*Gimpl::twist.adjomega[mmu];
                else vtmp = chi*Gimpl::twist.omega[mmu];
            }
            mult(&phi(), &U(mu), &vtmp());
            
        } else { 
            mult(&phi(), &U(mu), &chi());
        }
        
    }
      
    inline void multLinkProp(SitePropagator &phi, const SiteDoubledGaugeField &U,
                             const SitePropagator &chi, int mu)
    {
          assert(0);
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
        if(istwisted(mu)) tmp = where(coor == Lmu, phase * U * Gimpl::twist.omega[mu], U);
        else tmp = where(coor == Lmu, phase * U, U);
        PokeIndex<LorentzIndex>(Uds, tmp, mu);

        U = adj(Cshift(U, mu, -1));
        if(istwisted(mu)) U = where(coor == 0, conjugate(phase) * Gimpl::twist.adjomega[mu] * U, U);
        else U = where(coor == 0, conjugate(phase) * U, U);
        PokeIndex<LorentzIndex>(Uds, U, mu + 4);
      }
    }

    inline void InsertForce4D(GaugeField &mat, FermionField &Btilde, FermionField &A,int mu){
        assert(0);
    }
      
    inline void InsertForce5D(GaugeField &mat, FermionField &Btilde, FermionField &Atilde,int mu){
        assert(0);
    }
  };


/////////////////////////////////////////////////////////////////////////////
// Staggered fermion with smell
/////////////////////////////////////////////////////////////////////////////
template <class S, class Representation = FundamentalRepresentation >
class StaggeredSmellImpl : public TwistedGaugeImpl<GaugeImplTypes_ptscalar<S, Representation::Dimension > > {

    public:

    typedef RealD  _Coeff_t ;
    static const int Dimension = Representation::Dimension;
    static const bool isFundamental = Representation::isFundamental;
    static const bool LsVectorised=false;
    typedef TwistedGaugeImpl<GaugeImplTypes_ptscalar<S, Dimension > > Gimpl;
      
    //Necessary?
    constexpr bool is_fundamental() const{return Dimension == Nc ? 1 : 0;}
    
    typedef _Coeff_t Coeff_t;

    INHERIT_GIMPL_TYPES(Gimpl);
      
    template <typename vtype> using iImplSpinor            = iScalar<iScalar<iScalar<iMatrix<vtype, Dimension> > > >;
    template <typename vtype> using iImplHalfSpinor        = iScalar<iScalar<iScalar<iMatrix<vtype, Dimension> > > >;
    template <typename vtype> using iImplDoubledGaugeField = iVector<iScalar<iScalar<iMatrix<vtype, Dimension> > >, Nds>;
    template <typename vtype> using iImplPropagator        = iScalar<iScalar<iScalar<iMatrix<vtype, Dimension> > > >;
    
    typedef iImplSpinor<Simd>            SiteSpinor;
    typedef iImplHalfSpinor<Simd>        SiteHalfSpinor;
    typedef iImplDoubledGaugeField<Simd> SiteDoubledGaugeField;
    typedef iImplPropagator<Simd>        SitePropagator;
    
    typedef Lattice<SiteSpinor>            FermionField;
    typedef Lattice<SiteDoubledGaugeField> DoubledGaugeField;
    typedef Lattice<SitePropagator> PropagatorField;
    
    typedef SimpleCompressor<SiteSpinor> Compressor;
    typedef WilsonImplParams ImplParams; // because boundary phases are needed
    typedef CartesianStencil<SiteSpinor, SiteSpinor> StencilImpl;
    
    typedef Lattice<QCDpt::iSinglet<typename FermionField::vector_type> > LatticeSinglet;
    
    ImplParams Params;
    
    StaggeredSmellImpl(const ImplParams &p = ImplParams()) : Params(p){};
      
    inline void multLink(SiteSpinor &phi,
			 const SiteDoubledGaugeField &U,
			 const SiteSpinor &chi,
			 int mu,
             StencilEntry *SE,
             StencilImpl &St) {
        
        typedef SiteSpinor vobj;
        typedef typename SiteSpinor::scalar_object sobj;
        
        vobj vtmp;
        
        GridBase *grid = St._grid;
        
        const int Nsimd = grid->Nsimd();
        
        int direction = St._directions[mu];
        int distance = St._distances[mu];
        int ptype = St._permute_type[mu];
        int sl = St._grid->_simd_layout[direction];
        
        // Fixme X.Y.Z.T hardcode in stencil
        int mmu = mu % Nd;
        
        // assert our assumptions
        assert((distance == 1) || (distance == -1));  // nearest neighbour stencil hard code
        assert((sl == 1) || (sl == 2));
        
        std::vector<int> icoor;
        
        if ( SE->_around_the_world && istwisted(mmu) ) {
            
            if ( sl == 2 ) {
                
                std::vector<sobj> vals(Nsimd);
                
                extract(chi,vals);
                for(int s=0;s<Nsimd;s++){
                    
                    grid->iCoorFromIindex(icoor,s);
                    
                    assert((icoor[direction]==0)||(icoor[direction]==1));
                    
                    int permute_lane;
                    if ( distance == 1) {
                        permute_lane = icoor[direction]?1:0;
                    } else {
                        permute_lane = icoor[direction]?0:1;
                    }
                    
                    if ( permute_lane ) {
                        // distance = +1  -->  (U*omega) (psi*omegadag)
                        // distance = -1  -->  (omegadag*U) (psi*omega)
                        if(distance == 1) vals[s] = vals[s]*Gimpl::twist.adjomega[mmu];
                        else vals[s] = vals[s]*Gimpl::twist.omega[mmu];
                    }
                }
                merge(vtmp,vals);
                
            } else {
                // distance = +1  -->  (U*omega) (psi*omegadag)
                // distance = -1  -->  (omegadag*U) (psi*omega)
                if(distance == 1) vtmp = chi*Gimpl::twist.adjomega[mmu];
                else vtmp = chi*Gimpl::twist.omega[mmu];
            }
            mult(&phi(), &U(mu), &vtmp());
            
        } else { 
            mult(&phi(), &U(mu), &chi());
        }
        
    }
    inline void multLinkAdd(SiteSpinor &phi,
			    const SiteDoubledGaugeField &U,
			    const SiteSpinor &chi,
			    int mu,
                StencilEntry *SE,
                StencilImpl &St) {
        
        typedef SiteSpinor vobj;
        typedef typename SiteSpinor::scalar_object sobj;
        
        vobj vtmp;
        
        GridBase *grid = St._grid;
        
        const int Nsimd = grid->Nsimd();
        
        int direction = St._directions[mu];
        int distance = St._distances[mu];
        int ptype = St._permute_type[mu];
        int sl = St._grid->_simd_layout[direction];
        
        // Fixme X.Y.Z.T hardcode in stencil
        int mmu = mu % Nd;
        
        // assert our assumptions
        assert((distance == 1) || (distance == -1));  // nearest neighbour stencil hard code
        assert((sl == 1) || (sl == 2));
        
        std::vector<int> icoor;
        
        if ( SE->_around_the_world && istwisted(mmu) ) {
            
            if ( sl == 2 ) {
                
                std::vector<sobj> vals(Nsimd);
                
                extract(chi,vals);
                for(int s=0;s<Nsimd;s++){
                    
                    grid->iCoorFromIindex(icoor,s);
                    
                    assert((icoor[direction]==0)||(icoor[direction]==1));
                    
                    int permute_lane;
                    if ( distance == 1) {
                        permute_lane = icoor[direction]?1:0;
                    } else {
                        permute_lane = icoor[direction]?0:1;
                    }
                    
                    if ( permute_lane ) {
                        // distance = +1  -->  (U*omega) (psi*omegadag)
                        // distance = -1  -->  (omegadag*U) (psi*omega)
                        if(distance == 1) vals[s] = vals[s]*Gimpl::twist.adjomega[mmu];
                        else vals[s] = vals[s]*Gimpl::twist.omega[mmu];
                    }
                }
                merge(vtmp,vals);
                
            } else {
                // distance = +1  -->  (U*omega) (psi*omegadag)
                // distance = -1  -->  (omegadag*U) (psi*omega)
                if(distance == 1) vtmp = chi*Gimpl::twist.adjomega[mmu];
                else vtmp = chi*Gimpl::twist.omega[mmu];
            }
            mac(&phi(), &U(mu), &vtmp());
            
        } else { 
            mac(&phi(), &U(mu), &chi());
        }
        
    }
      
    template <class ref>
    inline void loadLinkElement(Simd &reg, ref &memory) {
      reg = memory;
    }
      
    inline void DoubleStore(GridBase *GaugeGrid, DoubledGaugeField &Uds, const GaugeField &Uthin) {
      conformable(Uds._grid, GaugeGrid);
      conformable(Uthin._grid, GaugeGrid);
      GaugeLinkField U(GaugeGrid);
      GaugeLinkField Udag(GaugeGrid);
      Lattice<iScalar<vInteger> > coor(GaugeGrid);
      typedef typename Simd::scalar_type scalar_type;
      
      for (int mu = 0; mu < Nd; mu++) {

	// Staggered Phase.
	Lattice<iScalar<vInteger> > coor(GaugeGrid);
	Lattice<iScalar<vInteger> > x(GaugeGrid); LatticeCoordinate(x,0);
	Lattice<iScalar<vInteger> > y(GaugeGrid); LatticeCoordinate(y,1);
	Lattice<iScalar<vInteger> > z(GaugeGrid); LatticeCoordinate(z,2);
	Lattice<iScalar<vInteger> > t(GaugeGrid); LatticeCoordinate(t,3);

	Lattice<iScalar<vInteger> > lin_z(GaugeGrid); lin_z=x+y;
	Lattice<iScalar<vInteger> > lin_t(GaugeGrid); lin_t=x+y+z;

	ComplexField phases(GaugeGrid);	phases=1.0;

	if ( mu == 1 ) phases = where( mod(x    ,2)==(Integer)0, phases,-phases);
	if ( mu == 2 ) phases = where( mod(lin_z,2)==(Integer)0, phases,-phases);
	if ( mu == 3 ) phases = where( mod(lin_t,2)==(Integer)0, phases,-phases);

	U      = PeekIndex<LorentzIndex>(Uthin, mu);
	Udag   = adj( Cshift(U, mu, -1));

	U    = U    *phases;
	Udag = Udag *phases;


    // fermion phase and twist
    auto pha = Params.boundary_phases[mu];
    scalar_type phase( real(pha),imag(pha) );
    int Lmu = GaugeGrid->GlobalDimensions()[mu] - 1;
    LatticeCoordinate(coor, mu);
    if(istwisted(mu)){
        U = where(coor == Lmu, phase * U * Gimpl::twist.omega[mu], U);
        Udag = where(coor == 0, conjugate(phase) * Gimpl::twist.adjomega[mu] * Udag, Udag);
    }
    else{
        U = where(coor == Lmu, phase * U, U);
        Udag = where(coor == 0, conjugate(phase) * Udag, Udag);
    }

	PokeIndex<LorentzIndex>(Uds, U, mu);
	PokeIndex<LorentzIndex>(Uds, Udag, mu + 4);

      }
    }

    inline void InsertForce4D(GaugeField &mat, FermionField &Btilde, FermionField &A,int mu){
        GaugeLinkField link(mat._grid);
        // matrix product automatically performs the sum over smells
        // and gives a matrix in colour space.
        // Need to resolve nesting because mult(iVector,iVector) is not defined...
        FermionField tmp = adj(A);
        link = zero;
        parallel_for(int ss=0;ss<mat._grid->oSites();ss++){
                link._odata[ss]._internal._internal += Btilde._odata[ss]._internal._internal * tmp._odata[ss]._internal._internal;
        }
        PokeIndex<LorentzIndex>(mat,link,mu);
    }
      
    inline void InsertForce5D(GaugeField &mat, FermionField &Btilde, FermionField &Atilde,int mu){
      assert (0); 
      // Must never hit
    }
  };

  
/////////////////////////////////////////////////////////////////////////////
// Perturbative staggered fermion with smell
/////////////////////////////////////////////////////////////////////////////
template <class S, class Representation = FundamentalRepresentation >
class PStaggeredSmellImpl : public TwistedGaugeImpl<GaugeImplTypes_pt<S, Representation::Dimension > > {

    public:

    // single order implementation
    typedef typename QCDpt::StaggeredSmellImpl<S,Representation> SOimpl;
    
    typedef RealD  _Coeff_t ;
    static const int Dimension = Representation::Dimension;
    static const bool isFundamental = Representation::isFundamental;
    static const bool LsVectorised=false;
    typedef TwistedGaugeImpl<GaugeImplTypes_pt<S, Dimension > > Gimpl;
      
    //Necessary?
    constexpr bool is_fundamental() const{return Dimension == Nc ? 1 : 0;}
    
    typedef _Coeff_t Coeff_t;

    INHERIT_GIMPL_TYPES(Gimpl);
      
    template <typename vtype> using iImplSpinor            = iScalar<iScalar<iPert<iMatrix<vtype, Dimension>, Np> > >;
    template <typename vtype> using iImplHalfSpinor        = iScalar<iScalar<iPert<iVector<vtype, Dimension>, Np> > >;
    template <typename vtype> using iImplDoubledGaugeField = iVector<iScalar<iPert<iMatrix<vtype, Dimension>, Np> >, 2*Nds>;
    template <typename vtype> using iImplPropagator        = iScalar<iScalar<iPert<iMatrix<vtype, Dimension>, Np> > >;
    
    typedef iImplSpinor<Simd>            SiteSpinor;
    typedef iImplHalfSpinor<Simd>        SiteHalfSpinor;
    typedef iImplDoubledGaugeField<Simd> SiteDoubledGaugeField;
    typedef iImplPropagator<Simd>        SitePropagator;
    
    typedef Lattice<SiteSpinor>            FermionField;
    typedef Lattice<SiteDoubledGaugeField> DoubledGaugeField;
    typedef Lattice<SitePropagator> PropagatorField;
    
    typedef SimpleCompressor<SiteSpinor> Compressor;
    typedef WilsonImplParams ImplParams; // because boundary phases are needed
    typedef CartesianStencil<SiteSpinor, SiteSpinor> StencilImpl;
    
    typedef Lattice<QCDpt::iSinglet<typename FermionField::vector_type> > LatticeSinglet;
    
    ImplParams Params;
    
    PStaggeredSmellImpl(const ImplParams &p = ImplParams()) : Params(p){};
      
    inline void multLink(SiteSpinor &phi,
			 const SiteDoubledGaugeField &U,
			 const SiteSpinor &chi,
			 int mu,
             StencilEntry *SE,
                         StencilImpl &St) {
        
        typedef SiteSpinor vobj;
        typedef typename SiteSpinor::scalar_object sobj;
        
        vobj vtmp;
        
        GridBase *grid = St._grid;
        
        const int Nsimd = grid->Nsimd();
        
        int direction = St._directions[mu];
        int distance = St._distances[mu];
        int ptype = St._permute_type[mu];
        int sl = St._grid->_simd_layout[direction];
        
        // Fixme X.Y.Z.T hardcode in stencil
        int mmu = mu % Nd;
        
        // assert our assumptions
        assert((distance == 1) || (distance == -1));  // nearest neighbour stencil hard code
        assert((sl == 1) || (sl == 2));
        
        std::vector<int> icoor;
        
        if ( SE->_around_the_world && istwisted(mmu) ) {
            
            if ( sl == 2 ) {
                
                std::vector<sobj> vals(Nsimd);
                
                extract(chi,vals);
                for(int s=0;s<Nsimd;s++){
                    
                    grid->iCoorFromIindex(icoor,s);
                    
                    assert((icoor[direction]==0)||(icoor[direction]==1));
                    
                    int permute_lane;
                    if ( distance == 1) {
                        permute_lane = icoor[direction]?1:0;
                    } else {
                        permute_lane = icoor[direction]?0:1;
                    }
                    
                    if ( permute_lane ) {
                        // distance = +1  -->  (U*omega) (psi*omegadag)
                        // distance = -1  -->  (omegadag*U) (psi*omega)
                        if(distance == 1) vals[s] = vals[s]*Gimpl::twist.adjomega[mmu];
                        else vals[s] = vals[s]*Gimpl::twist.omega[mmu];
                    }
                }
                merge(vtmp,vals);
                
            } else {
                // distance = +1  -->  (U*omega) (psi*omegadag)
                // distance = -1  -->  (omegadag*U) (psi*omega)
                if(distance == 1) vtmp = chi*Gimpl::twist.adjomega[mmu];
                else vtmp = chi*Gimpl::twist.omega[mmu];
            }
            mult(&phi(), &U(mu), &vtmp());
            
        } else {
            mult(&phi(), &U(mu), &chi());
        }
        
    }
    inline void multLinkAdd(SiteSpinor &phi,
			    const SiteDoubledGaugeField &U,
			    const SiteSpinor &chi,
			    int mu,
                StencilEntry *SE,
                            StencilImpl &St) {
        
        typedef SiteSpinor vobj;
        typedef typename SiteSpinor::scalar_object sobj;
        
        vobj vtmp;
        
        GridBase *grid = St._grid;
        
        const int Nsimd = grid->Nsimd();
        
        int direction = St._directions[mu];
        int distance = St._distances[mu];
        int ptype = St._permute_type[mu];
        int sl = St._grid->_simd_layout[direction];
        
        // Fixme X.Y.Z.T hardcode in stencil
        int mmu = mu % Nd;
        
        // assert our assumptions
        assert((distance == 1) || (distance == -1));  // nearest neighbour stencil hard code
        assert((sl == 1) || (sl == 2));
        
        std::vector<int> icoor;
        
        if ( SE->_around_the_world && istwisted(mmu) ) {
            
            if ( sl == 2 ) {
                
                std::vector<sobj> vals(Nsimd);
                
                extract(chi,vals);
                for(int s=0;s<Nsimd;s++){
                    
                    grid->iCoorFromIindex(icoor,s);
                    
                    assert((icoor[direction]==0)||(icoor[direction]==1));
                    
                    int permute_lane;
                    if ( distance == 1) {
                        permute_lane = icoor[direction]?1:0;
                    } else {
                        permute_lane = icoor[direction]?0:1;
                    }
                    
                    if ( permute_lane ) {
                        // distance = +1  -->  (U*omega) (psi*omegadag)
                        // distance = -1  -->  (omegadag*U) (psi*omega)
                        if(distance == 1) vals[s] = vals[s]*Gimpl::twist.adjomega[mmu];
                        else vals[s] = vals[s]*Gimpl::twist.omega[mmu];
                    }
                }
                merge(vtmp,vals);
                
            } else {
                // distance = +1  -->  (U*omega) (psi*omegadag)
                // distance = -1  -->  (omegadag*U) (psi*omega)
                if(distance == 1) vtmp = chi*Gimpl::twist.adjomega[mmu];
                else vtmp = chi*Gimpl::twist.omega[mmu];
            }
            mac(&phi(), &U(mu), &vtmp());
            
        } else {
            mac(&phi(), &U(mu), &chi());
        }
        
    }
      
    template <class ref>
    inline void loadLinkElement(Simd &reg, ref &memory) {
      reg = memory;
    }
      
    inline void DoubleStore(GridBase *GaugeGrid, DoubledGaugeField &Uds, const GaugeField &Uthin) {
        conformable(Uds._grid, GaugeGrid);
        conformable(Uthin._grid, GaugeGrid);
        GaugeLinkField U(GaugeGrid);
        GaugeLinkField Udag(GaugeGrid);
        GaugeLinkField U_mirror(GaugeGrid);
        GaugeLinkField Udag_mirror(GaugeGrid);
        Lattice<iScalar<vInteger> > coor(GaugeGrid);
        typedef typename Simd::scalar_type scalar_type;
        
        for (int mu = 0; mu < Nd; mu++) {
            
            // Staggered Phase.
            Lattice<iScalar<vInteger> > coor(GaugeGrid);
            Lattice<iScalar<vInteger> > x(GaugeGrid); LatticeCoordinate(x,0);
            Lattice<iScalar<vInteger> > y(GaugeGrid); LatticeCoordinate(y,1);
            Lattice<iScalar<vInteger> > z(GaugeGrid); LatticeCoordinate(z,2);
            Lattice<iScalar<vInteger> > t(GaugeGrid); LatticeCoordinate(t,3);
            
            Lattice<iScalar<vInteger> > lin_z(GaugeGrid); lin_z=x+y;
            Lattice<iScalar<vInteger> > lin_t(GaugeGrid); lin_t=x+y+z;
            
            ComplexField phases(GaugeGrid);    phases=1.0;
            
            if ( mu == 1 ) phases = where( mod(x    ,2)==(Integer)0, phases,-phases);
            if ( mu == 2 ) phases = where( mod(lin_z,2)==(Integer)0, phases,-phases);
            if ( mu == 3 ) phases = where( mod(lin_t,2)==(Integer)0, phases,-phases);
            
            U      = PeekIndex<LorentzIndex>(Uthin, mu);
            Udag   = adj( Cshift(U, mu, -1));
            
            U_mirror    = adj(U);
            Udag_mirror = adj(Udag);
            U    = U    *phases;
            Udag = Udag *phases;
            
            
            // fermion phase and twist
            auto pha = Params.boundary_phases[mu];
            scalar_type phase( real(pha),imag(pha) );
            int Lmu = GaugeGrid->GlobalDimensions()[mu] - 1;
            LatticeCoordinate(coor, mu);
            if(istwisted(mu)){
                U        = where(coor == Lmu, phase * U * Gimpl::twist.omega[mu], U);
                U_mirror = where(coor == Lmu, Gimpl::twist.adjomega[mu] * U_mirror, U_mirror);
                Udag        = where(coor == 0, conjugate(phase) * Gimpl::twist.adjomega[mu] * Udag, Udag);
                Udag_mirror = where(coor == 0, Udag_mirror * Gimpl::twist.omega[mu], Udag_mirror);
            }
            else{
                U = where(coor == Lmu, phase * U, U);
                Udag = where(coor == 0, conjugate(phase) * Udag, Udag);
            }
            
            PokeIndex<LorentzIndex>(Uds, U,           mu     );
            PokeIndex<LorentzIndex>(Uds, Udag,        mu + 4 );
            PokeIndex<LorentzIndex>(Uds, U_mirror,    mu + 8 );
            PokeIndex<LorentzIndex>(Uds, Udag_mirror, mu + 12);
        }
    }

    inline void InsertForce4D(GaugeField &mat, FermionField &Btilde, FermionField &A,int mu){
      assert(0);
    }   
      
    inline void InsertForce5D(GaugeField &mat, FermionField &Btilde, FermionField &Atilde,int mu){
      assert (0); 
      // Must never hit
    }
  };
  
  /////////////////////////////////////////////////////////////////////////////
// Nonperturbative naive staggered fermion
/////////////////////////////////////////////////////////////////////////////
template <class S, class Representation = FundamentalRepresentation >
class NaiveStaggeredImpl : public PeriodicGaugeImpl<GaugeImplTypes<S, Representation::Dimension > > {

    public:
    
    typedef RealD  _Coeff_t ;
    static const int Dimension = Representation::Dimension;
    static const bool isFundamental = Representation::isFundamental;
    static const bool LsVectorised=false;
    typedef PeriodicGaugeImpl<GaugeImplTypes<S, Dimension > > Gimpl;
      
    //Necessary?
    constexpr bool is_fundamental() const{return Dimension == Nc ? 1 : 0;}
    
    typedef _Coeff_t Coeff_t;

    INHERIT_GIMPL_TYPES(Gimpl);
      
    template <typename vtype> using iImplSpinor            = iScalar<iScalar<iVector<vtype, Dimension> > >;
    template <typename vtype> using iImplHalfSpinor        = iScalar<iScalar<iVector<vtype, Dimension> > >;
    template <typename vtype> using iImplDoubledGaugeField = iVector<iScalar<iMatrix<vtype, Dimension> >, Nds>;
    template <typename vtype> using iImplPropagator        = iScalar<iScalar<iMatrix<vtype, Dimension> > >;
    
    typedef iImplSpinor<Simd>            SiteSpinor;
    typedef iImplHalfSpinor<Simd>        SiteHalfSpinor;
    typedef iImplDoubledGaugeField<Simd> SiteDoubledGaugeField;
    typedef iImplPropagator<Simd>        SitePropagator;
    
    typedef Lattice<SiteSpinor>            FermionField;
    typedef Lattice<SiteDoubledGaugeField> DoubledGaugeField;
    typedef Lattice<SitePropagator> PropagatorField;
    
    typedef SimpleCompressor<SiteSpinor> Compressor;
    typedef StaggeredImplParams ImplParams;
    typedef CartesianStencil<SiteSpinor, SiteSpinor> StencilImpl;
    
    typedef Lattice<iSinglet<typename FermionField::vector_type> > LatticeSinglet;
    
    ImplParams Params;
    
    NaiveStaggeredImpl(const ImplParams &p = ImplParams()) : Params(p){};
      
    inline void multLink(SiteSpinor &phi,
			 const SiteDoubledGaugeField &U,
			 const SiteSpinor &chi,
			 int mu,
             StencilEntry *SE,
             StencilImpl &St) {
        mult(&phi(), &U(mu), &chi());
    }
    
    inline void multLinkAdd(SiteSpinor &phi,
			    const SiteDoubledGaugeField &U,
			    const SiteSpinor &chi,
			    int mu,
                StencilEntry *SE,
                StencilImpl &St) {
        mac(&phi(), &U(mu), &chi());
    }
      
    template <class ref>
    inline void loadLinkElement(Simd &reg, ref &memory) {
      reg = memory;
    }
      
    inline void DoubleStore(GridBase *GaugeGrid, DoubledGaugeField &Uds, const GaugeField &Uthin) {
      conformable(Uds._grid, GaugeGrid);
      conformable(Uthin._grid, GaugeGrid);
      GaugeLinkField U(GaugeGrid);
      GaugeLinkField Udag(GaugeGrid);
      Lattice<iScalar<vInteger> > coor(GaugeGrid);
      typedef typename Simd::scalar_type scalar_type;
      
      for (int mu = 0; mu < Nd; mu++) {

	// Staggered Phase.
	Lattice<iScalar<vInteger> > coor(GaugeGrid);
	Lattice<iScalar<vInteger> > x(GaugeGrid); LatticeCoordinate(x,0);
	Lattice<iScalar<vInteger> > y(GaugeGrid); LatticeCoordinate(y,1);
	Lattice<iScalar<vInteger> > z(GaugeGrid); LatticeCoordinate(z,2);
	Lattice<iScalar<vInteger> > t(GaugeGrid); LatticeCoordinate(t,3);

	Lattice<iScalar<vInteger> > lin_z(GaugeGrid); lin_z=x+y;
	Lattice<iScalar<vInteger> > lin_t(GaugeGrid); lin_t=x+y+z;

	ComplexField phases(GaugeGrid);	phases=1.0;

	if ( mu == 1 ) phases = where( mod(x    ,2)==(Integer)0, phases,-phases);
	if ( mu == 2 ) phases = where( mod(lin_z,2)==(Integer)0, phases,-phases);
	if ( mu == 3 ) phases = where( mod(lin_t,2)==(Integer)0, phases,-phases);

	U      = PeekIndex<LorentzIndex>(Uthin, mu);
	Udag   = adj( Cshift(U, mu, -1));

	U    = U    *phases;
	Udag = Udag *phases;
    
	PokeIndex<LorentzIndex>(Uds, U, mu);
	PokeIndex<LorentzIndex>(Uds, Udag, mu + 4);

      }
    }

    inline void InsertForce4D(GaugeField &mat, FermionField &Btilde, FermionField &A,int mu){
        GaugeLinkField link(mat._grid);
        link = TraceIndex<SpinIndex>(outerProduct(Btilde,A));
        PokeIndex<LorentzIndex>(mat,link,mu);
    }
      
    inline void InsertForce5D(GaugeField &mat, FermionField &Btilde, FermionField &Atilde,int mu){
      assert (0); 
      // Must never hit
    }
  };
  
/////////////////////////////////////////////////////////////////////////////
// Staggered fermion in the adjoint representation
/////////////////////////////////////////////////////////////////////////////
template <class S, class Representation = FundamentalRepresentation >
class StaggeredAdjointImpl : public TwistedGaugeImpl<GaugeImplTypes_ptscalar<S, Representation::Dimension > > {
    
public:
    
    typedef RealD  _Coeff_t ;
    static const int Dimension = Representation::Dimension;
    static const bool isFundamental = Representation::isFundamental;
    static const bool LsVectorised=false;
    typedef TwistedGaugeImpl<GaugeImplTypes_ptscalar<S, Dimension > > Gimpl;
    
    //Necessary?
    constexpr bool is_fundamental() const{return Dimension == Nc ? 1 : 0;}
    
    typedef _Coeff_t Coeff_t;
    
    INHERIT_GIMPL_TYPES(Gimpl);
    
    template <typename vtype> using iImplSpinor            = iScalar<iScalar<iScalar<iMatrix<vtype, Dimension> > > >;
    template <typename vtype> using iImplHalfSpinor        = iScalar<iScalar<iScalar<iMatrix<vtype, Dimension> > > >;
    template <typename vtype> using iImplDoubledGaugeField = iVector<iScalar<iScalar<iMatrix<vtype, Dimension> > >, 2*Nds>;
    template <typename vtype> using iImplPropagator        = iScalar<iScalar<iScalar<iMatrix<vtype, Dimension> > > >;
    
    typedef iImplSpinor<Simd>            SiteSpinor;
    typedef iImplHalfSpinor<Simd>        SiteHalfSpinor;
    typedef iImplDoubledGaugeField<Simd> SiteDoubledGaugeField;
    typedef iImplPropagator<Simd>        SitePropagator;
    
    typedef Lattice<SiteSpinor>            FermionField;
    typedef Lattice<SiteDoubledGaugeField> DoubledGaugeField;
    typedef Lattice<SitePropagator> PropagatorField;
    
    typedef SimpleCompressor<SiteSpinor> Compressor;
    typedef WilsonImplParams ImplParams; // because boundary phases are needed
    typedef CartesianStencil<SiteSpinor, SiteSpinor> StencilImpl;
    
    typedef Lattice<QCDpt::iSinglet<typename FermionField::vector_type> > LatticeSinglet;
    
    ImplParams Params;
    
    StaggeredAdjointImpl(const ImplParams &p = ImplParams()) : Params(p){};
    
    inline void multLink(SiteSpinor &phi,
                         const SiteDoubledGaugeField &U,
                         const SiteSpinor &chi,
                         int mu,
                         StencilEntry *SE,
                         StencilImpl &St) {
        SiteSpinor phitmp;
        mult(&phitmp(), &U(mu), &chi());
        mult(&phi(), &phitmp(), &U(mu+Nds));
    }
    inline void multLinkAdd(SiteSpinor &phi,
                            const SiteDoubledGaugeField &U,
                            const SiteSpinor &chi,
                            int mu,
                            StencilEntry *SE,
                            StencilImpl &St) {
        SiteSpinor phitmp;
        mult(&phitmp(), &U(mu), &chi());
        mac(&phi(), &phitmp(), &U(mu+Nds));
    }
    
    template <class ref>
    inline void loadLinkElement(Simd &reg, ref &memory) {
        reg = memory;
    }
    
    inline void DoubleStore(GridBase *GaugeGrid, DoubledGaugeField &Uds, const GaugeField &Uthin) {
        conformable(Uds._grid, GaugeGrid);
        conformable(Uthin._grid, GaugeGrid);
        GaugeLinkField U(GaugeGrid);
        GaugeLinkField Udag(GaugeGrid);
        GaugeLinkField U_mirror(GaugeGrid);
        GaugeLinkField Udag_mirror(GaugeGrid);
        Lattice<iScalar<vInteger> > coor(GaugeGrid);
        typedef typename Simd::scalar_type scalar_type;
        
        for (int mu = 0; mu < Nd; mu++) {
            
            // Staggered Phase.
            Lattice<iScalar<vInteger> > coor(GaugeGrid);
            Lattice<iScalar<vInteger> > x(GaugeGrid); LatticeCoordinate(x,0);
            Lattice<iScalar<vInteger> > y(GaugeGrid); LatticeCoordinate(y,1);
            Lattice<iScalar<vInteger> > z(GaugeGrid); LatticeCoordinate(z,2);
            Lattice<iScalar<vInteger> > t(GaugeGrid); LatticeCoordinate(t,3);
            
            Lattice<iScalar<vInteger> > lin_z(GaugeGrid); lin_z=x+y;
            Lattice<iScalar<vInteger> > lin_t(GaugeGrid); lin_t=x+y+z;
            
            ComplexField phases(GaugeGrid);    phases=1.0;
            
            if ( mu == 1 ) phases = where( mod(x    ,2)==(Integer)0, phases,-phases);
            if ( mu == 2 ) phases = where( mod(lin_z,2)==(Integer)0, phases,-phases);
            if ( mu == 3 ) phases = where( mod(lin_t,2)==(Integer)0, phases,-phases);
            
            U      = PeekIndex<LorentzIndex>(Uthin, mu);
            Udag   = adj( Cshift(U, mu, -1));
            
            U_mirror    = adj(U);
            Udag_mirror = adj(Udag);
            U    = U    *phases;
            Udag = Udag *phases;
            
            
            // fermion phase and twist
            auto pha = Params.boundary_phases[mu];
            scalar_type phase( real(pha),imag(pha) );
            int Lmu = GaugeGrid->GlobalDimensions()[mu] - 1;
            LatticeCoordinate(coor, mu);
            if(istwisted(mu)){
                U        = where(coor == Lmu, phase * U * Gimpl::twist.omega[mu], U);
                U_mirror = where(coor == Lmu, Gimpl::twist.adjomega[mu] * U_mirror, U_mirror);
                Udag        = where(coor == 0, conjugate(phase) * Gimpl::twist.adjomega[mu] * Udag, Udag);
                Udag_mirror = where(coor == 0, Udag_mirror * Gimpl::twist.omega[mu], Udag_mirror);
            }
            else{
                U = where(coor == Lmu, phase * U, U);
                Udag = where(coor == 0, conjugate(phase) * Udag, Udag);
            }
            
            PokeIndex<LorentzIndex>(Uds, U,           mu     );
            PokeIndex<LorentzIndex>(Uds, Udag,        mu + 4 );
            PokeIndex<LorentzIndex>(Uds, U_mirror,    mu + 8 );
            PokeIndex<LorentzIndex>(Uds, Udag_mirror, mu + 12);
        }
    }
    
    inline void InsertForce4D(GaugeField &mat, FermionField &Btilde, FermionField &A,int mu){
        GaugeLinkField link(mat._grid);
        FermionField tmp = adj(A);
        link  = Btilde*tmp;
        link -= tmp*Btilde;
        PokeIndex<LorentzIndex>(mat,link,mu);
    }
    
    inline void InsertForce5D(GaugeField &mat, FermionField &Btilde, FermionField &Atilde,int mu){
        assert (0);
        // Must never hit
    }
};

/////////////////////////////////////////////////////////////////////////////
// Perturbative staggered fermion in the adjoint representation
/////////////////////////////////////////////////////////////////////////////
template <class S, class Representation = FundamentalRepresentation >
class PStaggeredAdjointImpl : public TwistedGaugeImpl<GaugeImplTypes_pt<S, Representation::Dimension > > {
    
public:
    
    // single order implementation
    typedef typename QCDpt::StaggeredAdjointImpl<S,Representation> SOimpl;
    
    typedef RealD  _Coeff_t ;
    static const int Dimension = Representation::Dimension;
    static const bool isFundamental = Representation::isFundamental;
    static const bool LsVectorised=false;
    typedef TwistedGaugeImpl<GaugeImplTypes_pt<S, Dimension > > Gimpl;
    
    //Necessary?
    constexpr bool is_fundamental() const{return Dimension == Nc ? 1 : 0;}
    
    typedef _Coeff_t Coeff_t;
    
    INHERIT_GIMPL_TYPES(Gimpl);
    
    template <typename vtype> using iImplSpinor            = iScalar<iScalar<iPert<iMatrix<vtype, Dimension>, Np> > >;
    template <typename vtype> using iImplHalfSpinor        = iScalar<iScalar<iPert<iVector<vtype, Dimension>, Np> > >;
    template <typename vtype> using iImplDoubledGaugeField = iVector<iScalar<iPert<iMatrix<vtype, Dimension>, Np> >, Nds>;
    template <typename vtype> using iImplPropagator        = iScalar<iScalar<iPert<iMatrix<vtype, Dimension>, Np> > >;
    
    typedef iImplSpinor<Simd>            SiteSpinor;
    typedef iImplHalfSpinor<Simd>        SiteHalfSpinor;
    typedef iImplDoubledGaugeField<Simd> SiteDoubledGaugeField;
    typedef iImplPropagator<Simd>        SitePropagator;
    
    typedef Lattice<SiteSpinor>            FermionField;
    typedef Lattice<SiteDoubledGaugeField> DoubledGaugeField;
    typedef Lattice<SitePropagator> PropagatorField;
    
    typedef SimpleCompressor<SiteSpinor> Compressor;
    typedef WilsonImplParams ImplParams; // because boundary phases are needed
    typedef CartesianStencil<SiteSpinor, SiteSpinor> StencilImpl;
    
    typedef Lattice<QCDpt::iSinglet<typename FermionField::vector_type> > LatticeSinglet;
    
    ImplParams Params;
    
    PStaggeredAdjointImpl(const ImplParams &p = ImplParams()) : Params(p){};
    
    inline void multLink(SiteSpinor &phi,
                         const SiteDoubledGaugeField &U,
                         const SiteSpinor &chi,
                         int mu,
                         StencilEntry *SE,
                         StencilImpl &St){
        SiteSpinor phitmp;
        mult(&phitmp(), &U(mu), &chi());
        mult(&phi(), &phitmp(), &U(mu+Nds));
    }
    inline void multLinkAdd(SiteSpinor &phi,
                            const SiteDoubledGaugeField &U,
                            const SiteSpinor &chi,
                            int mu,
                            StencilEntry *SE,
                            StencilImpl &St){
        SiteSpinor phitmp;
        mult(&phitmp(), &U(mu), &chi());
        mac(&phi(), &phitmp(), &U(mu+Nds));
    }
    
    template <class ref>
    inline void loadLinkElement(Simd &reg, ref &memory) {
        reg = memory;
    }
    
    inline void DoubleStore(GridBase *GaugeGrid, DoubledGaugeField &Uds, const GaugeField &Uthin) {
        conformable(Uds._grid, GaugeGrid);
        conformable(Uthin._grid, GaugeGrid);
        GaugeLinkField U(GaugeGrid);
        GaugeLinkField Udag(GaugeGrid);
        GaugeLinkField U_mirror(GaugeGrid);
        GaugeLinkField Udag_mirror(GaugeGrid);
        Lattice<iScalar<vInteger> > coor(GaugeGrid);
        typedef typename Simd::scalar_type scalar_type;
        
        for (int mu = 0; mu < Nd; mu++) {
            
            // Staggered Phase.
            Lattice<iScalar<vInteger> > coor(GaugeGrid);
            Lattice<iScalar<vInteger> > x(GaugeGrid); LatticeCoordinate(x,0);
            Lattice<iScalar<vInteger> > y(GaugeGrid); LatticeCoordinate(y,1);
            Lattice<iScalar<vInteger> > z(GaugeGrid); LatticeCoordinate(z,2);
            Lattice<iScalar<vInteger> > t(GaugeGrid); LatticeCoordinate(t,3);
            
            Lattice<iScalar<vInteger> > lin_z(GaugeGrid); lin_z=x+y;
            Lattice<iScalar<vInteger> > lin_t(GaugeGrid); lin_t=x+y+z;
            
            ComplexField phases(GaugeGrid);    phases=1.0;
            
            if ( mu == 1 ) phases = where( mod(x    ,2)==(Integer)0, phases,-phases);
            if ( mu == 2 ) phases = where( mod(lin_z,2)==(Integer)0, phases,-phases);
            if ( mu == 3 ) phases = where( mod(lin_t,2)==(Integer)0, phases,-phases);
            
            U      = PeekIndex<LorentzIndex>(Uthin, mu);
            Udag   = adj( Cshift(U, mu, -1));
            
            U_mirror    = adj(U);
            Udag_mirror = adj(Udag);
            U    = U    *phases;
            Udag = Udag *phases;
            
            
            // fermion phase and twist
            auto pha = Params.boundary_phases[mu];
            scalar_type phase( real(pha),imag(pha) );
            int Lmu = GaugeGrid->GlobalDimensions()[mu] - 1;
            LatticeCoordinate(coor, mu);
            if(istwisted(mu)){
                U        = where(coor == Lmu, phase * U * Gimpl::twist.omega[mu], U);
                U_mirror = where(coor == Lmu, Gimpl::twist.adjomega[mu] * U_mirror, U_mirror);
                Udag        = where(coor == 0, conjugate(phase) * Gimpl::twist.adjomega[mu] * Udag, Udag);
                Udag_mirror = where(coor == 0, Udag_mirror * Gimpl::twist.omega[mu], Udag_mirror);
            }
            else{
                U = where(coor == Lmu, phase * U, U);
                Udag = where(coor == 0, conjugate(phase) * Udag, Udag);
            }
            
            PokeIndex<LorentzIndex>(Uds, U,           mu     );
            PokeIndex<LorentzIndex>(Uds, Udag,        mu + 4 );
            PokeIndex<LorentzIndex>(Uds, U_mirror,    mu + 8 );
            PokeIndex<LorentzIndex>(Uds, Udag_mirror, mu + 12);
        }
    }
    
    inline void InsertForce4D(GaugeField &mat, FermionField &Btilde, FermionField &A,int mu){
        assert(0);
    }
    
    inline void InsertForce5D(GaugeField &mat, FermionField &Btilde, FermionField &Atilde,int mu){
        assert (0);
        // Must never hit
    }
};
    
    /////////////////////////////////////////////////////////////////////////////
    // Nonperturbative naive staggered fermion in the adjoint representation
    /////////////////////////////////////////////////////////////////////////////
template <class S, class Representation = FundamentalRepresentation >
//class StaggeredAdjointImplNP : public PeriodicGaugeImpl<GaugeImplTypes<S, Representation::Dimension > > {
class StaggeredAdjointImplNP : public TwistedGaugeImplNP<GaugeImplTypes<S, Representation::Dimension > > {//$//
    
public:
    
    typedef RealD  _Coeff_t ;
    static const int Dimension = Representation::Dimension;
    static const bool isFundamental = Representation::isFundamental;
    static const bool LsVectorised=false;
    typedef TwistedGaugeImplNP<GaugeImplTypes<S, Dimension > > Gimpl;//$//
//    typedef PeriodicGaugeImpl<GaugeImplTypes<S, Dimension > > Gimpl;
    
    //Necessary?
    constexpr bool is_fundamental() const{return Dimension == Nc ? 1 : 0;}
    
    typedef _Coeff_t Coeff_t;
    
    INHERIT_GIMPL_TYPES(Gimpl);
    
    template <typename vtype> using iImplSpinor            = iScalar<iScalar<iMatrix<vtype, Dimension> > >;
    template <typename vtype> using iImplHalfSpinor        = iScalar<iScalar<iMatrix<vtype, Dimension> > >;
    template <typename vtype> using iImplDoubledGaugeField = iVector<iScalar<iMatrix<vtype, Dimension> >, 2*Nds>;
    template <typename vtype> using iImplPropagator        = iScalar<iScalar<iMatrix<vtype, Dimension> > >;
    
    typedef iImplSpinor<Simd>            SiteSpinor;
    typedef iImplHalfSpinor<Simd>        SiteHalfSpinor;
    typedef iImplDoubledGaugeField<Simd> SiteDoubledGaugeField;
    typedef iImplPropagator<Simd>        SitePropagator;
    
    typedef Lattice<SiteSpinor>            FermionField;
    typedef Lattice<SiteDoubledGaugeField> DoubledGaugeField;
    typedef Lattice<SitePropagator> PropagatorField;
    
    typedef SimpleCompressor<SiteSpinor> Compressor;
    typedef StaggeredImplParams ImplParams;
    typedef CartesianStencil<SiteSpinor, SiteSpinor> StencilImpl;
    
    typedef Lattice<iSinglet<typename FermionField::vector_type> > LatticeSinglet;
    
    ImplParams Params;
    
    StaggeredAdjointImplNP(const ImplParams &p = ImplParams()) : Params(p){};
    
    inline void multLink(SiteSpinor &phi,
                         const SiteDoubledGaugeField &U,
                         const SiteSpinor &chi,
                         int mu,
                         StencilEntry *SE,
                         StencilImpl &St) {
        SiteSpinor phitmp;
        mult(&phitmp(), &U(mu), &chi());
        mult(&phi(), &phitmp(), &U(mu+Nds));
    }
    
    inline void multLinkAdd(SiteSpinor &phi,
                            const SiteDoubledGaugeField &U,
                            const SiteSpinor &chi,
                            int mu,
                            StencilEntry *SE,
                            StencilImpl &St) {
        SiteSpinor phitmp;
        mult(&phitmp(), &U(mu), &chi());
        mac(&phi(), &phitmp(), &U(mu+Nds));
    }
    
    template <class ref>
    inline void loadLinkElement(Simd &reg, ref &memory) {
        reg = memory;
    }
    
    inline void DoubleStore(GridBase *GaugeGrid, DoubledGaugeField &Uds, const GaugeField &Uthin){
        conformable(Uds._grid, GaugeGrid);
        conformable(Uthin._grid, GaugeGrid);
        GaugeLinkField U(GaugeGrid);
        GaugeLinkField Udag(GaugeGrid);
        GaugeLinkField U_mirror(GaugeGrid);
        GaugeLinkField Udag_mirror(GaugeGrid);
        Lattice<iScalar<vInteger> > coor(GaugeGrid);
        typedef typename Simd::scalar_type scalar_type;
        
        for (int mu = 0; mu < Nd; mu++) {
            
            // Staggered Phase.
            Lattice<iScalar<vInteger> > coor(GaugeGrid);
            Lattice<iScalar<vInteger> > x(GaugeGrid); LatticeCoordinate(x,0);
            Lattice<iScalar<vInteger> > y(GaugeGrid); LatticeCoordinate(y,1);
            Lattice<iScalar<vInteger> > z(GaugeGrid); LatticeCoordinate(z,2);
            Lattice<iScalar<vInteger> > t(GaugeGrid); LatticeCoordinate(t,3);
            
            Lattice<iScalar<vInteger> > lin_z(GaugeGrid); lin_z=x+y;
            Lattice<iScalar<vInteger> > lin_t(GaugeGrid); lin_t=x+y+z;
            
            ComplexField phases(GaugeGrid);    phases=1.0;
            
            if ( mu == 1 ) phases = where( mod(x    ,2)==(Integer)0, phases,-phases);
            if ( mu == 2 ) phases = where( mod(lin_z,2)==(Integer)0, phases,-phases);
            if ( mu == 3 ) phases = where( mod(lin_t,2)==(Integer)0, phases,-phases);
            
            U      = PeekIndex<LorentzIndex>(Uthin, mu);
            Udag   = adj( Cshift(U, mu, -1));
            
            U_mirror    = adj(U);
            Udag_mirror = adj(Udag);
            U    = U    *phases;
            Udag = Udag *phases;
            
            
            // twist//$//
            int Lmu = GaugeGrid->GlobalDimensions()[mu] - 1;
            LatticeCoordinate(coor, mu);
            if(istwisted(mu)){
                U        = where(coor == Lmu, U * Gimpl::twist.omega[mu], U);
                U_mirror = where(coor == Lmu, Gimpl::twist.adjomega[mu] * U_mirror, U_mirror);
                Udag        = where(coor == 0, Gimpl::twist.adjomega[mu] * Udag, Udag);
                Udag_mirror = where(coor == 0, Udag_mirror * Gimpl::twist.omega[mu], Udag_mirror);
            }
            
            PokeIndex<LorentzIndex>(Uds, U,           mu     );
            PokeIndex<LorentzIndex>(Uds, Udag,        mu + 4 );
            PokeIndex<LorentzIndex>(Uds, U_mirror,    mu + 8 );
            PokeIndex<LorentzIndex>(Uds, Udag_mirror, mu + 12);
        }
    }
    
    inline void InsertForce4D(GaugeField &mat, FermionField &Btilde, FermionField &A,int mu){
        GaugeLinkField link(mat._grid);
        FermionField tmp = adj(A);
        link  = Btilde*tmp;
        link -= tmp*Btilde;
        PokeIndex<LorentzIndex>(mat,link,mu);
    }
    
    inline void InsertForce5D(GaugeField &mat, FermionField &Btilde, FermionField &Atilde,int mu){
        assert (0);
        // Must never hit
    }
};
    
typedef PWilsonSmellImpl<vComplex,  FundamentalRepresentation, CoeffReal > PWilsonSmellImplR;  // Real.. whichever prec
typedef PWilsonSmellImpl<vComplexF, FundamentalRepresentation, CoeffReal > PWilsonSmellImplF;  // Float
typedef PWilsonSmellImpl<vComplexD, FundamentalRepresentation, CoeffReal > PWilsonSmellImplD;  // Double

typedef WilsonSmellImpl<vComplex,  FundamentalRepresentation, CoeffReal > WilsonSmellImplR;  // Real.. whichever prec
typedef WilsonSmellImpl<vComplexF, FundamentalRepresentation, CoeffReal > WilsonSmellImplF;  // Float
typedef WilsonSmellImpl<vComplexD, FundamentalRepresentation, CoeffReal > WilsonSmellImplD;  // Double


typedef PStaggeredSmellImpl<vComplex,  FundamentalRepresentation > PStaggeredSmellImplR;  // Real.. whichever prec
typedef PStaggeredSmellImpl<vComplexF, FundamentalRepresentation > PStaggeredSmellImplF;  // Float
typedef PStaggeredSmellImpl<vComplexD, FundamentalRepresentation > PStaggeredSmellImplD;  // Double

typedef StaggeredSmellImpl<vComplex,  FundamentalRepresentation > StaggeredSmellImplR;  // Real.. whichever prec
typedef StaggeredSmellImpl<vComplexF, FundamentalRepresentation > StaggeredSmellImplF;  // Float
typedef StaggeredSmellImpl<vComplexD, FundamentalRepresentation > StaggeredSmellImplD;  // Double

typedef NaiveStaggeredImpl<vComplex,  FundamentalRepresentation > NaiveStaggeredImplR;  // Real.. whichever prec
typedef NaiveStaggeredImpl<vComplexF, FundamentalRepresentation > NaiveStaggeredImplF;  // Float
typedef NaiveStaggeredImpl<vComplexD, FundamentalRepresentation > NaiveStaggeredImplD;  // Double
    
typedef StaggeredAdjointImpl<vComplex,  FundamentalRepresentation > StaggeredAdjointImplR;  // Real.. whichever prec
typedef StaggeredAdjointImpl<vComplexF, FundamentalRepresentation > StaggeredAdjointImplF;  // Float
typedef StaggeredAdjointImpl<vComplexD, FundamentalRepresentation > StaggeredAdjointImplD;  // Double
typedef PStaggeredAdjointImpl<vComplex,  FundamentalRepresentation > PStaggeredAdjointImplR;  // Real.. whichever prec
typedef PStaggeredAdjointImpl<vComplexF, FundamentalRepresentation > PStaggeredAdjointImplF;  // Float
typedef PStaggeredAdjointImpl<vComplexD, FundamentalRepresentation > PStaggeredAdjointImplD;  // Double

typedef StaggeredAdjointImplNP<vComplex,  FundamentalRepresentation > StaggeredAdjointImplNPR;  // Real.. whichever prec
typedef StaggeredAdjointImplNP<vComplexF, FundamentalRepresentation > StaggeredAdjointImplNPF;  // Float
typedef StaggeredAdjointImplNP<vComplexD, FundamentalRepresentation > StaggeredAdjointImplNPD;  // Double
}}}

#endif
