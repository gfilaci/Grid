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
         // and gives a matrix in colour space
         //$// A is the noise, so it should not be perturbative...
         //$// need to resolve nesting because mult(iVector,iVector) is not defined...
        FermionField tmp = adj(A);
        parallel_for(int ss=0;ss<mat._grid->oSites();ss++){
            for (int alpha=0; alpha<Ns; alpha++)
                tmp._odata[ss]._internal._internal[alpha] = Btilde._odata[ss]._internal._internal[alpha] * tmp._odata[ss]._internal._internal[alpha];
        }
        link = TraceIndex<SpinIndex>(tmp);
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
         // and gives a matrix in colour space
         //$// A is the noise, so it should not be perturbative...
         //$// need to resolve nesting because mult(iVector,iVector) is not defined...
        FermionField tmp = adj(A);
        parallel_for(int ss=0;ss<mat._grid->oSites();ss++){
            for (int alpha=0; alpha<Ns; alpha++)
                tmp._odata[ss]._internal._internal[alpha] = Btilde._odata[ss]._internal._internal[alpha] * tmp._odata[ss]._internal._internal[alpha];
        }
        link = TraceIndex<SpinIndex>(tmp);
        PokeIndex<LorentzIndex>(mat,link,mu);
    }
      
    inline void InsertForce5D(GaugeField &mat, FermionField &Btilde, FermionField &Atilde,int mu){
        assert(0);
    }
  };


  
typedef PWilsonSmellImpl<vComplex,  FundamentalRepresentation, CoeffReal > PWilsonSmellImplR;  // Real.. whichever prec
typedef PWilsonSmellImpl<vComplexF, FundamentalRepresentation, CoeffReal > PWilsonSmellImplF;  // Float
typedef PWilsonSmellImpl<vComplexD, FundamentalRepresentation, CoeffReal > PWilsonSmellImplD;  // Double

typedef WilsonSmellImpl<vComplex,  FundamentalRepresentation, CoeffReal > WilsonSmellImplR;  // Real.. whichever prec
typedef WilsonSmellImpl<vComplexF, FundamentalRepresentation, CoeffReal > WilsonSmellImplF;  // Float
typedef WilsonSmellImpl<vComplexD, FundamentalRepresentation, CoeffReal > WilsonSmellImplD;  // Double

}}}

#endif
