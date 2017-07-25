    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/QCD.h

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>
Author: neo <cossu@post.kek.jp>
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

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#ifndef GRID_QCD_PT_BASE_H
#define GRID_QCD_PT_BASE_H
namespace Grid{
namespace QCD {
namespace QCDpt {
    
    static const int Np=3;
    
    #define PertIndex      2
    #define ColourIndexPT  3
    
    
    template<typename vtype> using iSinglet                       = iScalar<iScalar<iScalar<iScalar<vtype> > > >;
    template<typename vtype> using iPertSinglet                   = iScalar<iScalar<iPert<iScalar<vtype>, Np> > >;
    
    template<typename vtype> using iColourMatrix                  = iScalar<iScalar<iScalar<iMatrix<vtype, Nc> > > > ;
    template<typename vtype> using iLorentzColourMatrix           = iVector<iScalar<iScalar<iMatrix<vtype, Nc> > >, Nd > ;
    template<typename vtype> using iPertColourMatrix              = iScalar<iScalar<iPert<iMatrix<vtype, Nc>, Np> > > ;
    template<typename vtype> using iLorentzPertColourMatrix       = iVector<iScalar<iPert<iMatrix<vtype, Nc>, Np> >, Nd > ;
    
    template<typename vtype> using iColourVector                  = iScalar<iScalar<iScalar<iVector<vtype, Nc> > > >;
    template<typename vtype> using iSpinColourVector              = iScalar<iVector<iScalar<iVector<vtype, Nc> >, Ns> >;
    template<typename vtype> using iPertColourVector              = iScalar<iScalar<iPert<iVector<vtype, Nc>, Np> > >;
    template<typename vtype> using iSpinPertColourVector          = iScalar<iVector<iPert<iVector<vtype, Nc>, Np>, Ns> >;
    
    template<typename vtype> using iSpinColourMatrix              = iScalar<iVector<iScalar<iMatrix<vtype, Nc> >, Ns> >;
    template<typename vtype> using iSpinPertColourMatrix          = iScalar<iVector<iPert<iMatrix<vtype, Nc>, Np>, Ns> >;
    
    
    // perturbative fundamental types (minimal nesting)
    
    typedef iPert<Complex, Np >                          PComplex;
    typedef iPert<ComplexF, Np>                          PComplexF;
    typedef iPert<ComplexD, Np>                          PComplexD;
    
    typedef iPert<vComplex, Np >                         vPComplex ;
    typedef iPert<vComplexF, Np>                         vPComplexF;
    typedef iPert<vComplexD, Np>                         vPComplexD;
    
    typedef iPert<Real, Np >                             PReal;
    typedef iPert<RealF, Np>                             PRealF;
    typedef iPert<RealD, Np>                             PRealD;
    
    typedef iPert<vReal, Np >                            vPReal;
    typedef iPert<vRealF, Np>                            vPRealF;
    typedef iPert<vRealD, Np>                            vPRealD;
    
    typedef iPert<Integer, Np >                          PInteger;
    typedef iPert<vInteger, Np>                          vPInteger;
    
    // singlets
    
    typedef iSinglet<Complex >                           TComplex;
    typedef iSinglet<ComplexF>                           TComplexF;
    typedef iSinglet<ComplexD>                           TComplexD;
    
    typedef iSinglet<vComplex >                          vTComplex ;
    typedef iSinglet<vComplexF>                          vTComplexF;
    typedef iSinglet<vComplexD>                          vTComplexD;
    
    typedef iSinglet<Real >                              TReal;
    typedef iSinglet<RealF>                              TRealF;
    typedef iSinglet<RealD>                              TRealD;
    
    typedef iSinglet<vReal >                             vTReal;
    typedef iSinglet<vRealF>                             vTRealF;
    typedef iSinglet<vRealD>                             vTRealD;
    
    typedef iSinglet<Integer >                           TInteger;
    typedef iSinglet<vInteger>                           vTInteger;
    
    typedef iPertSinglet<Complex >                       TPertComplex;
    typedef iPertSinglet<ComplexF>                       TPertComplexF;
    typedef iPertSinglet<ComplexD>                       TPertComplexD;
    
    typedef iPertSinglet<vComplex >                      vTPertComplex ;
    typedef iPertSinglet<vComplexF>                      vTPertComplexF;
    typedef iPertSinglet<vComplexD>                      vTPertComplexD;
    
    typedef iPertSinglet<Real >                          TPertReal;
    typedef iPertSinglet<RealF>                          TPertRealF;
    typedef iPertSinglet<RealD>                          TPertRealD;
    
    typedef iPertSinglet<vReal >                         vTPertReal;
    typedef iPertSinglet<vRealF>                         vTPertRealF;
    typedef iPertSinglet<vRealD>                         vTPertRealD;
    
    typedef iPertSinglet<Integer >                       TPertInteger;
    typedef iPertSinglet<vInteger>                       vTPertInteger;
    
    // gauge fields
    
    typedef iColourMatrix<Complex  >                     ColourMatrix;
    typedef iColourMatrix<ComplexF >                     ColourMatrixF;
    typedef iColourMatrix<ComplexD >                     ColourMatrixD;
    
    typedef iColourMatrix<vComplex  >                    vColourMatrix;
    typedef iColourMatrix<vComplexF >                    vColourMatrixF;
    typedef iColourMatrix<vComplexD >                    vColourMatrixD;
    
    typedef iLorentzColourMatrix<Complex  >              LorentzColourMatrix;
    typedef iLorentzColourMatrix<ComplexF >              LorentzColourMatrixF;
    typedef iLorentzColourMatrix<ComplexD >              LorentzColourMatrixD;
    
    typedef iLorentzColourMatrix<vComplex  >             vLorentzColourMatrix;
    typedef iLorentzColourMatrix<vComplexF >             vLorentzColourMatrixF;
    typedef iLorentzColourMatrix<vComplexD >             vLorentzColourMatrixD;
    
    typedef iPertColourMatrix<Complex  >                 PertColourMatrix;
    typedef iPertColourMatrix<ComplexF >                 PertColourMatrixF;
    typedef iPertColourMatrix<ComplexD >                 PertColourMatrixD;
    
    typedef iPertColourMatrix<vComplex  >                vPertColourMatrix;
    typedef iPertColourMatrix<vComplexF >                vPertColourMatrixF;
    typedef iPertColourMatrix<vComplexD >                vPertColourMatrixD;
    
    typedef iLorentzPertColourMatrix<Complex  >          LorentzPertColourMatrix;
    typedef iLorentzPertColourMatrix<ComplexF >          LorentzPertColourMatrixF;
    typedef iLorentzPertColourMatrix<ComplexD >          LorentzPertColourMatrixD;
    
    typedef iLorentzPertColourMatrix<vComplex  >         vLorentzPertColourMatrix;
    typedef iLorentzPertColourMatrix<vComplexF >         vLorentzPertColourMatrixF;
    typedef iLorentzPertColourMatrix<vComplexD >         vLorentzPertColourMatrixD;
    
    // fermions without smell
    
    typedef iColourVector<Complex  >                     ColourVector;
    typedef iColourVector<ComplexF >                     ColourVectorF;
    typedef iColourVector<ComplexD >                     ColourVectorD;
    
    typedef iColourVector<vComplex  >                    vColourVector;
    typedef iColourVector<vComplexF >                    vColourVectorF;
    typedef iColourVector<vComplexD >                    vColourVectorD;
    
    typedef iSpinColourVector<Complex  >                 SpinColourVector;
    typedef iSpinColourVector<ComplexF >                 SpinColourVectorF;
    typedef iSpinColourVector<ComplexD >                 SpinColourVectorD;
    
    typedef iSpinColourVector<vComplex  >                vSpinColourVector;
    typedef iSpinColourVector<vComplexF >                vSpinColourVectorF;
    typedef iSpinColourVector<vComplexD >                vSpinColourVectorD;
    
    typedef iPertColourVector<Complex  >                 PertColourVector;
    typedef iPertColourVector<ComplexF >                 PertColourVectorF;
    typedef iPertColourVector<ComplexD >                 PertColourVectorD;
    
    typedef iPertColourVector<vComplex  >                vPertColourVector;
    typedef iPertColourVector<vComplexF >                vPertColourVectorF;
    typedef iPertColourVector<vComplexD >                vPertColourVectorD;
    
    typedef iSpinPertColourVector<Complex  >             SpinPertColourVector;
    typedef iSpinPertColourVector<ComplexF >             SpinPertColourVectorF;
    typedef iSpinPertColourVector<ComplexD >             SpinPertColourVectorD;
    
    typedef iSpinPertColourVector<vComplex  >            vSpinPertColourVector;
    typedef iSpinPertColourVector<vComplexF >            vSpinPertColourVectorF;
    typedef iSpinPertColourVector<vComplexD >            vSpinPertColourVectorD;
    
    // fermions with smell
    
    typedef iSpinColourMatrix<Complex  >                 SpinColourMatrix;
    typedef iSpinColourMatrix<ComplexF >                 SpinColourMatrixF;
    typedef iSpinColourMatrix<ComplexD >                 SpinColourMatrixD;
    
    typedef iSpinColourMatrix<vComplex  >                vSpinColourMatrix;
    typedef iSpinColourMatrix<vComplexF >                vSpinColourMatrixF;
    typedef iSpinColourMatrix<vComplexD >                vSpinColourMatrixD;
    
    typedef iSpinPertColourMatrix<Complex  >             SpinPertColourMatrix;
    typedef iSpinPertColourMatrix<ComplexF >             SpinPertColourMatrixF;
    typedef iSpinPertColourMatrix<ComplexD >             SpinPertColourMatrixD;
    
    typedef iSpinPertColourMatrix<vComplex  >            vSpinPertColourMatrix;
    typedef iSpinPertColourMatrix<vComplexF >            vSpinPertColourMatrixF;
    typedef iSpinPertColourMatrix<vComplexD >            vSpinPertColourMatrixD;
    
    
    // lattice singlets
    
    typedef Lattice<vTComplex>                           LatticeComplex;
    typedef Lattice<vTComplexF>                          LatticeComplexF;
    typedef Lattice<vTComplexD>                          LatticeComplexD;
    
    typedef Lattice<vTReal>                              LatticeReal;
    typedef Lattice<vTRealF>                             LatticeRealF;
    typedef Lattice<vTRealD>                             LatticeRealD;
    
    typedef Lattice<vTInteger>                           LatticeInteger;

    typedef Lattice<vTPertComplex>                       LatticePertComplex;
    typedef Lattice<vTPertComplexF>                      LatticePertComplexF;
    typedef Lattice<vTPertComplexD>                      LatticePertComplexD;
    
    typedef Lattice<vTPertReal>                          LatticePertReal;
    typedef Lattice<vTPertRealF>                         LatticePertRealF;
    typedef Lattice<vTPertRealD>                         LatticePertRealD;
    
    typedef Lattice<vTPertInteger>                       LatticePertInteger;
    
    // lattice gauge fields
    
    typedef Lattice<vColourMatrix>                       LatticeColourMatrix;
    typedef Lattice<vColourMatrixF>                      LatticeColourMatrixF;
    typedef Lattice<vColourMatrixD>                      LatticeColourMatrixD;
    
    typedef Lattice<vLorentzColourMatrix>                LatticeLorentzColourMatrix;
    typedef Lattice<vLorentzColourMatrixF>               LatticeLorentzColourMatrixF;
    typedef Lattice<vLorentzColourMatrixD>               LatticeLorentzColourMatrixD;
    
    typedef Lattice<vPertColourMatrix>                   LatticePertColourMatrix;
    typedef Lattice<vPertColourMatrixF>                  LatticePertColourMatrixF;
    typedef Lattice<vPertColourMatrixD>                  LatticePertColourMatrixD;
    
    typedef Lattice<vLorentzPertColourMatrix>            LatticeLorentzPertColourMatrix;
    typedef Lattice<vLorentzPertColourMatrixF>           LatticeLorentzPertColourMatrixF;
    typedef Lattice<vLorentzPertColourMatrixD>           LatticeLorentzPertColourMatrixD;
    
    // lattice fermions without smell
    
    typedef Lattice<vColourVector>                       LatticeColourVector;
    typedef Lattice<vColourVectorF>                      LatticeColourVectorF;
    typedef Lattice<vColourVectorD>                      LatticeColourVectorD;
    
    typedef Lattice<vSpinColourVector>                   LatticeSpinColourVector;
    typedef Lattice<vSpinColourVectorF>                  LatticeSpinColourVectorF;
    typedef Lattice<vSpinColourVectorD>                  LatticeSpinColourVectorD;
    
    typedef Lattice<vPertColourVector>                   LatticePertColourVector;
    typedef Lattice<vPertColourVectorF>                  LatticePertColourVectorF;
    typedef Lattice<vPertColourVectorD>                  LatticePertColourVectorD;
    
    typedef Lattice<vSpinPertColourVector>               LatticeSpinPertColourVector;
    typedef Lattice<vSpinPertColourVectorF>              LatticeSpinPertColourVectorF;
    typedef Lattice<vSpinPertColourVectorD>              LatticeSpinPertColourVectorD;
    
    // lattice fermions with smell
    
    typedef Lattice<vSpinColourMatrix>                   LatticeSpinColourMatrix;
    typedef Lattice<vSpinColourMatrixF>                  LatticeSpinColourMatrixF;
    typedef Lattice<vSpinColourMatrixD>                  LatticeSpinColourMatrixD;
    
    typedef Lattice<vSpinPertColourMatrix>               LatticeSpinPertColourMatrix;
    typedef Lattice<vSpinPertColourMatrixF>              LatticeSpinPertColourMatrixF;
    typedef Lattice<vSpinPertColourMatrixD>              LatticeSpinPertColourMatrixD;
    
    
    
    // Physical names
    
    typedef LatticeLorentzPertColourMatrix               LatticeGaugeField;
    typedef LatticeLorentzPertColourMatrixF              LatticeGaugeFieldF;
    typedef LatticeLorentzPertColourMatrixD              LatticeGaugeFieldD;
    
    typedef LatticeSpinPertColourMatrix                  LatticeFermion;
    typedef LatticeSpinPertColourMatrixF                 LatticeFermionF;
    typedef LatticeSpinPertColourMatrixD                 LatticeFermionD;
    
    
    //////////////////////////////////////////////////////////////////////////////
    // Redefine templates leaning on ColourIndex attribute
    //////////////////////////////////////////////////////////////////////////////
    
    template<class vobj> auto peekColour(const vobj &rhs,int i) -> decltype(PeekIndex<ColourIndexPT>(rhs,0))
    {
      return PeekIndex<ColourIndexPT>(rhs,i);
    }
    template<class vobj> auto peekColour(const vobj &rhs,int i,int j) -> decltype(PeekIndex<ColourIndexPT>(rhs,0,0))
    {
      return PeekIndex<ColourIndexPT>(rhs,i,j);
    }
    template<class vobj> auto peekColour(const Lattice<vobj> &rhs,int i) -> decltype(PeekIndex<ColourIndexPT>(rhs,0))
    {
      return PeekIndex<ColourIndexPT>(rhs,i);
    }
    template<class vobj> auto peekColour(const Lattice<vobj> &rhs,int i,int j) -> decltype(PeekIndex<ColourIndexPT>(rhs,0,0))
    {
      return PeekIndex<ColourIndexPT>(rhs,i,j);
    }
    
    template<class vobj> auto peekPert(const vobj &rhs,int i) -> decltype(PeekIndex<PertIndex>(rhs,0))
    {
      return PeekIndex<PertIndex>(rhs,i);
    }
    
    template<class vobj> auto peekPert(const Lattice<vobj> &rhs,int i) -> decltype(PeekIndex<PertIndex>(rhs,0))
    {
      return PeekIndex<PertIndex>(rhs,i);
    }
    
    
    
    
    template<class vobj>
      void pokeColour(Lattice<vobj> &lhs,
              const Lattice<decltype(peekIndex<ColourIndexPT>(lhs._odata[0],0))> & rhs,
              int i)
    {
      PokeIndex<ColourIndexPT>(lhs,rhs,i);
    }
    template<class vobj> 
      void pokeColour(Lattice<vobj> &lhs,
              const Lattice<decltype(peekIndex<ColourIndexPT>(lhs._odata[0],0,0))> & rhs,
              int i,int j)
    {
      PokeIndex<ColourIndexPT>(lhs,rhs,i,j);
    }
    
    template<class vobj> void pokeColour(vobj &lhs,const decltype(peekIndex<ColourIndexPT>(lhs,0)) & rhs,int i)
    {
      pokeIndex<ColourIndexPT>(lhs,rhs,i);
    }
    template<class vobj> void pokeColour(vobj &lhs,const decltype(peekIndex<ColourIndexPT>(lhs,0,0)) & rhs,int i,int j)
    {
      pokeIndex<ColourIndexPT>(lhs,rhs,i,j);
    }
    
    template<class vobj>
      void pokePert(Lattice<vobj> &lhs,
              const Lattice<decltype(peekIndex<PertIndex>(lhs._odata[0],0))> & rhs,
              int i)
    {
      PokeIndex<PertIndex>(lhs,rhs,i);
    }
    template<class vobj> void pokePert(vobj &lhs,const decltype(peekIndex<PertIndex>(lhs,0)) & rhs,int i)
    {
      pokeIndex<PertIndex>(lhs,rhs,i);
    }
    
    
    
    
    template<int Index,class vobj> inline Lattice<vobj> transposeColour(const Lattice<vobj> &lhs){
      return transposeIndex<ColourIndexPT>(lhs);
    }
    template<int Index,class vobj> inline vobj transposeColour(const vobj &lhs){
      return transposeIndex<ColourIndexPT>(lhs);
    }
    
    template<int Index,class vobj>
    inline auto traceColour(const Lattice<vobj> &lhs) -> Lattice<decltype(traceIndex<ColourIndexPT>(lhs._odata[0]))>
    {
      return traceIndex<ColourIndexPT>(lhs);
    }
    template<int Index,class vobj>
    inline auto traceColour(const vobj &lhs) -> Lattice<decltype(traceIndex<ColourIndexPT>(lhs))>
    {
      return traceIndex<ColourIndexPT>(lhs);
    }

}   //namespace QCDpt
}   //namespace QCD
} // Grid

#endif
