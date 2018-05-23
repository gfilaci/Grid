/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/fermion/NaiveStaggeredKernels.cc

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
#include <Grid/qcd/action/fermion/FermionCore.h>

namespace Grid {
namespace QCD {

int NaiveStaggeredKernelsStatic::Opt= NaiveStaggeredKernelsStatic::OptGeneric;

template <class Impl>
NaiveStaggeredKernels<Impl>::NaiveStaggeredKernels(const ImplParams &p) : Base(p){};

////////////////////////////////////////////
// Generic implementation; move to different file?
////////////////////////////////////////////

template <class Impl>
void NaiveStaggeredKernels<Impl>::DhopSiteDepth(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U,
					   SiteSpinor *buf, int sF,
					   int sU, const FermionField &in, SiteSpinor &out,int threeLink) {
  const SiteSpinor *chi_p;
  SiteSpinor chi;
  SiteSpinor Uchi;
  StencilEntry *SE;
  int ptype;
  int skew = 0;
  
  ///////////////////////////
  // Xp
  ///////////////////////////

  SE = st.GetEntry(ptype, Xp+skew, sF);
  if (SE->_is_local) {
    if (SE->_permute) {
      chi_p = &chi;
      permute(chi,  in._odata[SE->_offset], ptype);
    } else {
      chi_p = &in._odata[SE->_offset];
    }
  } else {
    chi_p = &buf[SE->_offset];
  }
  Impl::multLink(Uchi, U._odata[sU], *chi_p, Xp, SE, st);

  ///////////////////////////
  // Yp
  ///////////////////////////
  SE = st.GetEntry(ptype, Yp+skew, sF);
  if (SE->_is_local) {
    if (SE->_permute) {
      chi_p = &chi;
      permute(chi,  in._odata[SE->_offset], ptype);
    } else {
      chi_p = &in._odata[SE->_offset];
    }
  } else {
    chi_p = &buf[SE->_offset];
  }
  Impl::multLinkAdd(Uchi, U._odata[sU], *chi_p, Yp, SE, st);

  ///////////////////////////
  // Zp
  ///////////////////////////
  SE = st.GetEntry(ptype, Zp+skew, sF);
  if (SE->_is_local) {
    if (SE->_permute) {
      chi_p = &chi;
      permute(chi,  in._odata[SE->_offset], ptype);
    } else {
      chi_p = &in._odata[SE->_offset];
    }
  } else {
    chi_p = &buf[SE->_offset];
  }
  Impl::multLinkAdd(Uchi, U._odata[sU], *chi_p, Zp, SE, st);

  ///////////////////////////
  // Tp
  ///////////////////////////
  SE = st.GetEntry(ptype, Tp+skew, sF);
  if (SE->_is_local) {
    if (SE->_permute) {
      chi_p = &chi;
      permute(chi,  in._odata[SE->_offset], ptype);
    } else {
      chi_p = &in._odata[SE->_offset];
    }
  } else {
    chi_p = &buf[SE->_offset];
  }
  Impl::multLinkAdd(Uchi, U._odata[sU], *chi_p, Tp, SE, st);

  ///////////////////////////
  // Xm
  ///////////////////////////
  SE = st.GetEntry(ptype, Xm+skew, sF);
  if (SE->_is_local) {
    if (SE->_permute) {
      chi_p = &chi;
      permute(chi,  in._odata[SE->_offset], ptype);
    } else {
      chi_p = &in._odata[SE->_offset];
    }
  } else {
    chi_p = &buf[SE->_offset];
  }
  Impl::multLinkAdd(Uchi, U._odata[sU], *chi_p, Xm, SE, st);

  ///////////////////////////
  // Ym
  ///////////////////////////
  SE = st.GetEntry(ptype, Ym+skew, sF);
  if (SE->_is_local) {
    if (SE->_permute) {
      chi_p = &chi;
      permute(chi,  in._odata[SE->_offset], ptype);
    } else {
      chi_p = &in._odata[SE->_offset];
    }
  } else {
    chi_p = &buf[SE->_offset];
  }
  Impl::multLinkAdd(Uchi, U._odata[sU], *chi_p, Ym, SE, st);

  ///////////////////////////
  // Zm
  ///////////////////////////
  SE = st.GetEntry(ptype, Zm+skew, sF);
  if (SE->_is_local) {
    if (SE->_permute) {
      chi_p = &chi;
      permute(chi,  in._odata[SE->_offset], ptype);
    } else {
      chi_p = &in._odata[SE->_offset];
    }
  } else {
    chi_p = &buf[SE->_offset];
  }
  Impl::multLinkAdd(Uchi, U._odata[sU], *chi_p, Zm, SE, st);

  ///////////////////////////
  // Tm
  ///////////////////////////
  SE = st.GetEntry(ptype, Tm+skew, sF);
  if (SE->_is_local) {
    if (SE->_permute) {
      chi_p = &chi;
      permute(chi,  in._odata[SE->_offset], ptype);
    } else {
      chi_p = &in._odata[SE->_offset];
    }
  } else {
    chi_p = &buf[SE->_offset];
  }
  Impl::multLinkAdd(Uchi, U._odata[sU], *chi_p, Tm, SE, st);

  vstream(out, Uchi);
};

template <class Impl>
void NaiveStaggeredKernels<Impl>::DhopSiteDag(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, SiteSpinor *buf, int LLs, int sU, const FermionField &in, FermionField &out) {
  SiteSpinor naive;
  int oneLink =0;
  switch(Opt) {
  case OptGeneric:
    for(int s=0;s<LLs;s++){
       int sF=s+LLs*sU;
       DhopSiteDepth(st,lo,U,buf,sF,sU,in,naive,oneLink);
       out._odata[sF] =-naive; // antihermitian
     }
    break;
  default:
    std::cout<<"Oops Opt = "<<Opt<<std::endl;
    assert(0);
    break;
  }
};

template <class Impl>
void NaiveStaggeredKernels<Impl>::DhopSite(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, SiteSpinor *buf, int LLs, int sU, const FermionField &in, FermionField &out)
{
  int oneLink  =0;
  SiteSpinor naive;
  switch(Opt) {
  case OptGeneric:
    for(int s=0;s<LLs;s++){
      int sF=LLs*sU+s;
      DhopSiteDepth(st,lo,U,buf,sF,sU,in,naive,oneLink);
      out._odata[sF] =naive;
    }
    break;
  default:
    std::cout<<"Oops Opt = "<<Opt<<std::endl;
    assert(0);
    break;
  }
};

template <class Impl>
void NaiveStaggeredKernels<Impl>::DhopDir( StencilImpl &st, DoubledGaugeField &U, SiteSpinor *buf, int sF,
				      int sU, const FermionField &in, FermionField &out, int dir) 
{
  
  const SiteSpinor *chi_p;
  SiteSpinor chi;
  SiteSpinor Uchi;
  StencilEntry *SE;
  int ptype;
  
  SE = st.GetEntry(ptype, dir, sF);
  if (SE->_is_local) {
    if (SE->_permute) {
      chi_p = &chi;
      permute(chi,  in._odata[SE->_offset], ptype);
    } else {
      chi_p = &in._odata[SE->_offset];
    }
  } else {
    chi_p = &buf[SE->_offset];
  }
  Impl::multLink(Uchi, U._odata[sU], *chi_p, dir, SE, st);
  vstream(out._odata[sF], Uchi);
  
};

template class NaiveStaggeredKernels<QCDpt::PStaggeredSmellImplF>;
template class NaiveStaggeredKernels<QCDpt::PStaggeredSmellImplD>;
template class NaiveStaggeredKernels<QCDpt::StaggeredSmellImplF>;
template class NaiveStaggeredKernels<QCDpt::StaggeredSmellImplD>;
template class NaiveStaggeredKernels<QCDpt::PStaggeredAdjointImplF>;
template class NaiveStaggeredKernels<QCDpt::PStaggeredAdjointImplD>;
template class NaiveStaggeredKernels<QCDpt::StaggeredAdjointImplF>;
template class NaiveStaggeredKernels<QCDpt::StaggeredAdjointImplD>;
template class NaiveStaggeredKernels<QCDpt::NaiveStaggeredImplF>;
template class NaiveStaggeredKernels<QCDpt::NaiveStaggeredImplD>;
template class NaiveStaggeredKernels<QCDpt::StaggeredAdjointImplNPF>;
template class NaiveStaggeredKernels<QCDpt::StaggeredAdjointImplNPD>;
}}

