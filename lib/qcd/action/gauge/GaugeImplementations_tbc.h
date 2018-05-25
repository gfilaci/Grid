/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/gauge/GaugeImplementations_tbc.h

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
#ifndef GRID_QCD_GAUGE_IMPLEMENTATIONS_TBC_H
#define GRID_QCD_GAUGE_IMPLEMENTATIONS_TBC_H

#include "GaugeImplTypes.h"
#include "GaugeImplTypes_pt.h"
#include "twistmatrices_pt.h"

namespace Grid {
namespace QCD {
namespace QCDpt {


template <class GimplTypes> class TwistedGaugeImpl : public GimplTypes {
public:
  INHERIT_GIMPL_TYPES(GimplTypes);
  
  static twistmatrices<Nc,typename Simd::scalar_type> twist;
  
  template <class covariant>
  static Lattice<covariant> CovShiftForward(const GaugeLinkField &Link, int mu,
                                            const Lattice<covariant> &field) {
      GridBase * grid = Link._grid;
      
      int Lmu = grid->GlobalDimensions()[mu]-1;
      
      conformable(field,Link);
      
      Lattice<iScalar<vInteger> > coor(grid);    LatticeCoordinate(coor,mu);
      
      Lattice<covariant> field_bc = Cshift(field,mu,1);
      
      if(istwisted(mu))
          field_bc = where(coor==Lmu,twist.forward(field_bc,mu),field_bc);
      
      return Link*field_bc;
  }

  template <class covariant>
  static Lattice<covariant> CovShiftBackward(const GaugeLinkField &Link, int mu,
                                             const Lattice<covariant> &field) {
    GridBase * grid = field._grid;

    int Lmu = grid->GlobalDimensions()[mu]-1;

    conformable(field,Link);

    Lattice<iScalar<vInteger> > coor(grid);    LatticeCoordinate(coor,mu);

    Lattice<covariant> tmp(grid);

    tmp = adj(Link)*field;
    
    if(istwisted(mu))
    tmp = where(coor==Lmu,twist.backward(tmp,mu),tmp);
    
    return Cshift(tmp,mu,-1);
  }

  static inline GaugeLinkField
  CovShiftIdentityBackward(const GaugeLinkField &Link, int mu) {
    GridBase *grid = Link._grid;
    int Lmu = grid->GlobalDimensions()[mu] - 1;

    Lattice<iScalar<vInteger>> coor(grid);
    LatticeCoordinate(coor, mu);
    
    GaugeLinkField tmp(grid);
    tmp = adj(Link);
    
    if(istwisted(mu))
    tmp = where(coor == Lmu, twist.backward(tmp,mu), tmp);
    
    return Cshift(tmp, mu, -1);
  }
  static inline GaugeLinkField
  CovShiftIdentityForward(const GaugeLinkField &Link, int mu) {
    return Link;
  }

  static inline GaugeLinkField ShiftStaple(const GaugeLinkField &Link, int mu) {
    GridBase *grid = Link._grid;
    int Lmu = grid->GlobalDimensions()[mu] - 1;

    Lattice<iScalar<vInteger>> coor(grid);
    LatticeCoordinate(coor, mu);

    GaugeLinkField tmp(grid);
    tmp = Cshift(Link, mu, 1);
    
    if(istwisted(mu))
    tmp = where(coor == Lmu, twist.forward(tmp,mu), tmp);
    
    return tmp;
  }
  
  static inline GaugeLinkField MoveForward(const GaugeLinkField &Link, int mu) {
    return ShiftStaple(Link,mu);
  }
  
  static inline GaugeLinkField MoveBackward(const GaugeLinkField &Link, int mu) {
    GridBase *grid = Link._grid;
    int Lmu = grid->GlobalDimensions()[mu] - 1;

    Lattice<iScalar<vInteger>> coor(grid);
    LatticeCoordinate(coor, mu);

    GaugeLinkField tmp = Link;
    
    if(istwisted(mu))
    tmp = where(coor == Lmu, twist.backward(tmp,mu), tmp);
    
    return Cshift(tmp, mu, -1);
  }
  
  static inline bool isPeriodicGaugeField(void) { return false; }
};

    
template <class GimplTypes> class TwistedGaugeImplNP : public GimplTypes {
public:
    INHERIT_GIMPL_TYPES(GimplTypes);
    
    static twistmatricesNP<Nc,typename Simd::scalar_type> twist;
    
    template <class covariant>
    static Lattice<covariant> CovShiftForward(const GaugeLinkField &Link, int mu,
                                              const Lattice<covariant> &field) {
        GridBase * grid = Link._grid;
        
        int Lmu = grid->GlobalDimensions()[mu]-1;
        
        conformable(field,Link);
        
        Lattice<iScalar<vInteger> > coor(grid);    LatticeCoordinate(coor,mu);
        
        Lattice<covariant> field_bc = Cshift(field,mu,1);
        
        if(istwisted(mu))
            field_bc = where(coor==Lmu,twist.forward(field_bc,mu),field_bc);
        
        return Link*field_bc;
    }
    
    template <class covariant>
    static Lattice<covariant> CovShiftBackward(const GaugeLinkField &Link, int mu,
                                               const Lattice<covariant> &field) {
        GridBase * grid = field._grid;
        
        int Lmu = grid->GlobalDimensions()[mu]-1;
        
        conformable(field,Link);
        
        Lattice<iScalar<vInteger> > coor(grid);    LatticeCoordinate(coor,mu);
        
        Lattice<covariant> tmp(grid);
        
        tmp = adj(Link)*field;
        
        if(istwisted(mu))
            tmp = where(coor==Lmu,twist.backward(tmp,mu),tmp);
        
        return Cshift(tmp,mu,-1);
    }
    
    static inline GaugeLinkField
    CovShiftIdentityBackward(const GaugeLinkField &Link, int mu) {
        GridBase *grid = Link._grid;
        int Lmu = grid->GlobalDimensions()[mu] - 1;
        
        Lattice<iScalar<vInteger>> coor(grid);
        LatticeCoordinate(coor, mu);
        
        GaugeLinkField tmp(grid);
        tmp = adj(Link);
        
        if(istwisted(mu))
            tmp = where(coor == Lmu, twist.backward(tmp,mu), tmp);
        
        return Cshift(tmp, mu, -1);
    }
    static inline GaugeLinkField
    CovShiftIdentityForward(const GaugeLinkField &Link, int mu) {
        return Link;
    }
    
    static inline GaugeLinkField ShiftStaple(const GaugeLinkField &Link, int mu) {
        GridBase *grid = Link._grid;
        int Lmu = grid->GlobalDimensions()[mu] - 1;
        
        Lattice<iScalar<vInteger>> coor(grid);
        LatticeCoordinate(coor, mu);
        
        GaugeLinkField tmp(grid);
        tmp = Cshift(Link, mu, 1);
        
        if(istwisted(mu))
            tmp = where(coor == Lmu, twist.forward(tmp,mu), tmp);
        
        return tmp;
    }
    
    static inline GaugeLinkField MoveForward(const GaugeLinkField &Link, int mu) {
        return ShiftStaple(Link,mu);
    }
    
    static inline GaugeLinkField MoveBackward(const GaugeLinkField &Link, int mu) {
        GridBase *grid = Link._grid;
        int Lmu = grid->GlobalDimensions()[mu] - 1;
        
        Lattice<iScalar<vInteger>> coor(grid);
        LatticeCoordinate(coor, mu);
        
        GaugeLinkField tmp = Link;
        
        if(istwisted(mu))
            tmp = where(coor == Lmu, twist.backward(tmp,mu), tmp);
        
        return Cshift(tmp, mu, -1);
    }
    
    static inline bool isPeriodicGaugeField(void) { return false; }
};
    
template <class GimplTypes>
twistmatrices<Nc,typename GimplTypes::Simd::scalar_type> TwistedGaugeImpl<GimplTypes>::twist;

template <class GimplTypes>
twistmatricesNP<Nc,typename GimplTypes::Simd::scalar_type> TwistedGaugeImplNP<GimplTypes>::twist;
    
typedef TwistedGaugeImpl<GimplTypes_ptR> TwistedGimpl_ptR;
typedef PeriodicGaugeImpl<GimplTypes_ptR> PeriodicGimpl_ptR;

}
}
}

#endif
