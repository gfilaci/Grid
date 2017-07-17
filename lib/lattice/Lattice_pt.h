    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/lattice/Lattice_trace.h

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

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#ifndef GRID_LATTICE_PT_H
#define GRID_LATTICE_PT_H

namespace Grid {
    
  template<class vobj>
    inline auto Exponentiate(const Lattice<vobj> &lhs)
    -> Lattice<decltype(Exponentiate(lhs._odata[0]))>
    {
        Lattice<decltype(Exponentiate(lhs._odata[0]))> ret(lhs._grid);
        parallel_for(int ss=0;ss<lhs._grid->oSites();ss++){
            ret._odata[ss] = Exponentiate(lhs._odata[ss]);
        }
        return ret;
    }

template<class vobj>
    inline auto Logarithm(const Lattice<vobj> &lhs)
    -> Lattice<decltype(Logarithm(lhs._odata[0]))>
    {
        Lattice<decltype(Logarithm(lhs._odata[0]))> ret(lhs._grid);
        parallel_for(int ss=0;ss<lhs._grid->oSites();ss++){
            ret._odata[ss] = Logarithm(lhs._odata[ss]);
        }
        return ret;
    }

}
#endif

