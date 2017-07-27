/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/utils/twistmatrices_pt.cc

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

#include <Grid.h>
#include "twistmatrices_pt.h"

namespace Grid {
namespace QCD {
namespace QCDpt {
    
// specialisation for SU(2) twist matrices
template<>
void twistmatrices<2>::initialisetwist(){
    Complex im(0.,1.);
    
    identity()()()(0,0) = 1.;
    identity()()()(0,1) = 0.;
    identity()()()(1,0) = 0.;
    identity()()()(1,1) = 1.;
    
    omega1()()()(0,0) = -im;
    omega1()()()(0,1) = 0.;
    omega1()()()(1,0) = 0.;
    omega1()()()(1,1) = im;
    
    omega2()()()(0,0) = 0.;
    omega2()()()(0,1) = 1.;
    omega2()()()(1,0) = -1.;
    omega2()()()(1,1) = 0.;
}

// specialisation for SU(3) twist matrices
template<>
void twistmatrices<3>::initialisetwist(){
    Complex im(0.,1.);
    double tmp = 2. * M_PI / 3.;
    double tmpcos = std::cos(tmp);
    double tmpsin = std::sin(tmp);
    
    identity()()()(0,0) = 1.;
    identity()()()(0,1) = 0.;
    identity()()()(0,2) = 0.;
    identity()()()(1,0) = 0.;
    identity()()()(1,1) = 1.;
    identity()()()(1,2) = 0.;
    identity()()()(2,0) = 0.;
    identity()()()(2,1) = 0.;
    identity()()()(2,2) = 1.;
    
    omega1()()()(0,0) = tmpcos - im*tmpsin;
    omega1()()()(0,1) = 0.;
    omega1()()()(0,2) = 0.;
    omega1()()()(1,0) = 0.;
    omega1()()()(1,1) = 1.;
    omega1()()()(1,2) = 0.;
    omega1()()()(2,0) = 0.;
    omega1()()()(2,1) = 0.;
    omega1()()()(2,2) = tmpcos + im*tmpsin;
    
    omega2()()()(0,0) = 0.;
    omega2()()()(0,1) = 1.;
    omega2()()()(0,2) = 0.;
    omega2()()()(1,0) = 0.;
    omega2()()()(1,1) = 0.;
    omega2()()()(1,2) = 1.;
    omega2()()()(2,0) = 1.;
    omega2()()()(2,1) = 0.;
    omega2()()()(2,2) = 0.;
}
    
    
}
}
}
