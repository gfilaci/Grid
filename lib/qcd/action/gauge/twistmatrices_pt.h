/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/utils/twistmatrices_pt.h

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
#ifndef QCD_UTIL_TWISTMATRICES_PT_H
#define QCD_UTIL_TWISTMATRICES_PT_H

namespace Grid {
namespace QCD {
namespace QCDpt {
  
template <int ncolour>
class twistmatrices {

public:
    template <typename vtype>
    using iTwistMatrix = iScalar<iScalar<iScalar<iMatrix<vtype, ncolour> > > >;
    
    typedef iTwistMatrix<vComplex>  vMatrix;
    typedef iTwistMatrix<vComplexF> vMatrixF;
    typedef iTwistMatrix<vComplexD> vMatrixD;
    
private:
    vMatrix omega1, omega2, identity;
    vMatrix omega[Nd], adjomega[Nd];
    
public:
  
  explicit twistmatrices()
  {
      initialisetwist();

#define istwisted(mu) (mu==0 || mu==1)
      
      omega[0] = omega1;
      omega[1] = omega2;
      omega[2] = identity;
      omega[3] = identity;

//      omega[0] = identity;
//      omega[1] = omega1;
//      omega[2] = omega2;
//      omega[3] = identity;
      
      for (int i=0; i<Nd; i++) {
          adjomega[i] = adj(omega[i]);
      }
      
  }
  
  void initialisetwist(){
      std::cout<<"Need to define twist matrices for Nc = "<<ncolour<<std::endl;
      exit(EXIT_FAILURE);
  }
  
  // manually taking the first complex object from SIMD vector... (portable?)
  iMatrix<Complex,Nd> twist_tensor(){
      iMatrix<Complex,Nd> eta;
      zeroit(eta);
      Complex im(0.,1.);
      for (int mu=0; mu<Nd; mu++) {
          for (int nu=0; nu<Nd; nu++) {
              eta(mu,nu) = TensorRemove(trace(omega[nu]*omega[mu]*adjomega[nu]*adjomega[mu])).v[0] / (double)ncolour;
              eta(mu,nu) += im * TensorRemove(trace(omega[nu]*omega[mu]*adj(omega[nu])*adj(omega[mu]))).v[1] / (double)ncolour;
              eta(mu,nu) = -im * (double)ncolour / (2*M_PI) * std::log(eta(mu,nu));
          }
      }
      return eta;
  }
  
  template<class T>
  T forward(const T &field, const int &mu)
  {
      return omega[mu] * field * adjomega[mu];
  }
  template<class T>
  T backward(const T &field, const int &mu)
  {
      return adjomega[mu] * field * omega[mu];
  }
    
    
};


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
#endif
