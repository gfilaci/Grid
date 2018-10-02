    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/tensors/Tensor_functions_pt.h

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
#ifndef GRID_FUNCTIONS_PT_H
#define GRID_FUNCTIONS_PT_H

namespace Grid {

  /////////////////////////////////////////////// 
  // Exponentiate perturbative series
  ///////////////////////////////////////////////
   
template<class vtype> inline iScalar<vtype> Exponentiate(const iScalar<vtype>&r)
    {
      iScalar<vtype> ret;
      ret._internal = Exponentiate(r._internal);
      return ret;
    }

template<class vtype, int N> inline iVector<vtype, N> Exponentiate(const iVector<vtype,N>&r)
    {
      iVector<vtype, N> ret;
      for (int i = 0; i < N; i++)
        ret._internal[i] = Exponentiate(r._internal[i]);
      return ret;
    }

template<class vtype, int N> inline iMatrix<vtype, N> Exponentiate(const iMatrix<vtype,N>&r)
    {
      iMatrix<vtype, N> ret;
      for (int i = 0; i < N; i++)
          for (int j = 0; j < N; j++)
              ret._internal[i][j] = Exponentiate(r._internal[i][j]);
      return ret;
    }

template<class vtype, int N> inline iPert<vtype, N> Exponentiate(const iPert<vtype,N> &P)
  {
      // ASSUMING P(0) = 0
      iPert<vtype, N> ret(P), newtmp(P), tmp;
      
      typedef vtype mytype;
      mytype unit(1.0);
      
      ret._internal[0] = unit;
      
      for(int k=2; k<N; k++){
          zeroit(tmp);
          // i runs from 1 to N-1 (interval where P is defined)
          // but newtemp starts at order k-1
          // so there is a contribution only for i<N+1-k
          for(int i=1; i<N+1-k; i++){
              // j runs from k-1 to N-1 (interval where newtmp is defined)
              // but imposing i+j<N leads to j<N-i
              for(int j=k-1; j<N-i; j++){
                  // now k<=i+k<N
                  tmp._internal[i+j] += newtmp._internal[j] * P._internal[i];
              }
          }
          
          newtmp = (1./(double)k) * tmp;
          ret += newtmp;
      }
      
      return ret;
  }
  
  ///////////////////////////////////////////////
  // Logarithm of perturbative series
  ///////////////////////////////////////////////
  
  template<class vtype> inline iScalar<vtype> Logarithm(const iScalar<vtype>&r)
    {
      iScalar<vtype> ret;
      ret._internal = Logarithm(r._internal);
      return ret;
    }

template<class vtype, int N> inline iVector<vtype, N> Logarithm(const iVector<vtype,N>&r)
    {
      iVector<vtype, N> ret;
      for (int i = 0; i < N; i++)
        ret._internal[i] = Logarithm(r._internal[i]);
      return ret;
    }

template<class vtype, int N> inline iMatrix<vtype, N> Logarithm(const iMatrix<vtype,N>&r)
    {
      iMatrix<vtype, N> ret;
      for (int i = 0; i < N; i++)
          for (int j = 0; j < N; j++)
              ret._internal[i][j] = Logarithm(r._internal[i][j]);
      return ret;
    }

template<class vtype, int N> inline iPert<vtype, N> Logarithm(const iPert<vtype,N> &P)
  {
      // ASSUMING P(0) = 1
      iPert<vtype, N> ret(P), newtmp(P), tmp;
      double factor, sign = 1.;
      
      zeroit(ret._internal[0]);
      
      for(int k=2; k<N; k++){
          zeroit(tmp);
          // i runs from 1 to N-1 (interval where P is defined)
          // but newtemp starts at order k-1
          // so there is a contribution only for i<N+1-k
          for(int i=1; i<N+1-k; i++){
              // j runs from k-1 to N-1 (interval where newtmp is defined)
              // but imposing i+j<N leads to j<N-i
              for(int j=k-1; j<N-i; j++){
                  // now k<=i+k<N
                  tmp._internal[i+j] += newtmp._internal[j] * P._internal[i];
              }
          }
          
          newtmp = tmp;
          sign = -sign;
          factor = sign/(double)k;
          ret += factor * newtmp;
      }
      
      return ret;
  }

  ///////////////////////////////////////////////
  // Add to a given order
  ///////////////////////////////////////////////
  
  template<class vtype1, class vtype2> inline iScalar<vtype1> AddToOrd(const int &ord, const iScalar<vtype1> &r, const iScalar<vtype2> &Q)
    {
      iScalar<vtype1> ret;
      ret._internal = AddToOrd(ord,r._internal,Q._internal);
      return ret;
    }

template<class vtype1, class vtype2, int N> inline iVector<vtype1, N> AddToOrd(const int &ord, const iVector<vtype1,N> &r, const iScalar<vtype2> &Q)
    {
      iVector<vtype1, N> ret;
      for (int i = 0; i < N; i++)
        ret._internal[i] = AddToOrd(ord,r._internal[i],Q._internal);
      return ret;
    }

template<class vtype1, class vtype2, int N> inline iVector<vtype1, N> AddToOrd(const int &ord, const iVector<vtype1,N> &r, const iVector<vtype2,N> &Q)
    {
      iVector<vtype1, N> ret;
      for (int i = 0; i < N; i++)
        ret._internal[i] = AddToOrd(ord,r._internal[i],Q._internal[i]);
      return ret;
    }
    
template<class vtype1, class vtype2, int N> inline iMatrix<vtype1, N> AddToOrd(const int &ord, const iMatrix<vtype1,N> &r, const iScalar<vtype2> &Q)
    {
      iMatrix<vtype1, N> ret;
      for (int i = 0; i < N; i++)
          for (int j = 0; j < N; j++)
              ret._internal[i][j] = AddToOrd(ord,r._internal[i][j],Q._internal);
      return ret;
    }
    
template<class vtype, int N> inline iPert<vtype, N> AddToOrd(const int &ord, const iPert<vtype,N> &P, const iScalar<vtype> &Q)
  {
      iPert<vtype, N> ret(P);
      
      ret._internal[ord] += Q._internal;
      
      return ret;
  }
  
  ///////////////////////////////////////////////
  // Add to a given order (void)
  ///////////////////////////////////////////////
  
  template<class vtype1, class vtype2> inline void AddToOrdVoid(const int &ord, iScalar<vtype1> &r, const iScalar<vtype2> &Q)
    {
      AddToOrdVoid(ord,r._internal,Q._internal);
    }

template<class vtype1, class vtype2, int N> inline void AddToOrdVoid(const int &ord, iVector<vtype1,N> &r, const iScalar<vtype2> &Q)
    {
      for (int i = 0; i < N; i++)
        AddToOrdVoid(ord,r._internal[i],Q._internal);
    }

template<class vtype1, class vtype2, int N> inline void AddToOrdVoid(const int &ord, iVector<vtype1,N> &r, const iVector<vtype2,N> &Q)
    {
      for (int i = 0; i < N; i++)
        AddToOrdVoid(ord,r._internal[i],Q._internal[i]);
    }
    
template<class vtype1, class vtype2, int N> inline void AddToOrdVoid(const int &ord, iMatrix<vtype1,N> &r, const iScalar<vtype2> &Q)
    {
      for (int i = 0; i < N; i++)
          for (int j = 0; j < N; j++)
              AddToOrdVoid(ord,r._internal[i][j],Q._internal);
    }
    
template<class vtype, int N> inline void AddToOrdVoid(const int &ord, iPert<vtype,N> &P, const iScalar<vtype> &Q)
  {
      P._internal[ord] += Q._internal;
  }

  
template<class vtype1, class vtype2> inline void AddToOrdVoid(const int &ord, iScalar<vtype1> &r, const iScalar<vtype2> &Q, const RealD factor)
    {
      AddToOrdVoid(ord,r._internal,Q._internal,factor);
    }

template<class vtype1, class vtype2, int N> inline void AddToOrdVoid(const int &ord, iVector<vtype1,N> &r, const iScalar<vtype2> &Q, const RealD factor)
    {
      for (int i = 0; i < N; i++)
        AddToOrdVoid(ord,r._internal[i],Q._internal,factor);
    }

template<class vtype1, class vtype2, int N> inline void AddToOrdVoid(const int &ord, iVector<vtype1,N> &r, const iVector<vtype2,N> &Q, const RealD factor)
    {
      for (int i = 0; i < N; i++)
        AddToOrdVoid(ord,r._internal[i],Q._internal[i],factor);
    }
    
template<class vtype1, class vtype2, int N> inline void AddToOrdVoid(const int &ord, iMatrix<vtype1,N> &r, const iScalar<vtype2> &Q, const RealD factor)
    {
      for (int i = 0; i < N; i++)
          for (int j = 0; j < N; j++)
              AddToOrdVoid(ord,r._internal[i][j],Q._internal,factor);
    }
    
template<class vtype, int N> inline void AddToOrdVoid(const int &ord, iPert<vtype,N> &P, const iScalar<vtype> &Q, const RealD factor)
  {
      P._internal[ord] += Q._internal;
      P = factor * P;
  }
  
  ///////////////////////////////////////////////
  // Shifted sum
  ///////////////////////////////////////////////
  
  template<class vtype1, class vtype2> inline iScalar<vtype1> ShiftedSum(const int &ord, const iScalar<vtype1> &r, const iScalar<vtype2> &Q)
    {
      iScalar<vtype1> ret;
      ret._internal = ShiftedSum(ord,r._internal,Q._internal);
      return ret;
    }

template<class vtype1, class vtype2, int N> inline iVector<vtype1, N> ShiftedSum(const int &ord, const iVector<vtype1,N> &r, const iScalar<vtype2> &Q)
    {
      iVector<vtype1, N> ret;
      for (int i = 0; i < N; i++)
        ret._internal[i] = ShiftedSum(ord,r._internal[i],Q._internal);
      return ret;
    }

template<class vtype1, class vtype2, int N> inline iVector<vtype1, N> ShiftedSum(const int &ord, const iVector<vtype1,N> &r, const iVector<vtype2,N> &Q)
    {
      iVector<vtype1, N> ret;
      for (int i = 0; i < N; i++)
        ret._internal[i] = ShiftedSum(ord,r._internal[i],Q._internal[i]);
      return ret;
    }
    
template<class vtype1, class vtype2, int N> inline iMatrix<vtype1, N> ShiftedSum(const int &ord, const iMatrix<vtype1,N> &r, const iScalar<vtype2> &Q)
    {
      iMatrix<vtype1, N> ret;
      for (int i = 0; i < N; i++)
          for (int j = 0; j < N; j++)
              ret._internal[i][j] = ShiftedSum(ord,r._internal[i][j],Q._internal);
      return ret;
    }
    
template<class vtype, int N> inline iPert<vtype, N> ShiftedSum(const int &ord, const iPert<vtype,N> &P, const iPert<vtype, N> &Q)
  {
      iPert<vtype, N> ret(P);
      
      for (int i=ord; i<N; i++) {
          ret._internal[i] += Q._internal[i-ord];
      }
      return ret;
  }

  ///////////////////////////////////////////////
  // Shifted sum (void)
  ///////////////////////////////////////////////
  
template<class vtype1, class vtype2> inline void ShiftedSumVoid(const int &ord, iScalar<vtype1> &r, const iScalar<vtype2> &Q, const RealD factor)
    {
        ShiftedSumVoid(ord,r._internal,Q._internal,factor);
    }

template<class vtype1, class vtype2, int N> inline void ShiftedSumVoid(const int &ord, iVector<vtype1,N> &r, const iScalar<vtype2> &Q, const RealD factor)
    {
      for (int i = 0; i < N; i++)
        ShiftedSumVoid(ord,r._internal[i],Q._internal,factor);
    }

template<class vtype1, class vtype2, int N> inline void ShiftedSumVoid(const int &ord, iVector<vtype1,N> &r, const iVector<vtype2,N> &Q, const RealD factor)
    {
      for (int i = 0; i < N; i++)
        ShiftedSumVoid(ord,r._internal[i],Q._internal[i],factor);
    }
    
template<class vtype1, class vtype2, int N> inline void ShiftedSumVoid(const int &ord, iMatrix<vtype1,N> &r, const iScalar<vtype2> &Q, const RealD factor)
    {
      for (int i = 0; i < N; i++)
          for (int j = 0; j < N; j++)
              ShiftedSumVoid(ord,r._internal[i][j],Q._internal,factor);
    }
    
template<class vtype, int N> inline void ShiftedSumVoid(const int &ord, iPert<vtype,N> &P, const iPert<vtype, N> &Q, const RealD factor)
  {
      for (int i=ord; i<N; i++) {
          P._internal[i] += factor * Q._internal[i-ord];
      }
  }


/************************************************/
/*           QUAD PRECISION FUNCTIONS           */
/************************************************/
#ifdef USE_QUADPREC
////////////
// overloading for exp in quad precision
///////////
template<int N, int M> inline iPert<iMatrix<vComplexD,M>,N> Exponentiate(const iPert<iMatrix<vComplexD,M>,N> &P)
{
	typedef iMatrix<   vComplexD,M> vtype;
	typedef iMatrix<    ComplexD,M> stype;
	typedef iMatrix<__complex128,M> qtype;
	
	// extract simd vector
	int Nsimd = sizeof(vtype::vector_type) / sizeof(vtype::scalar_type);
	std::vector<iPert<stype,N>> dbuf(Nsimd);
	std::vector<iPert<qtype,N>> qbuf(Nsimd);
	extract(P,dbuf);
	
	// cast to quad precision
	for(int i=0; i<Nsimd; i++){
		for(int j=0; j<N; j++){
			for(int l=0; l<M; l++){
				for(int m=0; m<M; m++){
					__real__ qbuf[i](j)(l,m) = dbuf[i](j)(l,m).real();
					__imag__ qbuf[i](j)(l,m) = dbuf[i](j)(l,m).imag();
				}
			}
		}
	}
	
	// exponentiate
	for(int z=0; z<Nsimd; z++) qbuf[z] = QuadExponentiate(qbuf[z]);
	
	// cast to double precision
	for(int i=0; i<Nsimd; i++){
		for(int j=0; j<N; j++){
			for(int l=0; l<M; l++){
				for(int m=0; m<M; m++){
					dbuf[i](j)(l,m) = ComplexD(__real__ qbuf[i](j)(l,m), __imag__ qbuf[i](j)(l,m));
				}
			}
		}
	}
	
	// merge simd vector
	iPert<vtype, N> result;
	merge(result,dbuf);
	
	return result;
}

 template<int N, int M> inline iPert<iMatrix<__complex128,M>,N> QuadExponentiate(const iPert<iMatrix<__complex128,M>,N> &qbuf)
{
	typedef iMatrix<__complex128,M> qtype;
	auto outbuf = qbuf;
	
	iPert<qtype, N> newtmp, tmp;
	qtype unit(1.0);
			
		newtmp = outbuf;
		
		// ASSUMING P(0) = 0
		outbuf._internal[0] = unit;
		
		for(int k=2; k<N; k++){
		        tmp = zero;
			// i runs from 1 to N-1 (interval where P is defined)
			// but newtemp starts at order k-1
			// so there is a contribution only for i<N+1-k
			for(int i=1; i<N+1-k; i++){
				// j runs from k-1 to N-1 (interval where newtmp is defined)
				// but imposing i+j<N leads to j<N-i
				for(int j=k-1; j<N-i; j++){
					// now k<=i+k<N
					tmp._internal[i+j] += newtmp._internal[j] * qbuf._internal[i];
				}
			}
			
			newtmp = (1./(__float128)k) * tmp;
			outbuf += newtmp;
		}
	
	return outbuf;
}

////////////
// overloading for log in quad precision
///////////
template<int N, int M> inline iPert<iMatrix<vComplexD,M>,N> Logarithm(const iPert<iMatrix<vComplexD,M>,N> &P)
{
	typedef iMatrix<   vComplexD,M> vtype;
	typedef iMatrix<    ComplexD,M> stype;
	typedef iMatrix<__complex128,M> qtype;
	
	// extract simd vector
	int Nsimd = sizeof(vtype::vector_type) / sizeof(vtype::scalar_type);
	std::vector<iPert<stype,N>> dbuf(Nsimd);
	std::vector<iPert<qtype,N>> qbuf(Nsimd);
	extract(P,dbuf);
	
	// cast to quad precision
	for(int i=0; i<Nsimd; i++){
		for(int j=0; j<N; j++){
			for(int l=0; l<M; l++){
				for(int m=0; m<M; m++){
					__real__ qbuf[i](j)(l,m) = dbuf[i](j)(l,m).real();
					__imag__ qbuf[i](j)(l,m) = dbuf[i](j)(l,m).imag();
				}
			}
		}
	}
	
	// exponentiate
	for(int z=0; z<Nsimd; z++) qbuf[z] = QuadLogarithm(qbuf[z]);
	
	// cast to double precision
	for(int i=0; i<Nsimd; i++){
		for(int j=0; j<N; j++){
			for(int l=0; l<M; l++){
				for(int m=0; m<M; m++){
					dbuf[i](j)(l,m) = ComplexD(__real__ qbuf[i](j)(l,m), __imag__ qbuf[i](j)(l,m));
				}
			}
		}
	}
	
	// merge simd vector
	iPert<vtype, N> result;
	merge(result,dbuf);
	
	return result;
}

 template<int N, int M> inline iPert<iMatrix<__complex128,M>,N> QuadLogarithm(const iPert<iMatrix<__complex128,M>,N> &qbuf)
{
	typedef iMatrix<__complex128,M> qtype;
	auto outbuf = qbuf;
	
	iPert<qtype, N> newtmp, tmp;
	__float128 factor, sign;
		
		newtmp = outbuf;
		sign = 1.;
		
		// ASSUMING P(0) = 1
		outbuf._internal[0] = zero;
		
		for(int k=2; k<N; k++){
		        tmp = zero;
			// i runs from 1 to N-1 (interval where P is defined)
			// but newtemp starts at order k-1
			// so there is a contribution only for i<N+1-k
			for(int i=1; i<N+1-k; i++){
				// j runs from k-1 to N-1 (interval where newtmp is defined)
				// but imposing i+j<N leads to j<N-i
				for(int j=k-1; j<N-i; j++){
					// now k<=i+k<N
					tmp._internal[i+j] += newtmp._internal[j] * qbuf._internal[i];
				}
			}
			
			newtmp = tmp;
			sign = -sign;
			factor = sign/(__float128)k;
			outbuf += factor * newtmp;
		}
	
	return outbuf;
}
#endif

}

#endif
