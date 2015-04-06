#ifndef GRID_SIMD_H
#define GRID_SIMD_H

////////////////////////////////////////////////////////////////////////
// Define scalar and vector floating point types
//
// Scalar:   RealF, RealD, ComplexF, ComplexD
//
// Vector:  vRealF, vRealD, vComplexF, vComplexD
//
// Vector types are arch dependent
////////////////////////////////////////////////////////////////////////
    // TODO
    //
    // Base class to share common code between vRealF, VComplexF etc...
    //
    // lattice Broad cast assignment
    //
    // where() support
    // implement with masks, and/or? Type of the mask & boolean support?
    //
    // Unary functions
    // cos,sin, tan, acos, asin, cosh, acosh, tanh, sinh, // Scalar<vReal> only arg
    // exp, log, sqrt, fabs
    //
    // transposeColor, transposeSpin,
    // adjColor, adjSpin,
    // traceColor, traceSpin.
    // peekColor, peekSpin + pokeColor PokeSpin
    //
    // copyMask.
    //
    // localMaxAbs
    //
    // norm2,
    // sumMulti equivalent.
    // Fourier transform equivalent.
    //
    
  ////////////////////////////////////////////////////////////
  // SIMD Alignment controls
  ////////////////////////////////////////////////////////////
#ifdef HAVE_VAR_ATTRIBUTE_ALIGNED
#define ALIGN_DIRECTIVE(A) __attribute__ ((aligned(A)))
#else
#define ALIGN_DIRECTIVE(A) __declspec(align(A))
#endif

#ifdef SSE2
#include <pmmintrin.h>
#define SIMDalign ALIGN_DIRECTIVE(16)
#endif

#if defined(AVX1) || defined (AVX2)
#include <immintrin.h>
#define SIMDalign ALIGN_DIRECTIVE(32)
#endif

#ifdef AVX512
#include <immintrin.h>
#define SIMDalign ALIGN_DIRECTIVE(64)
#endif

namespace Grid {

  typedef  float  RealF;
  typedef  double RealD;
  typedef  RealF  Real;
  
  typedef std::complex<RealF> ComplexF;
  typedef std::complex<RealD> ComplexD;
  typedef std::complex<Real>  Complex;




  inline RealF adj(const RealF  & r){ return r; }
  inline RealF conj(const RealF  & r){ return r; }
  inline ComplexD localInnerProduct(const ComplexD & l, const ComplexD & r) { return conj(l)*r; }
  inline ComplexF localInnerProduct(const ComplexF & l, const ComplexF & r) { return conj(l)*r; }
  inline RealD localInnerProduct(const RealD & l, const RealD & r) { return l*r; }
  inline RealF localInnerProduct(const RealF & l, const RealF & r) { return l*r; }

    ////////////////////////////////////////////////////////////////////////////////
    //Provide support functions for basic real and complex data types required by Grid
    //Single and double precision versions. Should be able to template this once only.
    ////////////////////////////////////////////////////////////////////////////////
    inline void mac (ComplexD * __restrict__ y,const ComplexD * __restrict__ a,const ComplexD *__restrict__ x){ *y = (*a) * (*x)+(*y); };
    inline void mult(ComplexD * __restrict__ y,const ComplexD * __restrict__ l,const ComplexD *__restrict__ r){ *y = (*l) * (*r);}
    inline void sub (ComplexD * __restrict__ y,const ComplexD * __restrict__ l,const ComplexD *__restrict__ r){ *y = (*l) - (*r);}
    inline void add (ComplexD * __restrict__ y,const ComplexD * __restrict__ l,const ComplexD *__restrict__ r){ *y = (*l) + (*r);}
    inline ComplexD adj(const ComplexD& r){ return(conj(r)); }
    // conj already supported for complex
    
    inline void mac (ComplexF * __restrict__ y,const ComplexF * __restrict__ a,const ComplexF *__restrict__ x){ *y = (*a) * (*x)+(*y); }
    inline void mult(ComplexF * __restrict__ y,const ComplexF * __restrict__ l,const ComplexF *__restrict__ r){ *y = (*l) * (*r); }
    inline void sub (ComplexF * __restrict__ y,const ComplexF * __restrict__ l,const ComplexF *__restrict__ r){ *y = (*l) - (*r); }
    inline void add (ComplexF * __restrict__ y,const ComplexF * __restrict__ l,const ComplexF *__restrict__ r){ *y = (*l) + (*r); }
    inline Complex  adj(const Complex& r ){ return(conj(r)); }
    //conj already supported for complex
    
    inline void mac (RealD * __restrict__ y,const RealD * __restrict__ a,const RealD *__restrict__ x){  *y = (*a) * (*x)+(*y);}
    inline void mult(RealD * __restrict__ y,const RealD * __restrict__ l,const RealD *__restrict__ r){ *y = (*l) * (*r);}
    inline void sub (RealD * __restrict__ y,const RealD * __restrict__ l,const RealD *__restrict__ r){ *y = (*l) - (*r);}
    inline void add (RealD * __restrict__ y,const RealD * __restrict__ l,const RealD *__restrict__ r){ *y = (*l) + (*r);}
    inline RealD adj(const RealD & r){ return r; }  // No-op for real
    inline RealD conj(const RealD & r){ return r; }
    
    inline void mac (RealF * __restrict__ y,const RealF * __restrict__ a,const RealF *__restrict__ x){  *y = (*a) * (*x)+(*y); }
    inline void mult(RealF * __restrict__ y,const RealF * __restrict__ l,const RealF *__restrict__ r){ *y = (*l) * (*r); }
    inline void sub (RealF * __restrict__ y,const RealF * __restrict__ l,const RealF *__restrict__ r){ *y = (*l) - (*r); }
    inline void add (RealF * __restrict__ y,const RealF * __restrict__ l,const RealF *__restrict__ r){ *y = (*l) + (*r); }
    


  class Zero{};
  static Zero zero;
  template<class itype> inline void ZeroIt(itype &arg){ arg=zero;};
  template<>            inline void ZeroIt(ComplexF &arg){ arg=0; };
  template<>            inline void ZeroIt(ComplexD &arg){ arg=0; };
  template<>            inline void ZeroIt(RealF &arg){ arg=0; };
  template<>            inline void ZeroIt(RealD &arg){ arg=0; };



#if defined (SSE2)
    typedef __m128 fvec;
    typedef __m128d dvec;
    typedef __m128 cvec;
    typedef __m128d zvec;
    typedef __m128i ivec;
#endif
#if defined (AVX1) || defined (AVX2)
    typedef __m256 fvec;
    typedef __m256d dvec;
    typedef __m256  cvec;
    typedef __m256d zvec;
    typedef __m256i ivec;
#endif
#if defined (AVX512)
    typedef __m512  fvec;
    typedef __m512d dvec;
    typedef __m512  cvec;
    typedef __m512d zvec;
    typedef __m512i ivec;
#endif
#if defined (QPX)
    typedef float  fvec __attribute__ ((vector_size (16))); // QPX has same SIMD width irrespective of precision
    typedef float  cvec __attribute__ ((vector_size (16)));
    
    typedef vector4double dvec;
    typedef vector4double zvec;
#endif
#if defined (AVX1) || defined (AVX2) || defined (AVX512)
    inline void v_prefetch0(int size, const char *ptr){
          for(int i=0;i<size;i+=64){ //  Define L1 linesize above// What about SSE?
            _mm_prefetch(ptr+i+4096,_MM_HINT_T1);
            _mm_prefetch(ptr+i+512,_MM_HINT_T0);
          }
    }
#else 
    inline void v_prefetch0(int size, const char *ptr){};
#endif

};


/////////////////////////////////////////////////////////////////
// Generic extract/merge/permute
/////////////////////////////////////////////////////////////////
template<class vsimd,class scalar,int Nsimd>
inline void Gextract(vsimd &y,std::vector<scalar *> &extracted){
  // Bounce off stack is painful
  // temporary hack while I figure out the right interface
  scalar buf[Nsimd]; 
  vstore(y,buf);
  for(int i=0;i<Nsimd;i++){
    *extracted[i] = buf[i];
    extracted[i]++;
  }
};
template<class vsimd,class scalar,int Nsimd>
inline void Gmerge(vsimd &y,std::vector<scalar *> &extracted){
  scalar buf[Nsimd]; 
  for(int i=0;i<Nsimd;i++){
    buf[i]=*extracted[i];
    extracted[i]++;
  }
  vset(y,buf); 
};

//////////////////////////////////////////////////////////
// Permute
// Permute 0 every ABCDEFGH -> BA DC FE HG
// Permute 1 every ABCDEFGH -> CD AB GH EF
// Permute 2 every ABCDEFGH -> EFGH ABCD
// Permute 3 possible on longer iVector lengths (512bit = 8 double = 16 single)
// Permute 4 possible on half precision @512bit vectors.
//////////////////////////////////////////////////////////
// Should be able to make the permute/extract/merge independent of the
// vector subtype and reduce the volume of code.
template<class vsimd>
inline void Gpermute(vsimd &y,vsimd b,int perm){
      switch (perm){
#if defined(AVX1)||defined(AVX2)
      // 8x32 bits=>3 permutes
      case 2: y.v = _mm256_shuffle_ps(b.v,b.v,_MM_SHUFFLE(2,3,0,1)); break;
      case 1: y.v = _mm256_shuffle_ps(b.v,b.v,_MM_SHUFFLE(1,0,3,2)); break;
      case 0: y.v = _mm256_permute2f128_ps(b.v,b.v,0x01); break;
#endif
#ifdef SSE2
      case 1: y.v = _mm_shuffle_ps(b.v,b.v,_MM_SHUFFLE(2,3,0,1)); break;
      case 0: y.v = _mm_shuffle_ps(b.v,b.v,_MM_SHUFFLE(1,0,3,2));break;
#endif
#ifdef AVX512
	// 16 floats=> permutes
        // Permute 0 every abcd efgh ijkl mnop -> badc fehg jilk nmpo 
        // Permute 1 every abcd efgh ijkl mnop -> cdab ghef jkij opmn 
        // Permute 2 every abcd efgh ijkl mnop -> efgh abcd mnop ijkl
        // Permute 3 every abcd efgh ijkl mnop -> ijkl mnop abcd efgh
      case 3: y.v = _mm512_swizzle_ps(b.v,_MM_SWIZ_REG_CDAB); break;
      case 2: y.v = _mm512_swizzle_ps(b.v,_MM_SWIZ_REG_BADC); break;
      case 1: y.v = _mm512_permute4f128_ps(b.v,(_MM_PERM_ENUM)_MM_SHUFFLE(2,3,0,1)); break;
      case 0: y.v = _mm512_permute4f128_ps(b.v,(_MM_PERM_ENUM)_MM_SHUFFLE(1,0,3,2)); break;
#endif
#ifdef QPX
#error not implemented
#endif
      default: assert(0); break;
      }
    };

#include <Grid_vRealF.h>
#include <Grid_vRealD.h>
#include <Grid_vComplexF.h>
#include <Grid_vComplexD.h>
#include <Grid_vInteger.h>

#endif