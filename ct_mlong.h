#ifndef __BIPOLY_CT_MLONG
#define __BIPOLY_CT_MLONG 1

#define HAVE_CT 1
#define MODULAR 1

#if 1
// Smallest prime under 2^32
#define P 4294967291U
#define Pdef 5			// = 2^32-P
#define P_digits 10		// # of decimal digits in P (for atomodp)
#define P_power10 1000000000	// largest power of 10 less than P (for atomodp)

// Next smallest:
// P:    4294967279, 4294967231, 4294967197, 4294967189
// Pdef: 17,         65,         99,         107
#else
extern unsigned int P;			// Modulus
extern unsigned int Pdef;		// 2^w-P, where w is the number of bits in an unsigned int
#endif

// Data type for the polynomial coefficients.
typedef unsigned int ct;

// Arithmetic
//  These make a number of assumptions
//   #1: int is 32 bits
//   #2: long long is 64 bits
//   #3: P is nearly 2^32 (at least we assume P > (2^32 - 2^16))


#define _hi(a) (a >> 32)
#define _lo(a) (a & 0xFFFFFFFFULL)

#define __m_normalize(a) ((unsigned int) ((((unsigned int) a) >= P) ? ((a)-P) : (a)))
#define __m_normalize_general(ll) ((unsigned long long) (_lo(ll) + Pdef * _hi(ll)))

// a+b
#define A(a,b) (__m_normalize((unsigned long long) (a) + (unsigned long long) (b)))
// a-b
#define S(a,b) (A(a,P-b))
// a*b
#define M(a,b) (__m_normalize(__m_normalize_general(__m_normalize_general((unsigned long long) (a) * (unsigned long long) (b)))))
/*

// Above #define is equivalent to this easier-to-read version,
//  but optimizes better (for me)

inline ct M(ct a, ct b)
{
    unsigned long long m;
    ct z;
    
    m = (unsigned long long) (a) * (unsigned long long) (b);
    m = __m_normalize_general(m);
    m = __m_normalize_general(m);
    z = __m_normalize(m);
    return z;
}
*/
// a += b
#define AE(a,b) (a = A(a,b))
// a -= b
#define SE(a,b) (a = S(a,b))
// a *= b
#define ME(a,b) (a = M(a,b))

// a * b, where b is an integer, assumed to be +/-1
#define MS(a,b) ((b == -1) ? ((P-a) % P) : a)
//#define MS(a,b) ((a * b) + ((((int) b) >> 1) & P))

// a *= b, where b is as above
#define MSE(a,b) (a = MS(a,b))
    
// Useful constant values, of type ct.
#define COEFF_0 0
#define COEFF_NEG_1 (P-1)
#define COEFF_1 1

// Printf and scanf format characters for a signed and unsigned ct
#define CTF "%u"
#define CTFS "%u"
#define CTFX "%x"
#define CTFSX "%x"

// How to pass a ct to printf/scanf
#define I(ct) (&ct)
#define O(ct) (ct)

// String conversion for ct
//  (the function atomodp() is defined in bipoly.c)
#define atoct(s) atomodp(s)

// Largest positive value representable by a ct. (not applicable here)
//#define COEFF_MAX LONG_MAX
#define HAVE_COEFF_MAX 0

// Conversions to/from ct

#define from_i(i) (__m_normalize(i))

#define to_d(a) ((double) a)
#define to_i(a) ((int) a)
#define to_ht(a) ((ht) a)			// ht is a long long
#define sign(a) ((int) a > 0 ? 1 : -1)

// Misc
//inline ct ct_abs(ct a) { return ((int) a > 0) ? a : P-a; }
#define ct_abs(a) (a)
#define ct_sign(a) (1)

#endif
