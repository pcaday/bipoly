#ifndef __BIPOLY_CT_MSHORT
#define __BIPOLY_CT_MSHORT 1

#define HAVE_CT 1
#define MODULAR 1

#if 1
#define P 65521
#define P_digits 5			// # of decimal digits in P (for atomodp)
#define P_power10 10000			// largest power of 10 less than P (for atomodp)
#else
extern unsigned int P;			// Modulus
#endif

// Data type for the polynomial coefficients.
typedef unsigned int ct;

// Arithmetic
//  These assume int is at least P^2


#define __m_normalize(a) ((unsigned int) ((((unsigned int) a) >= P) ? ((a)-P) : (a)))

// a+b
#define A(a,b) (__m_normalize(a+b))
// a-b
#define S(a,b) (A(a,P-b))
// a*b
#define M(a,b) ((a*b)%P)

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
#define sign(a) (+1)

// Misc
#define ct_abs(a) (a)
#define ct_sign(a) (1)

#endif