#ifndef __BIPOLY_CT_NEWGMP
#define __BIPOLY_CT_NEWGMP 1

// Warning: this header allows the use of GMP only in a very unoptimized fashion. It is intended
//  only for use in checking the output of the program when using modular arithmetic.

#include <gmp.h>
#include <gmpxx.h>

#define HAVE_CT 1
#define IS_GMP
#define MAX_DEGREE 400
#define MAX_TERMS 62500
#undef MODULAR

// Data type for the polynomial coefficients.
typedef mpz_class ct;
typedef mpz_class uct;

// Arithmetic
#define A(a,b) ((a)+(b))
#define S(a,b) ((a)-(b))
#define M(a,b) ((a)*(b))
#define AE(a,b) ((a)+=(b))
#define SE(a,b) ((a)-=(b))
#define ME(a,b) ((a)*=(b))
#define MS(a,b) M(a,b)
#define MSE(a,b) ME(a,b)

// Useful constant values, of type ct.
#define COEFF_0 0
#define COEFF_NEG_1 (-1)
#define COEFF_1 1

// Printf and scanf format characters for a signed and unsigned ct
#define CTF "%Zd"
#define CTFS "%Zd"
#define CTFX "%Zx"
#define CTFSX "%Zx"

#define printf gmp_printf
#define fprintf gmp_fprintf
#define scanf gmp_scanf
#define fscanf gmp_fscanf

// How to pass a ct to printf/scanf
#define I(a) (a.get_mpz_t())
#define O(a) (a.get_mpz_t())

// String conversion for ct
inline ct atoct(const char *nptr) {
    mpz_class outp;
    
    gmp_sscanf(nptr, "%Zd", outp.get_mpz_t());
    
    return outp;
}

// Largest positive value representable by a ct.
#undef COEFF_MAX 
#define HAVE_COEFF_MAX 0

// Conversions from/to ct
#define from_i(i) (i)
#define to_i(a) (a.get_si())
#define to_d(a) (a.get_d())
#define to_ht(a) (mpz_to_ll(a.get_mpz_t()))

// BITS_PER_LIMB was renamed!
#ifndef BITS_PER_LIMB
	#define BITS_PER_LIMB GMP_LIMB_BITS
#endif

#if BITS_PER_LIMB >= 64
inline long long mpz_to_ll(const mpz_t z)
{
    return mpz_get_si(z);
}
#else
#if BITS_PER_LIMB == 32
inline long long mpz_to_ll(const mpz_t z)
{
    long long a;
    int sz;
    
    sz = z[0]._mp_size;
    if (sz == 0) return 0;
    a = z[0]._mp_d[0];
    if (!(sz == 1 || sz == -1)) a |= ((long long) (z[0]._mp_d[1])) << 32;
    if (sz < 0) a = -a;
    
    return a;
}
#else
    #error _mp_limb_t not 32 or 64 bits wide
#endif
#endif

// Misc
#define ct_sign(a) (a > 0 ? 1 : -1)
#define ct_abs abs
#ifndef __GNUC__
    inline ct abs(ct a) { return (a > 0) ? a : -a; }
#endif

#endif