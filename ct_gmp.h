#ifndef __BIPOLY_CT_GMP
#define __BIPOLY_CT_GMP 1

#include <gmpxx.h>

// Data type for the polynomial coefficients.
typedef mpz_class ct;
typedef mpz_class uct;
//typedef ct uct;


// Useful constant values, of type ct.
#define COEFF_0 0
#define COEFF_NEG_1 (-1)
#define COEFF_1 1

// Printf and scanf format characters for a ct
#define CTF "%Zd"
#define CTFS "%Zd"
#define CTFX "%Zx"
#define CTFSX "%Zx"

// How to pass a ct to printf/scanf
#define I(a) (a.get_mpz_t())
#define O(a) (a.get_mpz_t())

// String conversion for ct
//inline ct atoct(const char *nptr) { return atol(nptr); }
inline ct atoct(const char *nptr) {
    mpz_class outp;
    
    gmp_sscanf(nptr, "%Zd", outp.get_mpz_t());
    
    return outp;
}


// Largest positive value representable by a ct.
#define HAVE_COEFF_MAX 0
#define COEFF_MAX (#error !)

#define printf gmp_printf
#define fprintf gmp_fprintf
#define scanf gmp_scanf
#define fscanf gmp_fscanf

// Conversions from ct
#define to_int(a) (a.get_si())
#define to_dbl(a) (a.get_d())
#define to_ht(a) (mpz_to_ll(a.get_mpz_t()))

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

inline void C_mpz_addmul(ct &a, ct &b, ct &c) {    
//    gmp_printf("AM: %Zd+=%Zd*%Zd\n", a.get_mpz_t(), b.get_mpz_t(), c.get_mpz_t());
    mpz_addmul(a.get_mpz_t(), b.get_mpz_t(), c.get_mpz_t());
//    gmp_printf("%Zd\n", a.get_mpz_t());
}
inline void C_mpz_submul(ct &a, ct &b, ct &c) {
//    gmp_printf("SM: %Zd-=%Zd*%Zd\n", a.get_mpz_t(), b.get_mpz_t(), c.get_mpz_t());
    mpz_submul(a.get_mpz_t(), b.get_mpz_t(), c.get_mpz_t());
}

#endif