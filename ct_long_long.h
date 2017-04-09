#ifndef __BIPOLY_CT_LONG_LONG
#define __BIPOLY_CT_LONG_LONG 1


#define HAVE_CT 1
#undef MODULAR

// Data type for the polynomial coefficients.
typedef long long ct;
typedef unsigned long long uct;

// Useful constant values, of type ct.
#define COEFF_0 0
#define COEFF_NEG_1 (-1)
#define COEFF_1 1

// Arithmetic
#define A(a,b) ((a)+(b))
#define S(a,b) ((a)-(b))
#define M(a,b) ((a)*(b))
#define AE(a,b) ((a)+=(b))
#define SE(a,b) ((a)-=(b))
#define ME(a,b) ((a)*=(b))
#define MS(a,b) M(a,b)
#define MSE(a,b) ME(a,b)

// Printf and scanf format characters for a signed and unsigned ct
#define CTF "%qd"
#define CTFS "%qd"
#define CTFX "%qx"
#define CTFSX "%qx"

// How to pass a ct to printf/scanf
#define I(ct) (&ct)
#define O(ct) (ct)

// String conversion for ct
inline ct atoct(const char *nptr) { return strtoll(nptr, (char **) NULL, 10); }

// Largest positive value representable by a ct.
#define COEFF_MAX LONG_LONG_MAX
#define HAVE_COEFF_MAX 1

// Conversions from/to ct
#define from_i(i) (i)
#define to_d(a) ((double) a)
#define to_i(a) ((int) a)
#define to_ht(a) ((ht) a)			// ht is a long long

// Misc
#define ct_sign(a) ((int) a > 0 ? 1 : -1)
#define ct_abs abs
#ifndef __GNUC__
    inline ct abs(ct a) { return (a > 0) ? a : -a; }
#endif

#endif