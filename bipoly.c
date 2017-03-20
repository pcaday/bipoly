#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <limits.h>
#include <math.h>

#ifdef CT_mshort
	#include "ct_mshort.h"
	#define CT_typestr "C 'int' modulo 2^16-15"
#endif
#ifdef CT_mlong
	#include "ct_mlong.h"
	#define CT_typestr "C 'long' modulo 2^32-5"
#endif
#ifdef CT_long
	#include "ct_long.h"
	#define CT_typestr "C 'long'"
#endif
#ifdef CT_ll
	#include "ct_long_long_osx.h"
	#define CT_typestr "C 'long long'"
#endif
#ifdef CT_gmp
	#include "ct_newgmp.h"
	#define CT_typestr "GMP integer"
#endif

#ifndef HAVE_CT
	#error no coefficient type (ct) defined (use -DCT_xxx)
#elif !HAVE_CT
	#error no coefficient type (ct) defined (use -DCT_xxx)
#endif


#ifndef INFINITY
	#define INFINITY 1e500			// bigger than a float
#endif
#ifdef BIPOLY_IS_SUN				// This is defined as something else in the Sun headers
	#undef RAND_MAX
#endif
#ifndef RAND_MAX
	#define RAND_MAX (0x7FFFFFFFL)		// Maximum value returned by random()
#endif

// VT100 sequences for various text styles
#if 1
	#define VT_RESET "\x1B[0m"
	#define VT_BOLD "\x1B[1m"
	#define VT_RED "\x1B[31m"
	#define VT_BLUE "\x1B[34m"
#else
	#define VT_RESET ""
	#define VT_BOLD ""
	#define VT_RED ""
	#define VT_BLUE ""
#endif


/*
	IMPORTANT: This code does not check for overflow of the coefficients, which are always assumed to fit in the data type ct.
	This may cause problems, even if the input polynomial satisfies these constraints.  You should always verify the
	output independently.
*/

// Rational data type corresponding to ct
typedef struct {
	ct n, d;
} ctq;

// Data type for polynomial hashes
typedef long long ht;

// Printf and scanf format characters for ht
#define HTF "%llx"
#define HTFS "%qx"


#define _min(a, b) ((a) < (b) ? (a) : (b))
#define _max(a, b) ((a) > (b) ? (a) : (b))

#define issign(c)		((c)=='+'||(c)=='-')
#define isnumeric(c)	(isdigit(c)||issign(c))

#define delta_msecs(s,t)		((long) (1000ULL*(t-s)/CLOCKS_PER_SEC))

#ifndef MAX_DEGREE
	#define MAX_DEGREE	400
#endif
#ifndef MAX_TERMS
	#define MAX_TERMS	(MAX_DEGREE+1)*(MAX_DEGREE+1)
#endif
#define MAX_MOVES	1000
#define MAX_SA_MOVES	1000
#define END_MOVES	4
#define MAX_STRING	(1<<22)

#define HIGH_PENALTY	200000		// Penalty used in avoidlists for previous minima.

// Command numbers
#define CMD_ADD1	 1
#define CMD_SUB1	 2
#define CMD_ADD2	 3
#define CMD_SUB2	 4 
#define CMD_FLIP1	 5
#define CMD_FLIP2	 6
#define CMD_SEP1 	 7
#define CMD_SEP2 	 8
#define CMD_TOP	 9			// only comands below CMD_TOP are used during search, the rest are cosmetic
#define CMD_NEG1   (CMD_TOP)
#define CMD_NEG2 (CMD_TOP+1)
#define CMD_SWAP (CMD_TOP+2)		// may also be accomplished via sep1,sep2,sep1 or sep2,sep1,sep2
#define CMD_END   (CMD_TOP+3)		// for the 'loop' function of the user interface

#define CMD_SKEW1P	 109
#define CMD_SKEW2P	 110
#define CMD_SKEW1M	 111
#define CMD_SKEW2M	 112

#ifndef MAX_SEARCH_DEPTH
	#define MAX_SEARCH_DEPTH 80
#endif
int SEARCH_DEPTH=8;

char var1, var2;
unsigned long random_seed=1;
float init_temp=1.;
float bias=0.;
float global_weightf=30.;           // old default was 25
int start_radius=40;
int start_rate=200;                     // old default was 100
int output_sorted=1;
int CD_Function=1;
ct CD_maxc_discourage, CD_maxc_disallow;
int CD_totd_ceiling;
int sait_max_trials=4, sait_no_expand=0;
int sait_runs=0, sa_temp_chg_iters=0;
int g_mark_accepts=0;
int compute_maps=1;

// we intentionally don't treat sep and flip as invertible here, they can take out a factor of x or y
static int inverse_cmds[CMD_END] = { 0, 2, 1, 4, 3, 0, 0, 0, 0, 9, 10, 11 };
//static int inverse_cmds[CMD_END] = { 0, 2, 1, 4, 3, 0, 0, 0, 0, 11, 12, 9, 10, 13, 14, 15 };
// In this table, sep and flip are treated as invertible
static int inverse_cmds_2[CMD_END] = { 0, 2, 1, 4, 3, 5, 6, 7, 8, 9, 10, 11 };
//static int inverse_cmds_2[CMD_END] = { 0, 2, 1, 4, 3, 5, 6, 7, 8, 11, 12, 9, 10, 13, 14, 15 };

static ct bc[MAX_DEGREE][MAX_DEGREE];
static ct coeff[MAX_DEGREE+1][MAX_DEGREE+1];

int sa_move_seq[MAX_SA_MOVES];
int sa_moves;


typedef struct {
	ct c;
	short e[2];
} biterm;

typedef struct {
	int n;
	biterm t[MAX_TERMS];
} bipoly;


// Indices into the vector output by bipoly_evaluate and printed by
//  bipoly_print_sorted2
enum {
	MIND = 0,
	MONIC = 1,			// 0=monic, 1=not monic
	MAXD = 2,
	TOTD = 3,
	TERMS = 4,
	TOTC = 5
};

// Avoidlist data types and declarations
typedef struct {
	ht hash;
	int data;
} al_entry;

typedef struct {
	al_entry *e;
	int size, alloc, exp;
} bpht;

bpht bph_tab = {NULL, 0, 0, -1};
bpht bph_itab = {NULL, 0, 0, -1};

void bipoly_print (bipoly *fs, char v1, char v2, int star, int printInfo, FILE *fp);
void bipoly_print_sorted2 (bipoly *fs, char v1, char v2, int star, int printInfo, FILE *fp);
void bipoly_scan (bipoly *fs, char *str, char v1, char v2);
void bipoly_scan_sorted2 (bipoly *fs, char *str, char v1, char v2);

int bipoly_scan_cmd (char *s);
int bipoly_exec_cmd (bipoly *g, bipoly *f, int cmd);

void bipoly_print_map (int seq[], int n, char v1, char v2, char w1, char w2);
void bipoly_print_map_reduced (int seq[], int n, char v1, char v2, char w1, char w2);
void bipoly_apply_map_to_point (int seq[], int n, const ctq *iv1, const ctq *iv2, ctq *ov1, ctq *ov2, int *undef);
int bipoly_simplify_map(bipoly *init, int move_seq[], int moves);

void bipoly_sm_zero_redundant_ranges(bipoly *init, int move_seq[], int moves);
int bipoly_sm_remove_redundant_ranges_pc(int move_seq[], int moves, bipoly *polys[]);

void bipoly_flip(bipoly *gs, bipoly *fs, int var);
void bipoly_sep(bipoly *gs, bipoly *fs, int var);
void bipoly_trans(bipoly *gs, bipoly *fs, int var, int sign);
void bipoly_skew(bipoly *gs, bipoly *fs, int var, int sign);

void bipoly_xchg(bipoly *fs, bipoly *gs);
void bipoly_add(bipoly *fs, bipoly *gs);
void bipoly_sub(bipoly *fs, bipoly *gs);
void bipoly_add_const(bipoly *fs, ct c);
void bipoly_mult(bipoly *fs, bipoly *gs);
void bipoly_scale(bipoly *fs, ct s);
void bipoly_scale_exp(bipoly *fs, int e0, int e1);
void bipoly_reduce_exps(bipoly *fs);
int bipoly_equal(bipoly *fs, bipoly *gs);
void bipoly_swap (bipoly *gs, bipoly *fs);
int bipoly_negcnt (bipoly *fs);
void bipoly_neg (bipoly *gs, bipoly *fs, int var);
void bipoly_sort (bipoly *fs, int var);
int bipoly_normalize (bipoly *fs, int *move_seq, int move_seq_len);

void bipolyrat_monomial_reduce(bipoly *num, bipoly *den);
void bipolyrat_mult_exp(bipoly *dstn, bipoly *dstd, bipoly *srcn, bipoly *srcd, int power);
void bipolyrat_apply_powers(bipoly *map1n, bipoly *map1d, bipoly *map2n, bipoly *map2d, int (*m1exp)[2], int (*m2exp)[2]);

void bipoly_loop();
void bipoly_auto();
void bipoly_thresh();
void bipoly_sa();
void bipoly_sa_loop();
void bipoly_sa_iter();
void bipoly_sa_maxrad(bipoly *fs, int max_radius, int rate);
void bipoly_samr_demo(bipoly *fs, int max_radius, int rate);
void bipoly_interactive_try_map(int move_seq[], int moves);
int bipoly_apply_map_ui(bipoly *fs, const char *mappath);

// Cost and comparison
void bipoly_evaluate_moniclike(bipoly *fs);
void bipoly_evaluate(bipoly *fs, ct evec[6]);
int bipoly_cmp(bipoly *fs, bipoly *gs, int *same);
int bipoly_cmp_evec(ct fevec[6], ct gevec[6]);
float bipoly_cost(bipoly *fs);
float bipoly_cost_diff(bipoly *fs, bipoly *gs, int *same);
void bipoly_cost_diff_config(bipoly *fs);
int bipoly_count_nt (int moves, int move_seq[]);

// Utilities
float randomfloat();
ct ct_gcd(ct a, ct b);
void bc_init ();
void fprint_evec(FILE *fp, ct evec[6]);
ht bipoly_hash(const bipoly *fs);
ct atoct_gen(const char *str);
#ifdef MODULAR
ct atomodp(const char *str);
#endif

// Avoidlist functions
void bph_add_poly(const bipoly *fs, int penalty);
void bph_add_poly_assocs(bipoly *fs, int penalty);
void bph_add_poly_nbhd3(bipoly *fs);

void bph_begin_insertions(int num_inserts);
void bph_finish_insertions();
void bph_insert(ht hash, int data);
void bph_alloc(bpht *t, int entries);
int bph_search(bpht *t, ht n, int *index, int *data);
int bph_search_first(bpht *t, ht n, int *index, int *data);
void bph_load(const char *name);
void bph_save(const char *name, int append);

int al_entry_compare(const void *a, const void *b);
void bph_sort(bpht *t);
void bph_append(bpht *t, ht hash, int data);


static inline void bipoly_copy(bipoly *fs, bipoly *gs)
	{ register int i, n; fs->n = n = gs->n; for ( i = 0 ; i < n ; i++ ) fs->t[i] = gs->t[i]; }

static inline int bipoly_deg0(bipoly *fs)
	{ register int i, n = fs->n, e, d = 0; for ( i = 0 ; i < n ; i++ ) { e = fs->t[i].e[0]; if (d < e) d = e; } return d; }

static inline int bipoly_deg1(bipoly *fs)
	{ register int i, n = fs->n, e, d = 0; for ( i = 0 ; i < n ; i++ ) { e = fs->t[i].e[1]; if (d < e) d = e; } return d; }

static inline int bipoly_deg(bipoly *fs, int var)
	{ register int i, n = fs->n, e, d = 0; for ( i = 0 ; i < n ; i++ ) { e = fs->t[i].e[var]; if (d < e) d = e; } return d; }

#define bipoly_degs(fs, d0, d1) \
	{ register int i, n = (fs)->n, e; d0 = d1 = 0; for ( i = 0 ; i < n ; i++ ) { e = (fs)->t[i].e[0]; if ((d0) < e) (d0) = e; \
						e = (fs)->t[i].e[1]; if ((d1) < e) (d1) = e; } }

// Find the largest monomial dividing fs and store its exponents in d0, d1
#define bipoly_mon_divider(fs, d0, d1) \
	{ register int i, n = (fs)->n, e; d0 = d1 = INT_MAX; for ( i = 0 ; i < n ; i++ ) { e = (fs)->t[i].e[0]; if ((d0) > e) (d0) = e; \
						e = (fs)->t[i].e[1]; if ((d1) > e) (d1) = e; } }
// Copy an evaluation vector
#define evec_copy(dest, src) \
	{ register int _i; for (_i = 0; _i < 6; _i++) dest[_i]=src[_i]; }

inline void bipoly_set_const(bipoly *fs, ct c)
{
	fs->n = 1;
	fs->t[0].e[0] = 0; fs->t[0].e[1] = 0; fs->t[0].c = c;
}


bipoly f, f0, f1, *g[MAX_SEARCH_DEPTH+1], h, h2;
bipoly bp_tmp[4];

//char bipoly_strbuf[MAX_STRING];
char bipoly_tmpbuf[MAX_STRING];

char *bipoly_get_strfile(const char *name)
{
	FILE *infp;
	char *buf;
	long len;
	
	infp = fopen (name,"r");
	if ( ! infp ) { fprintf (stderr, "bipoly: error opening file %s\n", name); return 0; }
	
	fseek(infp, 0, SEEK_END);
	len = ftell(infp);
	fseek(infp, 0, SEEK_SET);
	buf = (char *) malloc(len+1);
	if ( !buf ) { fprintf (stderr, "bipoly: could not allocate memory to read file\n"); return 0; }
	
	(void) fread (buf, 1, len, infp);
	fclose (infp);
	buf[len]='\0';

	return buf;
}

void usage()
{
	puts ("Usage: bipoly command [-mnp] [-I inputfile] [-O outputfile] [-C cmpfile] [-AI file] [-AO file] [-MI file] [-a runs] [-A func] [-b bias] [-c iters] [-d depth] [-R rate] [-r radius] [-s seed] [-T initial temperature] [-w weight] [poly]");
	puts ("   commands: help auto loop sa saloop samr saiter sarepeat threshold");
	puts ("             apply convmap stats saitertest");
	puts ("             addN subN flipN sepN negN swap skewN[m,p](where N is 1 or 2)");
	puts ("");
	puts ("      help:  print this message");
	puts ("      auto:  simplify input by repeated exhaustive search");
	puts ("      loop:  apply basic operations manually");
	puts ("        sa:  basic simulated annealing");
	puts ("    saloop:  step-by-step basic simulated annealing");
	puts ("      samr:  simulated annealing with maximum radius");
	puts ("    saiter:  repeated simulated annealing with maximum radius (typical command)");
  	puts ("  sarepeat:  try saiter multiple times (# specified by -a), output best result");
	puts (" threshold:  simple threshold search with exponentially decaying threshold");
	puts ("     apply:  apply move sequence to polynomial");
	puts ("   convmap:  convert move sequence to rational map");
	puts ("     stats:  generate performance statistics on 'saiter'");
	puts ("saitertest:  more statistics on 'saiter'");
	puts ("");
	puts (" The following commands are also available from within the 'loop' command:");
	puts ("      addN:  replace varN by varN+1");
	puts ("      subN:  replace varN by varN-1");
	puts ("     flipN:  replace varN by 1/varN");
	puts ("      sepN:  replace varN by 1/varN, multiply other variable by 1/varN");
	puts ("      negN:  replace varN by -varN");
	puts ("      swap:  switch var1 and var2");
	puts ("skew1[p,m]:  replace var1 by (var1+var2) or (var1-var2)");
	puts ("skew2[p,m]:  similarly for var2");
	puts ("");
	puts ("Parameters:");
	puts ("  -I:  input polynomial");
	puts ("  -O:  output file");
	puts ("  -C:  optimized polynomial to compare against for saitertest");
	puts ("         (must use the same variables as the input polynomial)");
	puts (" -AI:  load avoidlist from this file");
	puts (" -AO:  save new avoidlist to this file");
	puts (" -MI:  input map file for apply or convmap");
	puts ("         (series of whitespace-delimited numbers or names like addN)");
	puts ("");
	puts ("  -a:  number of simulated annealing runs (saitertest, stats, sarepeat)");
	puts ("  -A:  simulated annealing cost function (integer from 0 to 7: default 1)");
	puts ("  -b:  decimal amount to adjust energy cost difference (sa command)");
	puts ("  -c:  number of iterations between temperature changes for simulated annealing");
	puts ("  -d:  search depth (auto command; default: 8)");
	puts ("  -R:  cooling rate constant for simulated annealing (default: 200)");
	puts ("         (the higher the constant, the slower the cooling)");
	puts ("  -r:  initial radius for simulated annealing (default: 40)");
	puts ("  -s:  random seed (default: clock time)");
	puts ("  -T:  initial temperature for simulated annealing (default: 1.0)");
	puts ("  -w:  global weight factor for cost functions 5 - 7 (default: 30.0)");
	puts ("");
	puts ("  -m:  do not compute transformations (saves time)");
	puts ("  -n:  do not expand SA radius (saiter command)");
	puts ("  -p:  do not print output polynomial sorted on second variable");
	puts ("");
	puts ("Compiled-in coefficient type: " CT_typestr);
	exit(2);
}

int main (int argc, char *argv[])
{
	FILE *outfp = NULL;
	char *s, tmpc;
	int i;
	char *in_poly_str = NULL, *comp_poly_str = NULL;
	const char *bph_out = NULL, *map_in = NULL;
	const char *in_path = NULL;
	bipoly *gstore;

	random_seed = (long) time(NULL);
	
	if ( argc < 3 || strstr(argv[1], "help") ) usage();

	gstore = (bipoly *) malloc(sizeof(bipoly) * (MAX_DEGREE + 1));
	if (gstore == NULL) { fprintf(stderr, "bipoly: could not allocate memory\n"); return 3; }
	for (i = 0; i < MAX_DEGREE+1; i++) { g[i] = &gstore[i]; }

	// parse args
	for (i = 2; i < argc; i++) {
		if (i < (argc - 1)) {
			if (!strcmp(argv[i], "-c"))	{ sa_temp_chg_iters = atoi(argv[++i]); continue; }
			if (!strcmp(argv[i], "-a"))	{ sait_runs = atoi(argv[++i]); continue; }
			if (!strcmp(argv[i], "-A"))	{ CD_Function = atoi(argv[++i]); continue; }
			if (!strcmp(argv[i], "-b"))	{ bias = atof(argv[++i]); continue; }	
			if (!strcmp(argv[i], "-w"))	{ global_weightf = atof(argv[++i]); continue; }					if (!strcmp(argv[i], "-R"))	{ start_rate = atoi(argv[++i]); if (start_rate <= 0) { fprintf(stderr, "bipoly: invalid starting rate\n"); return 0; }
							continue; }
			if (!strcmp(argv[i], "-r"))	{ start_radius = atoi(argv[++i]); if (start_radius <= 0) { fprintf(stderr, "bipoly: invalid starting radius\n"); return 0; }
							continue; }
			if (!strcmp(argv[i], "-T"))	{ init_temp = atof(argv[++i]); continue; }
			if (!strcmp(argv[i], "-s"))	{ random_seed = atol(argv[++i]); if (random_seed == 0) { fprintf(stderr, "bipoly: invalid seed\n"); return 0; }
							continue; }
			if (!strcmp(argv[i], "-d"))	{ SEARCH_DEPTH = atoi(argv[++i]); if (SEARCH_DEPTH <= 0) { fprintf(stderr, "bipoly: invalid search depth\n"); return 0; }
							continue; }
			if (!strcmp(argv[i], "-AO"))	{ bph_out = argv[++i]; continue; }
			if (!strcmp(argv[i], "-AI"))	{ bph_load(argv[++i]); continue; }
			if (!strcmp(argv[i], "-MI"))	{ map_in = argv[++i]; continue; }
			if (!strcmp(argv[i], "-O"))	{ 
				if (outfp) { fprintf(stderr, "bipoly: more than one output polynomial file was given\n"); return 0; }
				i++;
				outfp = fopen (argv[i], "w");
				if (!outfp) { fprintf (stderr, "bipoly: error opening file %s for writing\n", argv[i]); return 0; }
				continue;
			}
			if (!strcmp(argv[i], "-o"))	{ 
				if (outfp) { fprintf(stderr, "bipoly: more than one output polynomial file was given\n"); return 0; }
				i++;
				outfp = fopen (argv[i], "a");
				if (!outfp) { fprintf (stderr, "bipoly: error opening file %s for writing\n", argv[i]); return 0; }
				continue;
			}
			if (!strcmp(argv[i], "-I"))	{
				if (in_poly_str) { fprintf(stderr, "bipoly: more than one input polynomial was given\n"); return 0; }
				i++;  in_poly_str = bipoly_get_strfile(argv[i]);
				in_path = argv[i];
				continue;
			}
			if (!strcmp(argv[i], "-C"))	{
				i++;  comp_poly_str = bipoly_get_strfile(argv[i]);
				continue;
			}
		}
		if (!strcmp(argv[i], "-p"))	{ output_sorted = 0; continue; }
		if (!strcmp(argv[i], "-n"))	{ sait_no_expand = 1; continue; }
		if (!strcmp(argv[i], "-m"))	{ compute_maps = 0; continue; }
		if (in_poly_str) usage();
		in_poly_str = argv[i];
	}

	if ( SEARCH_DEPTH > MAX_SEARCH_DEPTH ) { fprintf (stderr, "bipoly: search depth cannot be greater than MAX_SEARCH_DEPTH=%d\n", MAX_SEARCH_DEPTH); return 0; }
	if ( sa_temp_chg_iters == 0 ) {
		sa_temp_chg_iters = start_radius / 10;
		sa_temp_chg_iters = _min(sa_temp_chg_iters, 25);
	}

	printf("Search depth       %d\n", SEARCH_DEPTH);
	printf("Random seed value  %d\n", (int) random_seed);
	srandom(random_seed);

	if (!in_poly_str) { fprintf(stderr, "bipoly: no input polynomial\n"); return 0; }
	for ( s = in_poly_str ; *s && ! isalpha(*s) ; s++ );
	if ( ! *s ) fprintf (stderr, "Can't find two variables in specified poly\n");
	var1 = *s;
	for ( ; *s && (*s == var1 || ! isalpha(*s))  ; s++ );
	if ( ! *s ) fprintf (stderr, "Can't find two variables in specified poly\n");
	var2 = *s;
	for ( ; *s && (*s == var1 || *s == var2 || ! isalpha(*s))  ; s++ );
	if ( *s ) fprintf (stderr, "Please use single alphabetic characters for each the 2 variables\n");
	if (tolower(var1) > tolower(var2)) { tmpc = var1; var1 = var2; var2 = tmpc; }
	
	bipoly_scan_sorted2 (&f,in_poly_str,var1,var2);

	printf ("Input p(%c,%c) = ", var1, var2); bipoly_print (&f,var1,var2,0,1,stdout);
	
	bc_init();

	if ( strcmp (argv[1],"loop") == 0 ) {
		bipoly_loop();
		if ( outfp != NULL ) bipoly_print(&f,var1,var2,1,0,outfp);
	} else if ( strcmp (argv[1], "auto") == 0 ) {
		bipoly_auto();
		if ( outfp != NULL ) bipoly_print(&h,var1,var2,1,0,outfp);
	} else if ( strcmp (argv[1], "eval") == 0 ) {
		ct evec[6];
		bipoly_evaluate(&f, evec);
		fprint_evec(outfp ? outfp : stdout, evec);
	} else if ( strcmp(argv[1], "sa") == 0) {
		bipoly_sa();
		if ( outfp != NULL ) bipoly_print(&f,var1,var2,1,0,outfp);
	} else if ( strcmp(argv[1], "samr") == 0) {
		sa_moves = 0;
		bipoly_sa_maxrad(&f, start_radius, start_rate);
		if ( outfp != NULL ) bipoly_print(&f,var1,var2,1,0,outfp);
	} else if ( strcmp(argv[1], "samr-demo") == 0) {
		sa_moves = 0;
		bipoly_samr_demo(&f, start_radius, start_rate);
	} else if ( strcmp(argv[1], "saiter") == 0) {
		bipoly_sa_iter();
		if ( bph_out != NULL ) bph_add_poly_nbhd3(&h);
		if ( outfp != NULL ) bipoly_print(&f,var1,var2,1,0,outfp);
	} else if ( strcmp(argv[1], "saloop") == 0) {
		bipoly_sa_loop();
		if ( outfp != NULL ) bipoly_print(&f,var1,var2,1,0,outfp);
	} else if ( strcmp(argv[1], "sarepeat") == 0) {
		int i, j, cmp;
		int best_run = 0, best_moves = 0, best_move_seq[MAX_MOVES];
		ct bevec[6];
		clock_t start, end;
		long t, tmin = LONG_MAX, tmax = 0, ttot = 0;
		
		if (sait_runs == 0) sait_runs = 6;

		bipoly_evaluate(&f, bevec);
		for (i = 0; i < sait_runs; i++) {
			ct evec[6];

			start = clock();
			random_seed += 135 + (long) clock(); 
			srandom(random_seed);
			
			bipoly_sa_iter();
			bph_add_poly_nbhd3(&h);
			bipoly_evaluate(&f, evec);
			fprint_evec(stderr, evec); fprintf(stderr, ",  %d moves\n", sa_moves);
			cmp = bipoly_cmp_evec(bevec, evec);
			if (cmp > 0 || ((cmp == 0) && (sa_moves < best_moves))) {
				for (j = 0; j < 6; j++) bevec[j] = evec[j];
				best_moves = sa_moves;
				for (j = 0; j < sa_moves; j++) best_move_seq[j] = sa_move_seq[j];

			if (cmp > 0) best_run = i + 1;                              // record iteration where we first found the optimal polynomial
			}
			bipoly_copy(&f, &f0);

			end = clock();  t = delta_msecs(start, end);
			tmin = _min(tmin, t);  tmax = _max(tmax, t);  ttot += t;
		}
		
		for (j = 0; j < best_moves; j++) bipoly_exec_cmd(&f, &f, best_move_seq[j]);
		
		if (best_run > 0)
			printf("\n\nBest polynomial found at iteration %d", best_run);
		else
			printf("\n\nInput polynomial not improved after %d tries", sait_runs);

		printf ("\n\nMove sequence: "); for ( i = 0 ; i < best_moves ; i++ ) printf ("%d ", best_move_seq[i]);
		printf("\n\nMap:\n");
		bipoly_print_map_reduced(best_move_seq, best_moves, var1, var2, 'U', 'V');
		printf("\nFinal polynomial:\n");
		
		bipoly_sort(&f, 1);
		if (output_sorted)
			bipoly_print_sorted2(&f, var1, var2, 0, 1, stdout);
		else
			bipoly_print(&f, var1, var2, 0, 1, stdout);
		
		printf("Path length %d, reduction time (min/avg/max): %ld/%ld/%ld msecs\n", best_moves, tmin, ttot / sait_runs, tmax);
		
		if ( outfp != NULL ) bipoly_print(&f,var1,var2,1,0,outfp);
	} else if ( strcmp(argv[1], "saitertest") == 0 ) {
		#define MAX_SAIT_RUNS 400
		double avg[6];
		ct best[6], all[MAX_SAIT_RUNS][6], cmp[6];
		int j, aoa, better, oruns, ac[6], moves[MAX_SAIT_RUNS], cmoves, tmoves;
		clock_t start, end;
		
		if (sait_runs == 0) sait_runs = 20;
		if (sait_runs > MAX_SAIT_RUNS) { fprintf(stderr, "bipoly: requested more than MAX_SAIT_RUNS=%d", MAX_SAIT_RUNS); exit(0); }
		
		if (comp_poly_str != NULL) {
			fprintf(stderr, "Comparing against:\n");
			bipoly_scan(&f1,comp_poly_str,var1,var2);
			bipoly_print(&f1,var1,var2,0,0,stderr);
			bipoly_evaluate(&f1,cmp);
			putc('\n', stderr);
			fprint_evec(stderr, cmp);
			putc('\n', stderr);
		} else {
			fprintf(stderr, "No comparison polynomial specified, using original poly...\n");
			bipoly_evaluate(&f, cmp);
		}
		
		aoa = 0;  cmoves = 0;  tmoves = 0;  better = 0;  oruns = 0;
		for (j = 0; j < 6; j++) { best[j] = cmp[j];  ac[j] = 0;  avg[j] = 0.; }
		
		start = clock();
		for (i = 0; i < sait_runs; i++) {
			random_seed += 135 + (long) clock();
			srandom(random_seed);
			
			fprintf(stderr, "\r%4d/%d", i+1, sait_runs);
			bipoly_sa_iter();
			moves[i] = sa_moves;
			bipoly_evaluate(&f, all[i]);
			bipoly_copy(&f, &f0);
			
//			for (j = 0; j < 6; j++) { if (all[i][j] < best[j]) best[j] = all[i][j]; }
			if (bipoly_cmp_evec(best, all[i]) > 0) for (j=0; j < 6; j++) best[j] = all[i][j];
		}
		end = clock();
		fprintf(stderr, "\r");

		for (i = 0; i < sait_runs; i++) {
			// Count number of results at least as good as the comparison poly
			for (j = 0; j < 6; j++) {
				if (all[i][j] < cmp[j]) { better++; break; }
				else if (all[i][j] > cmp[j]) { break; }
			}
			if ( j==6 ) {
				aoa++;
			}
			tmoves += moves[i];
			if (bipoly_cmp_evec(all[i], best) == 0) {
				oruns++; cmoves += moves[i];
			}

			// Get the pertinent averages for each piece of information
			for (j = 0; j < 6; j++) {
				avg[j] += to_d(all[i][j]);  ac[j]++;
				if (all[i][j] > best[j]) break;
			}
			
//			fprintf(outfp ? outfp : stdout, "(%ld,%ld,%ld,%ld,%ld,%ld) - %d\n", all[i][0], all[i][1], all[i][2], all[i][3], all[i][4], all[i][5], moves[i]);
		}

		fprintf(outfp ? outfp : stdout, "%2.1f%%, %2.1f%%, (%.1f,%.1f,%.1f,%.1f,%.1f,%.1f), %.1f/%.1f, %.0f ms   ", (double) (100 * better) / sait_runs, (double) (100 * aoa) / sait_runs, avg[0] / ac[0], avg[1] / ac[1], avg[2] / ac[2], avg[3] / ac[3], avg[4] / ac[4], avg[5] / ac[5], (double) tmoves / sait_runs, (double) cmoves / oruns, (double) (delta_msecs(start,end)) / sait_runs);
		fprint_evec(outfp ? outfp : stdout, best);
		fprint_evec(outfp ? outfp : stdout, cmp);
		putc('\n', outfp ? outfp : stdout);
	} else if ( strcmp(argv[1], "threshold") == 0) {
		bipoly_thresh();
		if ( outfp != NULL ) bipoly_print(&f,var1,var2,1,0,outfp);
	} else if ( strcmp(argv[1], "apply") == 0) {
		if ( map_in == NULL ) { fprintf(stderr, "bipoly: no map was specified (use the -MI flag)\n"); exit(0); }
		(void) bipoly_apply_map_ui(&f, map_in);
		bipoly_sort(&f,1);
		if ( outfp != NULL ) {
			if (output_sorted)
				bipoly_print_sorted2(&f,var1,var2,1,1,outfp);
			else
				bipoly_print(&f,var1,var2,1,1,outfp);
		}
#ifdef IS_GMP
	} else if ( strcmp(argv[1], "papertab") == 0) {
		ct oevec[6], nevec[6], a;
		int log_omaxc, log_nmaxc, q;
		int i, moves;
		
		if ( map_in == NULL ) { fprintf(stderr, "bipoly: no map was specified (use the -MI flag)\n"); exit(0); }
		
		bipoly_evaluate(&f, oevec);
		log_omaxc = 0;
		for (i = 0; i < f.n; i++) {
			a = abs(f.t[i].c) - 1;
			q = mpz_sizeinbase(a.get_mpz_t(), 2);		// get number of bits in coefficient
			if (q > log_omaxc) log_omaxc = q;
		}
		
		moves = bipoly_apply_map_ui(&f, map_in);
		
		if ( outfp != NULL ) {
			bipoly_evaluate(&f, nevec);
			log_nmaxc = 0;
			for (i = 0; i < f.n; i++) {
				q = mpz_sizeinbase(f.t[i].c.get_mpz_t(), 2);
				if (q > log_nmaxc) log_nmaxc = q;
			}
			fprintf(outfp, "\n&& %d & %d & %d & %d & %d & %d & %d\n", to_i(oevec[MIND]), to_i(oevec[TERMS]), log_omaxc, moves, to_i(nevec[MIND]), to_i(nevec[TERMS]), log_nmaxc);
		}
#endif
	} else if ( strcmp(argv[1], "convmap") == 0) {
		char *nmp;
		static const char *sep = " \t\n\r";

		if ( map_in == NULL ) { fprintf(stderr, "bipoly: no map was specified (use the -MI flag)\n"); exit(0); }
		
		sa_moves = 0;
		nmp = strtok(bipoly_get_strfile(map_in), sep);
		while (nmp) {
			if (sa_moves >= MAX_MOVES) { fprintf(stderr, "bipoly: move sequence too long, increase MAX_MOVES\n"); exit(0); }
			sa_move_seq[sa_moves++] = bipoly_scan_cmd(nmp);
			if (sa_move_seq[sa_moves - 1] == 0)
				break;
			nmp = strtok(NULL, sep);
		}

		bipoly_print_map_reduced(sa_move_seq, sa_moves, var1, var2, 'U', 'V');
	} else if ( strcmp(argv[1], "stats") == 0) {
		int i, j;
		clock_t start, end;
		FILE *realout = outfp ? outfp : stdout;
		
		long max_ms = 0, min_ms = LONG_MAX, total_ms = 0;
		ct best_evec[6], worst_evec[6], total_evec[6];

		for (i = 0; i < 6; i++) {
			#if HAVE_COEFF_MAX
				best_evec[i] = COEFF_MAX;
			#else
				best_evec[i] = LONG_MAX;
			#endif
			worst_evec[i] = 0;
			total_evec[i] = 0;
		}
		if (sait_runs == 0) sait_runs = 6;

		for (i = 0; i < sait_runs; i++) {
			ct evec[6];
			
			random_seed += 135 + (long) clock();
			srandom(random_seed);
			
			start = clock();
			bipoly_sa_iter();
			end = clock();
			
			end = delta_msecs(start, end);
			if (end > max_ms) max_ms = end;
			if (end < min_ms) min_ms = end;
			total_ms += end;
			
//			bph_add_poly_nbhd3(&h);				/* don't want this for stats */
			bipoly_evaluate(&f, evec);
			fprint_evec(stderr, evec); fprintf(stderr, ",  %d moves\n", sa_moves);
			if (bipoly_cmp_evec(best_evec, evec) > 0) evec_copy(best_evec, evec);
			if (bipoly_cmp_evec(worst_evec, evec) < 0) evec_copy(worst_evec, evec);
			for (j = 0; j < 6; j++) total_evec[j] += evec[j];
			
			bipoly_copy(&f, &f0);		
		}
		
		fprintf(realout, "%s\n", in_path ? in_path : in_poly_str);
		fprintf(realout, "worst/avg/best time (ms): %ld/%ld/%ld\n", max_ms, total_ms / sait_runs, min_ms);
		fprintf(realout, "worst/avg/best evec:\n  (" CTF "," CTF "," CTF "," CTF "," CTF "," CTF ")", O(worst_evec[0]), O(worst_evec[1]), O(worst_evec[2]), O(worst_evec[3]), O(worst_evec[4]), O(worst_evec[5]));
		fprintf(realout, "\n  (%2f,%2f,%2f,%2f,%2f,%1f)", to_d(total_evec[0]) / sait_runs, to_d(total_evec[1]) / sait_runs, to_d(total_evec[2]) / sait_runs, to_d(total_evec[3]) / sait_runs, to_d(total_evec[4]) / sait_runs, to_d(total_evec[5]) / sait_runs);
		fprintf(realout, "\n  (" CTF "," CTF "," CTF "," CTF "," CTF "," CTF ")\n", O(best_evec[0]), O(best_evec[1]), O(best_evec[2]), O(best_evec[3]), O(best_evec[4]), O(best_evec[5]));		
	} else {
		bipoly_exec_cmd (&f,&f, bipoly_scan_cmd(argv[1]));
		puts ("Result: ");
		bipoly_print(&f,var1,var2,1,1,stdout);
		if ( outfp != NULL ) bipoly_print(&f,var1,var2,1,0,outfp);
	}
		
	if ( outfp != NULL ) fclose(outfp);
	if ( bph_out != NULL ) bph_save(bph_out, 0);
	
	return 0;
}

void bipoly_loop()
{
	int move_seq[MAX_MOVES], move_seq_len;
	int i;
	char buf[256];
	
	move_seq_len = 0;
	bipoly_copy(&f0,&f);
	for(;;) {
		printf ("cmd> ");  (void) fgets(buf,sizeof(buf),stdin); buf[strlen(buf)-1]='\0';
		i = bipoly_scan_cmd(buf);
		switch (bipoly_exec_cmd(&f,&f,i)) {
			case -1: continue;		// got an unrecognized command
			case 0: return;			// user typed 'end'
		}
		move_seq[move_seq_len++]=i;
		puts ("Result: ");
		bipoly_print(&f,'x','y',1,1,stdout);
		bipoly_print_map_reduced(move_seq, move_seq_len,var1,var2,'x','y');
	}
}

void bipoly_auto()
{
	clock_t start, restart, end = 0;
	int move_cmd[MAX_SEARCH_DEPTH], move_len;
	int move_seq[MAX_MOVES], move_seq_len, path_len;
	int i, j[MAX_SEARCH_DEPTH], m = 0, depth, max_depth, same;
	char w1, w2;
	
	start = restart = clock();

	bipoly_copy(&f0,&f);
// save swapped and sorted copy of initial curve for comparison later
//		bipoly_copy(&f0,&f);
//		bipoly_swap(&f0);
//		bipoly_sort(&f0,1);
	move_seq_len = max_depth = same = 0;
	bipoly_copy(&h,&f);
	for ( depth = 1 ;; depth++ ) {
		move_len = 0;
		bipoly_copy(g[0],&f);
		// hard wire search depth for speed
		for ( j[0] = 1 ; j[0] < CMD_TOP ; j[0]++ ) {
			bipoly_exec_cmd(g[1],g[0],j[0]);
			if ( !g[1]->n ) continue;
			if ( depth == 1 ) {
				if ( bipoly_cmp (g[1],&h,&same) < 0 ) { bipoly_copy (&h,g[1]); move_len = 1; for ( i=0 ; i < move_len ; i++ ) move_cmd[i]=j[i]; printf("Move (%d)\n", j[0]); bipoly_print(&h,var1,var2,0,1,stdout); } 
				continue;
			}
			for ( j[1] = 1 ; j[1] < CMD_TOP ; j[1]++ ) {
				if ( j[1] == inverse_cmds[j[0]] ) continue;
				bipoly_exec_cmd(g[2],g[1],j[1]);
				if ( !g[2]->n ) continue;
				if ( depth == 2 ) {
					if ( bipoly_cmp (g[2],&h,&same) < 0 ) { bipoly_copy (&h,g[2]); move_len = 2; for ( i=0 ; i < move_len ; i++ ) move_cmd[i]=j[i]; printf("Move (%d,%d)\n", j[0],j[1]);  bipoly_print(&h,var1,var2,0,1,stdout); }
					continue;
				}
				for ( j[2] = 1 ; j[2] < CMD_TOP ;  j[2]++ ) {
					if ( j[2] == inverse_cmds[j[1]] ) continue;
					bipoly_exec_cmd(g[3],g[2],j[2]);
					if ( !g[3]->n ) continue;
					if ( depth == 3 ) {
						if ( bipoly_cmp (g[3],&h,&same) < 0 ) { bipoly_copy (&h,g[3]); move_len = 3; for ( i=0 ; i < move_len ; i++ ) move_cmd[i]=j[i];printf("Move (%d,%d,%d)\n", j[0],j[1],j[2]);  bipoly_print(&h,var1,var2,0,1,stdout); }
						continue;
					}
					for ( j[3] = 1 ; j[3] < CMD_TOP ;  j[3]++ ) {
						if ( j[3] == inverse_cmds[j[2]] ) continue;
						bipoly_exec_cmd(g[4],g[3],j[3]);
						if ( !g[4]->n ) continue;
						if ( depth == 4 ) {
							if ( bipoly_cmp (g[4],&h,&same) < 0 ) { bipoly_copy (&h,g[4]); move_len = 4; for ( i=0 ; i < move_len ; i++ ) move_cmd[i]=j[i]; printf("Move (%d,%d,%d,%d)\n", j[0],j[1],j[2],j[3]); bipoly_print(&h,var1,var2,0,1,stdout); }
							continue;
						}
						for ( j[4] = 1 ; j[4] < CMD_TOP ; j[4]++ ) {
							if ( j[4] == inverse_cmds[j[3]] ) continue;
							bipoly_exec_cmd(g[5],g[4],j[4]);
							if ( !g[5]->n ) continue;
							if ( depth == 5 ) {
								if ( bipoly_cmp (g[5],&h,&same) < 0 ) { bipoly_copy (&h,g[5]); move_len = 5; for ( i=0 ; i < move_len ; i++ ) move_cmd[i]=j[i]; printf("Move (%d,%d,%d,%d,%d)\n", j[0],j[1],j[2],j[3],j[4]); bipoly_print(&h,var1,var2,0,1,stdout); }
								continue;
							}
							for ( j[5] = 1 ; j[5] < CMD_TOP ;  j[5]++ ) {
								if ( j[5] == inverse_cmds[j[4]] ) continue;
								bipoly_exec_cmd(g[6],g[5],j[5]);
								if ( !g[6]->n ) continue;
								if ( depth == 6 ) {
									if ( bipoly_cmp (g[6],&h,&same) < 0 ) { bipoly_copy (&h,g[6]);move_len = 6; for ( i=0 ; i < move_len ; i++ ) move_cmd[i]=j[i]; printf("Move (%d,%d,%d,%d,%d,%d)\n", j[0],j[1],j[2],j[3],j[4],j[5]); bipoly_print(&h,var1,var2,0,1,stdout); }
									continue;
								}
								for ( j[6] = 1 ; j[6] < CMD_TOP ;  j[6]++ ) {
									if ( j[6] == inverse_cmds[j[5]] ) continue;
									bipoly_exec_cmd(g[7],g[6],j[6]);
									if ( !g[7]->n ) continue;
									if ( depth == 7 ) {
										if ( bipoly_cmp (g[7],&h,&same) < 0 ) { bipoly_copy (&h,g[7]); move_len = 7; for ( i=0 ; i < move_len ; i++ ) move_cmd[i]=j[i]; printf("Move (%d,%d,%d,%d,%d,%d,%d)\n", j[0],j[1],j[2],j[3],j[4],j[5], j[6]); bipoly_print(&h,var1,var2,0,1,stdout); }
										continue;
									}		
									for ( j[7] = 1 ; j[7] < CMD_TOP ; j[7]++ ) {
										if ( j[7] == inverse_cmds[j[6]] ) continue;
										bipoly_exec_cmd(g[8],g[7],j[7]);
										if ( !g[8]->n ) continue;
										if ( depth == 8 ) {
											if ( bipoly_cmp (g[8],&h,&same) < 0 ) { bipoly_copy (&h,g[8]); move_len = 8; for ( i=0 ; i < move_len ; i++ ) move_cmd[i]=j[i];  printf("Move (%d,%d,%d,%d,%d,%d,%d,%d)\n", j[0],j[1],j[2],j[3],j[4],j[5],j[6],j[7]); bipoly_print(&h,var1,var2,0,1,stdout); }
											continue;
										}
										for ( j[8] = 1 ; j[8] < CMD_TOP ; j[8]++ ) {
											if ( j[8] == inverse_cmds[j[7]] ) continue;
											bipoly_exec_cmd(g[9],g[8],j[8]);
											if ( !g[9]->n ) continue;
											if ( depth == 9 ) {
												if ( bipoly_cmp (g[9],&h,&same) < 0 ) { bipoly_copy (&h,g[9]); move_len = 9; for ( i=0 ; i < move_len ; i++ ) move_cmd[i]=j[i];  printf("Move (%d,%d,%d,%d,%d,%d,%d,%d,%d)\n", j[0],j[1],j[2],j[3],j[4],j[5],j[6],j[7],j[8]); bipoly_print(&h,var1,var2,0,1,stdout); }
												continue;
											}
											for ( j[9] = 1 ; j[9] < CMD_TOP ; j[9]++ ) {
												if ( j[9] == inverse_cmds[j[7]] ) continue;
												bipoly_exec_cmd(g[10],g[9],j[9]);
												if ( !g[10]->n ) continue;
												if ( bipoly_cmp (g[10],&h,&same) < 0 ) { bipoly_copy (&h,g[10]); move_len = 10; for ( i=0 ; i < move_len ; i++ ) move_cmd[i]=j[i];  printf("Move (%d,%d,%d,%d,%d,%d,%d,%d,%d,%d)\n", j[0],j[1],j[2],j[3],j[4],j[5],j[6],j[7],j[8],j[9]); bipoly_print(&h,var1,var2,0,1,stdout); }
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		if ( ! move_len ) {
			end = clock();
			printf("depth=%d, time=%ld msecs\n", depth,delta_msecs(restart,end));
			restart = clock();
			if ( depth >= SEARCH_DEPTH ) break;
			continue;
		}
		if ( move_seq_len + move_len + END_MOVES > MAX_MOVES ) { fprintf (stderr, "Exceeded MAX_MOVES=%d!\n", MAX_MOVES); return; }
		printf ("*** Move %d ***  ", ++m); for ( i = 0 ; i < move_len ; i++ ) { printf ("%d ",move_cmd[i]);  move_seq[move_seq_len++] = move_cmd[i]; } 
		puts(""); 
		bipoly_copy (&f, &h);
		bipoly_sort (&f, 0);
		bipoly_print(&f,'x','y',0,1,stdout);
		bipoly_print_map_reduced(move_seq, move_seq_len, var1, var2, 'U', 'V');
		if ( depth > max_depth ) max_depth = depth;
		depth = 0; restart = clock();
	}
	printf ("Complete move sequence: "); for ( i = 0 ; i < move_seq_len ; i++ ) printf ("%d ", move_seq[i]); puts("");
	path_len = move_seq_len;

	move_seq_len = bipoly_normalize(&h, move_seq, move_seq_len);
	w1 = 'U'; w2 = 'V';
//		puts ("\nSwapped Initial Curve (for comparison only):");
//		bipoly_print_sorted2(&f0,'x','y',0,1);
	puts("\nFinal Result:");
	if (output_sorted)
		bipoly_print_sorted2(&h,'x','y',1,1,stdout);
	else
		bipoly_print(&h,'x','y',1,1,stdout);
	bipoly_print_map_reduced(move_seq, move_seq_len, var1, var2, w1, w2);
	if (output_sorted)
		bipoly_print_sorted2(&h,'x','y',0,1,stdout);
	else
		bipoly_print(&h,'x','y',0,1,stdout);
	printf("Path length %d, Max Successful Depth: %d  Total time: %ld msecs\n", path_len, max_depth, delta_msecs(start,end));
#if 0
	bipoly_copy(&f, &f0);
	for (i = 0; i < path_len; i++) {
		printf("%02d:  ", i);
		bipoly_exec_cmd(&f, move_seq[i]);
		bipoly_evaluate_moniclike(&f);
	}
#endif
}

int bipoly_apply_map_ui(bipoly *fs, const char *mappath)
{
	char *nmp;
	static const char *sep = " \t\n\r";
	int cmd;
	int moves = 0;
	
	nmp = strtok(bipoly_get_strfile(mappath), sep);
	while (nmp) {
		cmd = bipoly_scan_cmd(nmp);
		bipoly_exec_cmd(fs, fs, cmd);
		nmp = strtok(NULL, sep);
		moves++;
	}
	printf("Result:\n");
	bipoly_print(&f,var1,var2,1,1,stdout);
	
	return moves;
}

int bipoly_scan_cmd (char *cmd)
{
	int cn;
	
	cn = atoi(cmd);
	if ( cn != 0 ) return cn;
	
	if ( strcmp (cmd,"flip1") == 0 ) return CMD_FLIP1;
	else if ( strcmp (cmd,"flip2") == 0 ) return CMD_FLIP2;
	else if ( strcmp (cmd, "sep1") == 0 ) return CMD_SEP1;
	else if ( strcmp (cmd, "sep2") == 0 ) return CMD_SEP2;
	else if ( strcmp (cmd, "add1") == 0 ) return CMD_ADD1;
	else if ( strcmp (cmd, "add2") == 0 ) return CMD_ADD2;
	else if ( strcmp (cmd, "sub1") == 0 ) return CMD_SUB1;
	else if ( strcmp (cmd, "sub2") == 0 ) return CMD_SUB2;
	else if ( strcmp (cmd, "neg1") == 0 ) return CMD_NEG1;
	else if ( strcmp (cmd, "neg2") == 0 ) return CMD_NEG2;
	else if ( strcmp (cmd, "swap") == 0 ) return CMD_SWAP;
	else if ( strcmp (cmd, "end") == 0 ) return CMD_END;
	else if ( strcmp (cmd, "skew1p") == 0 ) return CMD_SKEW1P;
	else if ( strcmp (cmd, "skew2p") == 0 ) return CMD_SKEW2P;
	else if ( strcmp (cmd, "skew1m") == 0 ) return CMD_SKEW1M;
	else if ( strcmp (cmd, "skew2m") == 0 ) return CMD_SKEW2M;
	else return 0;
}


int bipoly_exec_cmd (bipoly *g, bipoly *f, int cmd)
{
	int result = 1;
	
	switch (cmd) {
	case CMD_FLIP1: bipoly_flip(g,f,0); break;
	case CMD_FLIP2: bipoly_flip(g,f,1); break;
	case CMD_SEP1: bipoly_sep(g,f,0); break;
	case CMD_SEP2: bipoly_sep(g,f,1); break;
	case CMD_ADD1: bipoly_trans(g,f,0,1); break;
	case CMD_ADD2: bipoly_trans(g,f,1,1); break;
	case CMD_SUB1: bipoly_trans(g,f,0,-1); break;
	case CMD_SUB2: bipoly_trans(g,f,1,-1); break;
	case CMD_NEG1: bipoly_neg(g,f,0); break;
	case CMD_NEG2: bipoly_neg(g,f,1); break;
	case CMD_SWAP: bipoly_swap(g,f); break;
	case CMD_END: result = 0; break;
	case CMD_SKEW1P: bipoly_skew(g,f,0,1); break;
	case CMD_SKEW2P: bipoly_skew(g,f,1,1); break;
	case CMD_SKEW1M: bipoly_skew(g,f,0,-1); break;
	case CMD_SKEW2M: bipoly_skew(g,f,1,-1); break;
	case 100: printf(HTF "\n", bipoly_hash(f)); break;
	case 101: bph_add_poly(f, 10); break;
	case 102: bph_begin_insertions(16); bph_add_poly_assocs(f, 10); bph_finish_insertions(); break;
	default: puts ("Unknown command"); result = -1; break;
	}
	
	return result;
}

void bipoly_evaluate(bipoly *fs, ct evec[6])
{
	int fmind, fmaxd, ftotd, fmonic, minv;
#ifndef MODULAR
	ct ftotc;
#endif
	register int i, j, d1, d2;
	biterm *f = &fs->t[0];
	int m = fs->n;
	
	d1 = d2 = ftotd = 0;
#ifndef MODULAR
	ftotc = from_i(0);
#endif
	for ( i = 0 ; i < m ; i++ ) {
		if ( f[i].e[0] > d1 ) d1 = f[i].e[0];
		if ( f[i].e[1] > d2 ) d2 = f[i].e[1];
		j = f[i].e[0]+f[i].e[1];
		if ( j > ftotd ) ftotd = j;
#ifndef MODULAR
		AE(ftotc, ct_abs(f[i].c));
#endif
	}
	fmind = ( d1 < d2 ? d1 : d2 );
	fmaxd = ( d1 > d2 ? d1 : d2 );
	minv = ( d1 < d2 ? 0 : 1);
	fmonic = 1;
	for ( i = 0 ; i < m ; i++ ) { if ( f[i].e[minv] == fmind && f[i].e[1-minv] ) fmonic = 0; }
	if ( ! fmonic && d1==d2 ) {
		minv = 1-minv;
		fmonic = 1;
		for ( i = 0 ; i < m ; i++ ) { if ( f[i].e[minv] == fmind && f[i].e[1-minv] ) fmonic = 0; }
		if ( ! fmonic ) minv = 1-minv;
	}
#if HAVE_COEFF_MAX
	if ( ftotc < 0 ) ftotc = COEFF_MAX;
#endif
	evec[MIND] = fmind;
	evec[MONIC] = 1 - fmonic;
	evec[MAXD] = fmaxd;
	evec[TOTD] = ftotd;
	evec[TERMS] = m;
#ifndef MODULAR
	evec[TOTC] = ftotc;
#else
	evec[TOTC] = 0;			// Not meaningful in this case...
#endif
}

void bipoly_evaluate_moniclike(bipoly *fs)
{
	int fmind, fmind_lp, minv;
	int lt_hi_power, lt_n_monom;
	int tt_hi_power, tt_n_monom;
	register int i, d1, d2, d1_lp, d2_lp;
	biterm *f = &fs->t[0];
	int m = fs->n;
	
	d1 = d2 = 0;
	d1_lp = d2_lp = MAX_DEGREE;
	for ( i = 0 ; i < m ; i++ ) {
		if ( f[i].e[0] > d1 ) d1 = f[i].e[0];
		if ( f[i].e[1] > d2 ) d2 = f[i].e[1];
		if ( f[i].e[0] < d1_lp ) d1_lp = f[i].e[0];
		if ( f[i].e[1] < d2_lp ) d2_lp = f[i].e[1];
	}
	fmind = ( d1 < d2 ? d1 : d2 );
	fmind_lp = ( d1 < d2 ? d1_lp : d2_lp );
	minv = ( d1 < d2 ? 0 : 1);
	lt_hi_power = 0, lt_n_monom = 0;
	tt_hi_power = 0, tt_n_monom = 0;
	for ( i = 0 ; i < m ; i++ ) {
		if ( f[i].e[minv] == fmind ) {
			lt_n_monom++;
			if ( f[i].e[1-minv] > lt_hi_power )
				lt_hi_power = f[i].e[1-minv];
		}
		if ( f[i].e[minv] == fmind_lp ) {
			tt_n_monom++;
			if ( f[i].e[1-minv] > tt_hi_power )
				tt_hi_power = f[i].e[1-minv];
		}
	}

	if (lt_n_monom > tt_n_monom) {
		lt_n_monom = tt_n_monom;
		lt_hi_power = tt_hi_power;
	} else if (lt_n_monom == tt_n_monom) {
		lt_hi_power = _min(lt_hi_power, tt_hi_power);
	}
	
	printf(VT_RED "moniclike: (%d,%d)" VT_RESET "\n", lt_hi_power, lt_n_monom);
}

/*
	If *same is TRUE then gs is assumed to be the same poly as in last call.
	If fs is not better than gs *same will be set to TRUE on the assumption that the caller will use the same gs next time.
*/
int bipoly_cmp (bipoly *fs, bipoly *gs, int *same)
{
	ct fevec[6];
	static ct gevec[6];

	bipoly_evaluate(fs, fevec);
	if (!*same)
		bipoly_evaluate(gs, gevec);
	
//printf("cmp: fmind=%d, gmind=%d, fmonic=%d, gmonic=%d, ftotd=%d, gtotd=%d, fmaxd=%d, gmaxd=%d, fterms=%d, gterms=%d, ftotc=%d, gtotc=%d\n",
//	  fmind, gmind, fmonic, gmonic, ftotd, gtotd, fmaxd, gmaxd, m, n, ftotc, gtotc);
	
	*same = 1;
	if ( fevec[MIND] < gevec[MIND] ) { printf("better mind %d vs %d\n", to_i(fevec[MIND]), to_i(gevec[MIND])); *same = 0; return -1; }
	if ( gevec[MIND] < fevec[MIND] ) return 1;
	if ( fevec[MONIC] < gevec[MONIC] ) { puts("monic"); *same = 0; return -1; }
	if ( gevec[MONIC] < fevec[MONIC] ) return 1;
//if ( m < n ) { printf ("fewer terms %d vs %d\n", m, n); return -1; }
//if ( m > n ) return 1;
	if ( fevec[MAXD] < gevec[MAXD] ) { printf ("better maxd %d vs %d\n", to_i(fevec[MAXD]), to_i(gevec[MAXD])); *same = 0; return -1; }
	if ( gevec[MAXD] < fevec[MAXD] ) return 1;
	if ( fevec[TOTD] < gevec[TOTD] ) { printf ("better totd %d vs %d\n", to_i(fevec[TOTD]), to_i(gevec[TOTD])); *same = 0; return -1; }
	if ( gevec[TOTD] < fevec[TOTD] ) return 1;
	if ( fevec[TERMS] < gevec[TERMS] ) { printf ("fewer terms %d vs %d\n", to_i(fevec[TERMS]), to_i(gevec[TERMS])); *same = 0; return -1; }
	if ( gevec[TERMS] < fevec[TERMS] ) return 1;
#ifndef MODULAR
	if ( fevec[TOTC] < gevec[TOTC] ) { printf ("smaller totc " CTF " vs " CTF "\n", O(fevec[TOTC]), O(gevec[TOTC])); *same = 0; return -1; }
	if ( gevec[TOTC] < fevec[TOTC] ) return 1;
#endif
	return 0;
}

// Compare two polynomials using their pre-computed evaluation vectors.
int bipoly_cmp_evec (ct fevec[6], ct gevec[6])
{
	int i;

	for (i = 0; i < 6; i ++) {
		if (fevec[i] < gevec[i]) return -1;
		if (fevec[i] > gevec[i]) return +1;
	}

	return 0;
}

// computes f(x+/-1,y) or f(x,y+/-1) depending on var and sign
void bipoly_trans(bipoly *gs, bipoly *fs, int var, int sign)
{
	register int i, j, k, e1, e2, m1, m2;
	register biterm *f = &fs->t[0];
	register biterm *g = &gs->t[0];
	int n = fs->n;
	register ct c;
	
	if ( var < 0 || var > 1 ) { fprintf (stderr, "bipoly_trans: var invalid, must be 0 or 1"); gs->n = 0; return; }
	m1 = m2 = 0;
	bipoly_degs(fs, m1, m2);

	if (m1 > MAX_DEGREE || m2 > MAX_DEGREE) {
		fprintf(stderr, "Exceeded MAX_DEGREE=%d in bipoly_trans\n", MAX_DEGREE);
		printf("Exceeded MAX_DEGREE=%d in bipoly_trans\n", MAX_DEGREE); fs->n = 0; return; }
	
	for ( i = 0 ; i <= m1 ; i++ ) for ( j = 0 ; j <= m2 ; j++ ) coeff[i][j] = from_i(0);
	// dup code to avoid branching in loop
	if ( var == 0 ) {
		if (sign < 0) {
			for ( k = 0 ; k < n ; k++ ) {
			// Have two versions of the subtract loop, one unrolled, one not. When ct is a long,
			//  these seem to make a significant time difference on my machine (20%), but not
			//  in a consistent manner (best if sub1 unrolled and sub2 not, why?)
#if 1
				e1 = f[k].e[0];  e2 = f[k].e[1]; c = f[k].c;
				for ( i = 0 ; i <= e1 ; i++ ) {
					if ( i&1 )
						SE(coeff[e1-i][e2], M(c,bc[e1][i]));
					else
						AE(coeff[e1-i][e2], M(c,bc[e1][i]));
				}
#else
				// unrolled version
				e1 = f[k].e[0];  e2 = f[k].e[1];  c = f[k].c;
				for ( i = 0 ; i < ((e1+1)>>1) ; i++ ) {
					AE(coeff[e1-(i<<1)][e2], M(c,bc[e1][i<<1]));
					SE(coeff[e1-((i<<1)+1)][e2], M(c,bc[e1][(i<<1)+1]));
				}
				if (!(e1&1)) AE(coeff[0][e2], c);
#endif
			}
		} else {
			for ( k = 0 ; k < n ; k++ ) {
				e1 = f[k].e[0];  e2 = f[k].e[1]; c = f[k].c;
				for ( i = 0 ; i <= e1 ; i++ )
					AE(coeff[e1-i][e2], M(c,bc[e1][i]));
			}
		}
	} else {
		if (sign < 0) {
			for ( k = 0 ; k < n ; k++ ) {
#if 1
				e1 = f[k].e[0];  e2 = f[k].e[1]; c = f[k].c;
				for ( i = 0 ; i <= e2 ; i++ ) {
					if ( i&1 )
						SE(coeff[e1][e2-i], M(c,bc[e2][i]));
					else
						AE(coeff[e1][e2-i], M(c,bc[e2][i]));
				}
#else
				// unrolled version
				e1 = f[k].e[0];  e2 = f[k].e[1];  c = f[k].c;
				for ( i = 0 ; i < ((e2+1)>>1) ; i++ ) {
					AE(coeff[e1][e2-(i<<1)], M(c,bc[e2][i<<1]));
					SE(coeff[e1][e2-((i<<1)+1)], M(c,bc[e2][(i<<1)+1]));
				}
				if (!(e2&1)) AE(coeff[e1][0], c);	
#endif
			}
		} else {
			for ( k = 0 ; k < n ; k++ ) {
				e1 = f[k].e[0];  e2 = f[k].e[1];  c = f[k].c;
				for ( i = 0 ; i <= e2 ; i++ )
					AE(coeff[e1][e2-i], M(c,bc[e2][i]));
			}		
		}
	}
	k = 0;
	for ( i = m1 ; i >= 0 ; i-- ) {
		for ( j = m2 ; j >= 0 ; j-- ) {
			if ( coeff[i][j] != 0 ) {
				g[k].e[0] = i;  g[k].e[1] = j;
				g[k].c = coeff[i][j];
				k++;
			}
		}
	}
	gs->n = k;
}

// computes f(x+/-y,y) or f(x,y+/-x) depending on var and sign
void bipoly_skew(bipoly *gs, bipoly *fs, int var, int sign)
{
	register int i, j, k, e1, e2, m;
	register biterm *f = &fs->t[0];
	register biterm *g = &gs->t[0];
	register ct c;
	int n = fs->n;
	
	if ( var < 0 || var > 1 ) { fprintf (stderr, "bipoly_skew: var invalid, must be 0 or 1"); gs->n = 0; return; }
	m = 0;
	for (i = 0; i < n; i++) {
		j = f[i].e[0] + f[i].e[1];
		if (j > m) m = j;
	}
	// If the result would exceed MAX_DEGREE, print warnings everywhere and return 0.
	if (m > MAX_DEGREE) { printf ("Exceeded MAX_DEGREE=%d in bipoly_skew\n", MAX_DEGREE);
			fprintf (stderr, "Exceeded MAX_DEGREE=%d in bipoly_skew\n", MAX_DEGREE); gs->n = 0; return; }

	for ( i = 0 ; i <= m ; i++ ) for ( j = 0 ; j <= m-i ; j++ ) coeff[i][j] = 0;

	if ( var == 0 ) {
		if (sign < 0) {
			for ( k = 0 ; k < n ; k++ ) {
				e1 = f[k].e[0];  e2 = f[k].e[1];  c = f[k].c;
				for ( i = 0 ; i <= e1 ; i++ ) {
					if ( i&1 )
						SE(coeff[e1-i][e2+i], M(c,bc[e1][i]));
					else
						AE(coeff[e1-i][e2+i], M(c,bc[e1][i]));
				}
			}
		} else {
			for ( k = 0 ; k < n ; k++ ) {
				e1 = f[k].e[0];  e2 = f[k].e[1];  c = f[k].c;
				for ( i = 0 ; i <= e1 ; i++ )
					AE(coeff[e1-i][e2+i], M(c,bc[e1][i]));
			}		
		}
	} else {
		if (sign < 0) {
			for ( k = 0 ; k < n ; k++ ) {
				e1 = f[k].e[0];  e2 = f[k].e[1];  c = f[k].c;
				for ( i = 0 ; i <= e2 ; i++ ) {
					if ( i&1 )
						SE(coeff[e1+i][e2-i], M(c,bc[e2][i]));
					else
						AE(coeff[e1+i][e2-i], M(c,bc[e2][i]));
				}
			}
		} else {
			for ( k = 0 ; k < n ; k++ ) {
				e1 = f[k].e[0];  e2 = f[k].e[1];  c = f[k].c;
				for ( i = 0 ; i <= e2 ; i++ )
					AE(coeff[e1+i][e2-i], M(c,bc[e2][i]));
			}		
		}
	}
	k = 0;
	for ( i = m ; i >= 0 ; i-- ) {
		for ( j = m-i ; j >= 0 ; j-- ) {
			if ( coeff[i][j] != 0 ) {
				g[k].e[0] = i;  g[k].e[1] = j;
				g[k].c = coeff[i][j];
				k++;
			}
		}
	}
	gs->n = k;
}

void bipoly_sort (bipoly *fs, int var)
{
	register int i, j, k, m1, m2;
	register biterm *f = &fs->t[0];
	int n = fs->n;
	
	m1 = m2 = 0;
	bipoly_degs(fs, m1, m2);
	if (m1 > MAX_DEGREE || m2 > MAX_DEGREE) {
		fprintf(stderr, "Exceeded MAX_DEGREE=%d in bipoly_sort\n", MAX_DEGREE);
		printf("Exceeded MAX_DEGREE=%d in bipoly_sort\n", MAX_DEGREE); fs->n = 0; return; }
	
	for ( i = 0 ; i <= m1 ; i++ ) for ( j = 0 ; j <= m2 ; j++ ) coeff[i][j] = from_i(0);
	for ( k = 0 ; k < n ; k++ ) coeff[f[k].e[0]][f[k].e[1]] = f[k].c;
	k = 0;
	if ( var == 0 ) {
		for ( i = m1 ; i >= 0 ; i-- ) {
			for ( j = m2 ; j >= 0 ; j-- ) {
				if ( coeff[i][j] != 0 ) {
					f[k].e[0] = i;  f[k].e[1] = j;
					f[k].c = coeff[i][j];
					k++;
				}
			}
		}
	} else {
		for ( j = m2 ; j >= 0 ; j-- ) {
			for ( i = m1 ; i >= 0 ; i-- ) {
				if ( coeff[i][j] != 0 ) {
					f[k].e[0] = i;  f[k].e[1] = j;
					f[k].c = coeff[i][j];
					k++;
				}
			}
		}
	}
}


void bipoly_flip(bipoly *gs, bipoly *fs, int var)
{
	int i, m;
	biterm *f = &fs->t[0];
	int n = fs->n;
	biterm *g = &gs->t[0];

	if ( var < 0 || var > 1 ) { fprintf (stderr, "bipoly_flip: var invalid, must be 0 or 1"); return; }
	
	if ( n <= 1 ) { gs->n = n; gs->t[0] = fs->t[0]; return; }
	m= 0;
	for ( i = 0 ; i < n ; i++ ) { if ( f[i].e[var] > m ) m = f[i].e[var]; }
	if (f == g) {
		for ( i = 0 ; i < n ; i++ ) { g[i].e[var] = m-g[i].e[var]; }
	} else {
		for ( i = 0 ; i < n ; i++ ) {
			g[i] = f[i];
			g[i].e[var] = m-f[i].e[var];
		}
		gs->n = n;
	}	
}


void bipoly_sep(bipoly *gs, bipoly *fs, int var)
{
	int i, j, m;
	biterm *f = &fs->t[0];
	biterm *g = &gs->t[0];
	int n = fs->n;
		
	if ( n <= 1 )  { gs->n = n; gs->t[0] = fs->t[0]; return; }
	if ( var < 0 || var > 1 ) { fprintf (stderr, "invalid var, must be 0 or 1"); return; }
	m= 0;
	for ( i = 0 ; i < n ; i++ ) { j = f[i].e[0]+f[i].e[1]; if ( j > m ) m = j; }
	if ( m > MAX_DEGREE ) { printf ("Exceeded MAX_DEGREE=%d in bipoly_sep\n", MAX_DEGREE);
				fprintf (stderr, "Exceeded MAX_DEGREE=%d in bipoly_sep\n", MAX_DEGREE); gs->n = 0; return; }
	if (f == g) {
		for ( i = 0 ; i < n ; i++ ) {
			j = g[i].e[0]+g[i].e[1];
			g[i].e[var] = m-j;
		}	
	} else {
		for ( i = 0 ; i < n ; i++ ) {
			g[i] = f[i];
			j = f[i].e[0]+f[i].e[1];
			g[i].e[var] = m-j;
		}
		gs->n = n;
	}
}

void bipoly_swap (bipoly *gs, bipoly *fs)
{
	register int d, i;
	register biterm *f = &fs->t[0], *g = &gs->t[0];
	register int n = fs->n;

	if (f == g) {
		for ( i = 0 ; i < n ; i++ ) { d = f[i].e[0]; f[i].e[0] = f[i].e[1]; f[i].e[1] = d; }	
	} else {
		for ( i = 0 ; i < n ; i++ ) {
			g[i].c = f[i].c;
			g[i].e[0] = f[i].e[1];
			g[i].e[1] = f[i].e[0];
		}
		gs->n = n;
	}
}

// assumes poly sorted on second variable
int bipoly_negcnt (bipoly *fs)
{
#if MODULAR
	return 0;
#else
	register int i, m, d, psign;
	biterm *f = &fs->t[0];
	int n = fs->n;
	
	d = -1;  psign = 1;
	for ( m = i = 0 ; i < n ; i++ ) {
		if ( f[i].e[1] == d ) {
			if ( psign*ct_sign(f[i].c) < 0 ) m++;
		} else {
			if ( ct_sign(f[i].c) < 0 ) m++;
			d = f[i].e[1];
			psign = ct_sign(f[i].c);
		}
	}
	return m;
#endif
}

void bipoly_neg (bipoly *gs, bipoly *fs, int var)
{
	register int i;
	register biterm *f = &fs->t[0], *g = &gs->t[0];
	int n = fs->n;
	
	if ( var < 0 || var > 1 ) { fprintf (stderr, "invalid var, must be 0 or 1"); return; }
	if (f == g) {
		for ( i = 0 ; i < n ; i++ ) if ( f[i].e[var]&1 ) { MSE(f[i].c, -1); }
	} else {
		for ( i = 0 ; i < n ; i++ ) {
//			g[i].c = (f[i].e[var]&1) ? -f[i].c : f[i].c;	// This doesn't work well with GMP
			if (f[i].e[var]&1) {
				g[i].c = MS(f[i].c, -1);
			} else {
				g[i].c = f[i].c;
			}
			g[i].e[0] = f[i].e[0];
			g[i].e[1] = f[i].e[1];
		}
		gs->n = n;
	}
}

// Tidy up a polynomial by making the second variable the one with minimum degree
//  and reducing the number of negative signs. Append any moves taken to the
//  given move sequence (cannot be NULL)
int bipoly_normalize(bipoly *fs, int *move_seq, int move_seq_len)
{
	int i, m0, m1, m2, m3, msl = move_seq_len;
	int swap;
	char w1, w2;
	biterm *f = &fs->t[0];
	int n = fs->n;
	
	// Make the second variable the one with minimum degree
bipoly_degs(fs, m1, m2);
	swap = 0;
	if ( m1 < m2 ) {
		swap = 1;
		m2 = m1;
	} else if ( m1 == m2 ) {	
	//  If the polynomial is non-monic and has the same degree in both variables, 
	//   try to make it monic by swapping variables.
		for ( i = 0 ; i < n ; i++ ) if ( f[i].e[1] == m2 && f[i].e[0] ) break;
		if ( i < n ) swap = 1;
	}
	w1 = 'U';  w2 = 'V';
	if ( swap ) { puts("Swapping variables");  move_seq[msl++] = CMD_SWAP; bipoly_swap (fs,fs); }
	bipoly_sort(fs, 1);
	if ( ct_sign(f[0].c) < 0 ) {
		puts("Multiplying by -1");
		for ( i = 0 ; i < n ; i++ ) { MSE(f[i].c, -1); }
	}
	
	// fiddle with signs - see if negating one or both signs reduces the number of negative signs
	bipoly_copy(&bp_tmp[1],fs);
	bipoly_neg(&bp_tmp[1],&bp_tmp[1],0);
	if ( ct_sign(bp_tmp[1].t[0].c) < 0 ) for ( i = 0 ; i < bp_tmp[1].n ; i++ ) { MSE(bp_tmp[1].t[i].c, -1); }
	bipoly_copy(&bp_tmp[2],fs);
	bipoly_neg(&bp_tmp[2],&bp_tmp[2],1);
	if ( ct_sign(bp_tmp[2].t[0].c) < 0 ) for ( i = 0 ; i < bp_tmp[2].n ; i++ ) { MSE(bp_tmp[2].t[i].c, -1); }
	bipoly_copy(&bp_tmp[3],fs);
	bipoly_neg(&bp_tmp[3],&bp_tmp[3],0);
	bipoly_neg(&bp_tmp[3],&bp_tmp[3],1);
	if ( ct_sign(bp_tmp[3].t[0].c) < 0 ) for ( i = 0 ; i < bp_tmp[3].n ; i++ ) { MSE(bp_tmp[3].t[i].c, -1); }
	m0 = bipoly_negcnt(fs);
	m1 = bipoly_negcnt(&bp_tmp[1]);
	m2 = bipoly_negcnt(&bp_tmp[2]);
	m3 = bipoly_negcnt(&bp_tmp[3]);
	if ( m1 < m0 && m1 <= m2 && m1 <= m3 ) {
		puts ("Negating x"); move_seq[msl++] = CMD_NEG1; bipoly_copy(fs,&bp_tmp[1]);
	} else if ( m2 < m0 && m2 <= m1 && m2 <= m3 ) {
		puts ("Negating y"); move_seq[msl++] = CMD_NEG2; bipoly_copy(fs,&bp_tmp[2]);
	} else if ( m3 < m0 && m3 <= m1 && m3 <= m2 ) {
		puts ("Negating x"); move_seq[msl++] = CMD_NEG1;
		puts ("Negating y"); move_seq[msl++] = CMD_NEG2;
		bipoly_copy(fs,&bp_tmp[3]);
	}
	
	return msl;
}

char str1n[MAX_STRING], str1d[MAX_STRING], str2n[MAX_STRING], str2d[MAX_STRING], tmp[MAX_STRING];

/*
	We make only a minimal attempt to simplify the rational functions obtained by composing the specified sequence.
	The resulting output will be much larger than necessary.  Use your favorite computer algebra system to
	simplify it (SAGE, Magma, Maple, Mathematica,...).
*/
void bipoly_print_map (int seq[], int n, char v1, char v2, char w1, char w2)
{
	register int i,j;
	register char *s, *t;
	
	str1n[0] = w1; str1n[1] = '\0'; str1d[0] = '\0'; str2n[0] = w2; str2n[1] = '\0'; str2d[0] = '\0';
	for ( i = n-1 ; i >= 0 ; i-- ) {
		switch ( seq[i] ) {
		case CMD_SEP1:
			if ( i && seq[i-1] == CMD_SEP1 ) { i--; continue; }	 		// skip undos
			if ( str1d[0] ) {
				if ( str2n[0] ) {
					for ( s = str2n ; *s && !issign(*s) ; s++ );  for ( t = str1d ; *t && !issign(*t) ; t++ );
					if ( *s ) {
						if ( *t ) sprintf(tmp, "(%s)*(%s)", str2n, str1d); else sprintf(tmp, "(%s)*%s",str2n, str1d);
					} else {
						if ( *t ) sprintf(tmp, "%s*(%s)", str2n, str1d); else sprintf(tmp, "%s*%s", str2n, str1d);
					}
					strcpy(str2n,tmp);
				} else {
					strcpy(str2n,str1d);
				}
			}
			if ( str1n[0] ) {
				if ( str2d[0] ) {
					for ( s = str2d ; *s && !issign(*s) ; s++ );  for ( t = str1n ; *t && !issign(*t) ; t++ );
					if ( *s ) {
						if ( *t ) sprintf(tmp, "(%s)*(%s)", str2d, str1n); else sprintf(tmp, "(%s)*%s",str2d, str1n);
					} else {
						if ( *t ) sprintf(tmp, "%s*(%s)", str2d, str1n); else sprintf(tmp, "%s*%s", str2d, str1n);
					}
					strcpy(str2d,tmp);
				} else {
					strcpy(str2d,str1n);
				}
			}
			strcpy(tmp,str1n); strcpy(str1n,str1d); strcpy(str1d,tmp);
			break;
		case CMD_SEP2:
			if ( i && seq[i-1] == CMD_SEP2 ) { i--; continue; }	 		// skip undos
			if ( str2d[0] ) {
				if ( str1n[0] ) {
					for ( s = str1n ; *s && !issign(*s) ; s++ );  for ( t = str2d ; *t && !issign(*t) ; t++ );
					if ( *s ) {
						if ( *t ) sprintf(tmp, "(%s)*(%s)", str1n, str2d); else sprintf(tmp, "(%s)*%s",str1n, str2d);
					} else {
						if ( *t ) sprintf(tmp, "%s*(%s)", str1n, str2d); else sprintf(tmp, "%s*%s", str1n, str2d);
					}
					strcpy(str1n,tmp);
				} else {
					strcpy(str1n,str2d);
				}
			}
			if ( str2n[0] ) {
				if ( str1d[0] ) {
					for ( s = str1d ; *s && !issign(*s) ; s++ );  for ( t = str2n ; *t && !issign(*t) ; t++ );
					if ( *s ) {
						if ( *t ) sprintf(tmp, "(%s)*(%s)", str1d, str2n); else sprintf(tmp, "(%s)*%s",str1d, str2n);
					} else {
						if ( *t ) sprintf(tmp, "%s*(%s)", str1d, str2n); else sprintf(tmp, "%s*%s", str1d, str2n);
					}
					strcpy(str1d,tmp);
				} else {
					strcpy(str1d,str2n);
				}
			}
			strcpy(tmp,str2n); strcpy(str2n,str2d); strcpy(str2d,tmp);
			break;
		case CMD_FLIP1:
			if ( i && seq[i-1] == CMD_FLIP1 ) { i--; continue; }	 			// skip undos
			strcpy(tmp,str1n); strcpy(str1n,str1d); strcpy(str1d,tmp);
			break;
		case CMD_FLIP2:
			if ( i && seq[i-1] == CMD_FLIP2 ) { i--; continue; }	 			// skip undos
			strcpy(tmp,str2n); strcpy(str2n,str2d); strcpy(str2d,tmp);
			break;
		case CMD_ADD1:
			if ( i && seq[i-1] == CMD_SUB1 ) { i--; continue; }	 		// skip undos
			if ( ! str1n[0] ) {
				sprintf(tmp,"1+%s", str1d);
			} else if ( ! str1d[0] ) {
				j=strlen(str1n);
				if ( j>1 && str1n[j-1]=='1' && str1n[j-2]=='-' ) {
					memcpy(tmp,str1n,j-2); tmp[j-2]='\0';
				} else {
					sprintf(tmp,"%s+1", str1n);
				}
			} else {
				sprintf(tmp,"%s+%s", str1n, str1d);
			}
			strcpy(str1n,tmp);
			break;
		case CMD_ADD2:
			if ( i && seq[i-1] == CMD_SUB2 ) { i--; continue; }	 		// skip undos
			if ( ! str2n[0] ) {
				sprintf(tmp,"1+%s", str2d);
			} else if ( ! str2d[0] ) {
				j=strlen(str2n);
				if ( j>1 && str2n[j-1]=='1' && str2n[j-2]=='-' ) {
					memcpy(tmp,str2n,j-2); tmp[j-2]='\0';
				} else {
					sprintf(tmp,"%s+1", str2n);
				}
			} else {
				sprintf(tmp,"%s+%s", str2n, str2d);
			}
			strcpy(str2n,tmp);
			break;
		case CMD_SUB1:
			if ( i && seq[i-1] == CMD_ADD1 ) { i--; continue; }	 		// skip undos
			for ( s = str1d ; *s && !issign(*s) ; s++ );
			if ( ! str1n[0] ) {
				if ( *s ) sprintf(tmp,"1-(%s)", str1d); else sprintf(tmp,"1-%s", str1d);
			} else if ( ! str1d[0] ) {
				j=strlen(str1n);
				if ( j>1 && str1n[j-1]=='1' && str1n[j-2]=='+' ) {
					memcpy(tmp,str1n,j-2); tmp[j-2]='\0';
				} else {
					sprintf(tmp,"%s-1", str1n);
				}
			} else {
				if ( *s ) sprintf(tmp,"%s-(%s)", str1n, str1d); else sprintf(tmp,"%s-%s", str1n, str1d);
			}
			strcpy(str1n,tmp);
			break;
		case CMD_SUB2:
			if ( i && seq[i-1] == CMD_ADD2 ) { i--; continue; }	 		// skip undos
			for ( s = str2d ; *s && !issign(*s) ; s++ );
			if ( ! str2n[0] ) {
				if ( *s ) sprintf(tmp,"1-(%s)", str2d); else sprintf(tmp,"1-%s", str2d);
			} else if ( ! str2d[0] ) {
				j=strlen(str2n);
				if ( j>1 && str2n[j-1]=='1' && str2n[j-2]=='+' ) {
					memcpy(tmp,str2n,j-2); tmp[j-2]='\0';
				} else {
					sprintf(tmp,"%s-1", str2n);
				}
			} else {
				if ( *s ) sprintf(tmp,"%s-(%s)", str2n, str2d); else sprintf(tmp,"%s-%s", str2n, str2d);
			}
			strcpy(str2n,tmp);
			break;
		case CMD_NEG1:
			if ( ! str1n[0] ) { strcpy(str1n,"-1"); } else { sprintf(tmp, "-(%s)", str1n); strcpy(str1n,tmp); }
			break;
		case CMD_NEG2:
			if ( ! str2n[0] ) { strcpy(str2n,"-1"); } else { sprintf(tmp, "-(%s)", str2n); strcpy(str2n,tmp); }
			break;
		case CMD_SWAP:
			strcpy (tmp, str1n); strcpy(str1n,str2n); strcpy(str2n,tmp);
			strcpy (tmp, str1d); strcpy(str1d,str2d); strcpy(str2d,tmp);
			break;
/*		case CMD_SKEW1P:
			sprintf(tmp, "(%s)*(%s)", str1d, str1n);
			sprintf(str1n, "((%s)*(%s))+((%s)*(%s))", str1n, str2d, str2n, str1d);
			strcpy(str1d, tmp);
			break; */
		default:
			fprintf(stderr, "Invalid cmd %d\n", seq[i]); return;
		}
//printf("cmd=%d\n", seq[i]);
//printf("%c = (%s) / (%s)\n", v1, str1n, str1d);
//printf("%c = (%s) / (%s)\n", v2, str2n, str2d);
	}
	if ( ! str1n[0] ) { str1n[0] = '1'; str1n[1] = '\0'; }
	if ( ! str1d[0] ) { str1d[0] = '1'; str1d[1] = '\0'; }	
	if ( ! str2n[0] ) { str2n[0] = '1'; str2n[1] = '\0'; }
	if ( ! str2d[0] ) { str2d[0] = '1'; str2d[1] = '\0'; }
	printf("%c = (%s) / (%s)\n", v1, str1n, str1d);
	printf("%c = (%s) / (%s)\n", v2, str2n, str2d);
}

void bipoly_print (bipoly *fs, char v1, char v2, int star, int printInfo, FILE *fp)
{
	int i, j;
	int m1, m2, m;
	char mstr[2];
	ct ac;
#ifndef MODULAR
	ct maxc, totc;
#endif
	biterm *f = &fs->t[0];
	int n = fs->n;
		
	if ( star ) { mstr[0]='*'; mstr[1]='\0'; } else mstr[0] = '\0';
	m = m1 = m2 = 0;
#ifndef MODULAR
	totc = maxc = 0;
#endif
	for ( i = 0 ; i < n ; i++ ) {
		j = f[i].e[0]+ f[i].e[1];
		if ( f[i].e[0] > m1 ) m1 = f[i].e[0];
		if ( f[i].e[1] > m2 ) m2 = f[i].e[1];
		if ( j > m ) m = j;
		if (i>0) fprintf (fp, " %c ", (ct_sign(f[i].c)<0?'-':'+')); else if (ct_sign(f[i].c)<0) fputc('-', fp);
		ac = ct_abs(f[i].c);
		if ( f[i].e[0] ) {
			if ( ac != 1 ) fprintf (fp, "" CTF "%s", O(ac), mstr); fprintf(fp, "%c", v1);  if ( f[i].e[0] > 1 ) fprintf(fp, "^%d", f[i].e[0]);
			if ( f[i].e[1] ) { fprintf(fp, "%s%c", mstr, v2);  if ( f[i].e[1] > 1 ) fprintf(fp, "^%d", f[i].e[1]); }
		} else if ( f[i].e[1] ) {
			if ( ac != 1 ) fprintf (fp, "" CTF "%s", O(ac),mstr); fprintf(fp, "%c", v2);  if ( f[i].e[1] > 1 ) fprintf(fp, "^%d", f[i].e[1]);
		} else {			
			fprintf(fp, "" CTF "", O(ac));
		}
#ifndef MODULAR
		totc += ac;
		if ( ac > maxc ) maxc = ac;
#endif
	}
	if (n == 0) fputc('0', fp);
	if (printInfo) {
		fputc('\n', fp);
#if HAVE_COEFF_MAX
		if ( totc < 0 ) totc = COEFF_MAX;
#endif
#ifndef MODULAR
		fprintf(fp, "Max %c degree %d, Max %c degree %d, Max degree %d, Terms %d, Totc " CTF ", Maxc " CTF "\n", v1, m1, v2, m2, m, n, O(totc), O(maxc));
#else
		fprintf(fp, "Max %c degree %d, Max %c degree %d, Max degree %d, Terms %d, Totc --, Maxc --\n", v1, m1, v2, m2, m, n);

#endif
//		bipoly_evaluate_moniclike(fs);
//		putchar('\n');
	}
}

// prints poly which is assumed to be sorted on the second variable
void bipoly_print_sorted2 (bipoly *fs, char v1, char v2, int star, int printInfo, FILE *fp)
{
	int i, j, d, inparen, psign = 0, startparen = 0;
	int m1, m2, m;
	char mstr[2];
	ct ac;
#ifndef MODULAR
	ct maxc, totc;
#endif
	biterm *f = &fs->t[0];
	int n = fs->n;
	
	if ( star ) { mstr[0]='*'; mstr[1]='\0'; } else mstr[0] = '\0';
	m = m1 = m2 = 0;
#ifndef MODULAR
	totc = maxc = 0;
#endif
	inparen = 0;  d = -1;
	for ( i = 0 ; i < n ; i++ ) {
		j = f[i].e[0]+ f[i].e[1];
		if ( f[i].e[0] > m1 ) m1 = f[i].e[0];
		if ( f[i].e[1] > m2 ) m2 = f[i].e[1];
		if ( j > m ) m = j;
		if ( f[i].e[1] != d ) {
			if ( inparen ) {
				if ( d > 1 ) {
					fprintf (fp, ")%s%c^%d", mstr, v2, d);
				} else if ( d== 1 ) {
					fprintf (fp, ")%s%c", mstr, v2);
				}
				inparen = 0;
			}
			if ( i < n-1 && f[i].e[1] && f[i+1].e[1] == f[i].e[1] ) {
				inparen = 1;  d = f[i].e[1];  startparen = 1;
				psign = ct_sign(f[i].c);
				if ( i > 0 ) fprintf (fp, " %c ", (psign<0?'-':'+')); else if ( psign<0 ) fputc('-',fp);
				fprintf(fp, "(");
			}
		}
		ac = ct_abs(f[i].c);
		if ( inparen ) {
			if ( ! startparen ) fprintf (fp, " %c ", (psign*ct_sign(f[i].c)<0?'-':'+')); 
			if ( f[i].e[0] ) {
				if ( ac != 1 ) fprintf (fp, "" CTF "%s", O(ac),mstr); fprintf(fp, "%c", v1);  if ( f[i].e[0] > 1 ) fprintf(fp, "^%d", f[i].e[0]);
			} else {
				fprintf(fp, "" CTF "",  O(ac));
			}
		} else {
			if ( i > 0 ) fprintf (fp, " %c ", (ct_sign(f[i].c)<0?'-':'+')); else if ( ct_sign(f[i].c) < 0 ) fputc('-',fp);
			if ( f[i].e[0] ) {
				if ( ac != 1 ) fprintf (fp, "" CTF "%s", O(ac),mstr); fprintf(fp, "%c", v1);  if ( f[i].e[0] > 1 ) fprintf(fp, "^%d", f[i].e[0]);
				if ( f[i].e[1] ) { fprintf(fp, "%s%c", mstr, v2);  if ( f[i].e[1] > 1 ) fprintf(fp, "^%d", f[i].e[1]); }
			} else if ( f[i].e[1] ) {
				if ( ac != 1 ) fprintf (fp, "" CTF "%s",O(ac),mstr); fprintf(fp, "%c", v2);  if ( f[i].e[1] > 1 ) fprintf(fp, "^%d", f[i].e[1]);
			} else {
				fprintf(fp, CTF, O(ac));
			}
		}
		startparen = 0;
#ifndef MODULAR
		totc += ac;
		if ( ac > maxc ) maxc = ac;
#endif
	}
	if ( inparen ) {
		if ( d > 1 ) {
			fprintf (fp, ")%s%c^%d", mstr, v2, d);
		} else if ( d== 1 ) {
			fprintf (fp, ")%s%c", mstr, v2);
		}
		inparen = 0;
	}
	if (n == 0) fputc('0', fp);
	if (printInfo) {
		fputc('\n', fp);
#if HAVE_COEFF_MAX
		if ( totc < 0 ) totc = COEFF_MAX;
#endif
#ifndef MODULAR
		fprintf(fp, "Max %c degree %d, %s, Max %c degree %d, Total degree %d, Terms %d, Totc " CTF ", Maxc " CTF "\n", v2, m2, (f[0].e[0]? "Not Monic": "Monic"), v1, m1, m, n, O(totc), O(maxc));
		fprintf (fp, "(%d,%d,%d,%d,%d," CTF ")\n", m2, (f[0].e[0]?1:0), m1, m, n, O(totc)); 
#else
		fprintf(fp, "Max %c degree %d, %s, Max %c degree %d, Total degree %d, Terms %d, Totc --, Maxc --\n", v2, m2, (f[0].e[0]? "Not Monic": "Monic"), v1, m1, m, n);
		fprintf (fp, "(%d,%d,%d,%d,%d,0)\n", m2, (f[0].e[0]?1:0), m1, m, n); 
#endif
	}
}


// No longer assumes v1 always preceeds v2 in each term
// we do not detect overflow if a coefficient is larger than the data type we've chosen for it (this could cause problems!)
void bipoly_scan (bipoly *fs, char *str, char v1, char v2)
{
	char (*termmapp)[MAX_DEGREE+1][MAX_DEGREE+1];
	char *s;
	ct c;
	int i, j, e1, e2, sign;
	biterm *f = &fs->t[0];
	
	termmapp = (char (*)[MAX_DEGREE+1][MAX_DEGREE+1]) calloc(1, (MAX_DEGREE + 1) * (MAX_DEGREE + 1));
	if (!termmapp) { fprintf(stderr, "Could not allocate memory for term map\n"); exit(1); }
        
	for ( i = 0 ; i <= MAX_DEGREE ; i++ ) for ( j = 0 ; j <= MAX_DEGREE ; j++ ) (*termmapp)[i][j] = 0;
	for ( s = str, i = 0 ; ; i++ ) {
		while ( *s && ! isnumeric(*s) && *s != v1 && *s != v2 ) s++;
		if ( ! *s ) break;
		if ( *s == v1 || *s == v2 ) {
			sign = 1;
			c = from_i(1);
		} else {
			sign = ( *s=='-' ? -1 : 1 );
			while ( *s && ! isdigit(*s) && *s != v1 && *s != v2 ) s++;
			if ( ! *s ) break;
			if ( isdigit(*s) ) {
				c = atoct(s); while (isdigit(*s)) s++;
			} else {
				c = from_i(1);
			}
		}
		while ( *s && *s != v1 && *s != v2 && ! issign(*s) ) s++;
//		if ( ! *s ) break;
		e1 = e2 = 0;
		if ( *s == v1 ) {
			if ( *(s+1) == '^' ) {
				s+=2;  e1 = atoi(s);  while (isdigit(*s)) s++;
			} else {
				e1 = 1;
				s++;
			}
			while ( *s && *s != v2 && ! issign(*s) ) s++;
		} else {
			e1 = 0;
		}
		if ( *s == v2 ) {
			if ( *(s+1) == '^' ) {
				s+=2;  e2 = atoi(s);  while (isdigit(*s)) s++;
			} else {
				e2 = 1;
				s++;
			}
			while ( *s && *s != v1 && ! issign(*s) ) s++;
		} else {
			e2 = 0;
		}
		if ( !e1 && (*s == v1) ) {
			if ( *(s+1) == '^' ) {
				s+=2;  e1 = atoi(s);  while (isdigit(*s)) s++;
			} else {
				e1 = 1;
				s++;
			}
			while ( *s && *s != v2 && ! issign(*s) ) s++;
		}
		if ( e1 > MAX_DEGREE || e2 > MAX_DEGREE ) { fprintf (stderr, "bipoly: input poly exceeded MAX_DEGREE = %d\n", MAX_DEGREE); exit(0); }
		if ( (*termmapp)[e1][e2] ) { fprintf (stderr, "bipoly: unsimplified input poly (or out of order variables): two terms with degree (%d,%d) detected\n", e1, e2); exit (0); }
		(*termmapp)[e1][e2] = 1;
		f[i].c = MS(c, sign);
		f[i].e[0] = (short)e1;
		f[i].e[1] = (short)e2;
	}
        
	free(termmapp);
        
	fs->n = i;
}

void bipoly_scan_sorted2(bipoly *fs, char *str, char v1, char v2)
{
	char (*termmapp)[MAX_DEGREE+1][MAX_DEGREE+1];
	char *s, *endp = NULL;
	ct c;
	int sign, psign=1;
	int i, j, e1, e2, pe2 = -1, inparen;
	biterm *f = &fs->t[0];
	
	termmapp = (char (*)[MAX_DEGREE+1][MAX_DEGREE+1]) calloc(1, (MAX_DEGREE + 1) * (MAX_DEGREE + 1));
	if (!termmapp) { fprintf(stderr, "bipoly: could not allocate memory for term map\n"); exit(1); }
        
	inparen = 0;
	
	for ( i = 0 ; i <= MAX_DEGREE ; i++ ) for ( j = 0 ; j <= MAX_DEGREE ; j++ ) (*termmapp)[i][j] = 0;
	for ( s = str, i = 0 ; ; ) {
		while ( *s && ! isnumeric(*s) && *s != v1 && *s != v2 && *s != '(' && *s != ')') s++;
		if ( ! *s ) break;
		
		if (*s == ')') {
			if (!inparen) { fprintf(stderr, "bipoly: extra ')'\n"); exit(0); }
			inparen = 0;
			s = endp; continue;
		}
		
		if ( *s == v1 || *s == v2 || *s == '(') {
			sign = 1;
			c = from_i(1);
		} else {
			sign = ( *s=='-' ? -1 : 1 );
			while ( *s && ! isdigit(*s) && *s != v1 && *s != v2 && *s != '(' && *s != ')') s++;
			if ( ! *s ) break;
			if ( isdigit(*s) ) {
				c = atoct(s); while (isdigit(*s)) s++;
			} else {
				c = from_i(1);
			}
		}
		while ( *s && *s != v1 && *s != v2 && ! issign(*s) && *s != '(' && *s != ')') s++;
		
		if (*s == '(') {
			if (inparen) { fprintf(stderr, "bipoly: nested parentheses in input poly not allowed\n"); exit(0); }
			psign = sign;
			endp = ++s;
			while ( *endp && *endp != ')' ) endp++;
			if ( !*endp ) { fprintf(stderr, "bipoly: unmatched parentheses in input poly\n"); exit(0); }
			
			while ( *endp && ! isdigit(*endp) && *endp != v1 && *endp != v2 && *endp != '(' && *endp != '(' && !issign(*endp)) endp++;
			pe2=0;
			if ( *endp == v1 ) { fprintf(stderr, "bipoly: input poly not sorted by second variable\n"); exit(0); }
			if ( *endp == v2 ) {
				if ( *(endp+1) == '^' ) {
					endp+=2;  pe2 = atoi(endp);  while (isdigit(*endp)) endp++;
				} else {
					pe2 = 1;
					endp++;
				}
			}
			inparen = 1;
			continue;
		}
		
		e1 = e2 = 0;

		if ( *s == v1 ) {
			if ( *(s+1) == '^' ) {
				s+=2;  e1 = atoi(s);  while (isdigit(*s)) s++;
			} else {
				e1 = 1;
				s++;
			}
			while ( *s && *s != v2 && ! issign(*s) && *s != '(' && *s != ')') s++;
		} else {
			e1 = 0;
		}
		if ( inparen ) {
			if ( *s == v2 ) { fprintf(stderr, "bipoly: Input poly not sorted by second variable\n"); exit(0); }
			e2 = pe2;
		} else if ( *s == v2 ) {
			if ( *(s+1) == '^' ) {
				s+=2;  e2 = atoi(s);  while (isdigit(*s)) s++;
			} else {
				e2 = 1;
				s++;
			}
			while ( *s && *s != v1 && ! issign(*s) ) s++;
		} else {
			e2 = 0;
		}
		if ( !inparen && !e1 && (*s == v1) ) {
			if ( *(s+1) == '^' ) {
				s+=2;  e1 = atoi(s);  while (isdigit(*s)) s++;
			} else {
				e1 = 1;
				s++;
			}
			while ( *s && *s != v2 && ! issign(*s) ) s++;
		}

		if ( inparen ) sign *= psign;
		if ( e1 > MAX_DEGREE || e2 > MAX_DEGREE ) { fprintf (stderr, "bipoly: Input poly exceeded MAX_DEGREE = %d\n", MAX_DEGREE); exit(0); }
		if ( (*termmapp)[e1][e2] ) { fprintf (stderr, "bipoly: Unsimplified input poly (or out of order variables) two terms with degree (%d,%d) detected\n", e1, e2); exit (0); }
		(*termmapp)[e1][e2] = 1;
		f[i].c = MS(c, sign);
		f[i].e[0] = (short)e1;
		f[i].e[1] = (short)e2;
		i++;
	}
        
	free(termmapp);
        
	fs->n = i;
}




bipoly map1n, map1d, map2n, map2d;

// Divide out the numerator and denominator by the largest monomial dividing both.
void bipolyrat_monomial_reduce(bipoly *num, bipoly *den)
{
	int dn[2], dd[2];
	
	bipoly_mon_divider(num, dn[0], dn[1]);  bipoly_mon_divider(den, dd[0], dd[1]);
	dn[0] = _min(dn[0], dd[0]);
	dn[1] = _min(dn[1], dd[1]);
	
	if (dn[0] > 0 || dn[1] > 0)
		{ bipoly_scale_exp(num, -dn[0], -dn[1]);  bipoly_scale_exp(den, -dn[0], -dn[1]); }
}

// Multiply (dstn/dstd) by (srcn/srcd)^power and store the result in (dstn/dstd)
//  Should be done more efficiently.
void bipolyrat_mult_exp(bipoly *dstn, bipoly *dstd, bipoly *srcn, bipoly *srcd, int power)
{
	int i;
	
	if (power > 0) {
		for (i = 0; i < power; i++) { bipoly_mult(dstn, srcn); bipoly_mult(dstd, srcd); }
	} else {
		for (i = 0; i < -power; i++) { bipoly_mult(dstn, srcd); bipoly_mult(dstd, srcn); }
	}
}

// Compute (map1n/map1d)^m1exp[J-1]*(map2n/map2d)^m2exp[J-1] and store in mapJn/mapJd for J = 1,2
//  Could be done more efficiently.
void bipolyrat_apply_powers(bipoly *m1n, bipoly *m1d, bipoly *m2n, bipoly *m2d, int (*m1exp)[2], int (*m2exp)[2])
{
	if (((*m1exp)[0] == 1) && ((*m1exp)[1] == 0) && ((*m2exp)[0] == 0) && ((*m2exp)[1] == 1))
		return;

	bipoly_copy(&bp_tmp[0], m1n);
	bipoly_copy(&bp_tmp[1], m1d);
	bipoly_copy(&bp_tmp[2], m2n);
	bipoly_copy(&bp_tmp[3], m2d);
	
	bipoly_set_const(m1n, COEFF_1);  bipoly_set_const(m1d, COEFF_1);
	bipoly_set_const(m2n, COEFF_1);  bipoly_set_const(m2d, COEFF_1);
	bipolyrat_mult_exp(m1n, m1d, &bp_tmp[0], &bp_tmp[1], (*m1exp)[0]);
	bipolyrat_mult_exp(m1n, m1d, &bp_tmp[2], &bp_tmp[3], (*m1exp)[1]);
	bipolyrat_mult_exp(m2n, m2d, &bp_tmp[0], &bp_tmp[1], (*m2exp)[0]);
	bipolyrat_mult_exp(m2n, m2d, &bp_tmp[2], &bp_tmp[3], (*m2exp)[1]);
}

void bipoly_print_map_reduced (int seq[], int n, char v1, char v2, char w1, char w2)
{
	register int i;
	int inSFSeq = 0, newISFS;
	int m1exp[2], m2exp[2];
	
	map1n.n = 1;
	map1n.t[0].e[0] = 1; map1n.t[0].e[1] = 0; map1n.t[0].c = COEFF_1;
	map2n.n = 1;
	map2n.t[0].e[0] = 0; map2n.t[0].e[1] = 1; map2n.t[0].c = COEFF_1;
	bipoly_set_const(&map1d, COEFF_1);
	bipoly_set_const(&map2d, COEFF_1);
	
	for ( i = n-1 ; i >= 0 ; i-- ) {
		newISFS = (seq[i] >= CMD_FLIP1) && (seq[i] <= CMD_SEP2);
		if (newISFS && !inSFSeq) {
			m1exp[0] = 1; m1exp[1] = 0;
			m2exp[0] = 0; m2exp[1] = 1;
		} else if (!newISFS && inSFSeq) {
			bipolyrat_apply_powers(&map1n, &map1d, &map2n, &map2d, &m1exp, &m2exp);
			if (!map1n.n || !map1d.n || !map2n.n || !map2d.n) break;    		// no terms means poly was too large
		}
				
		inSFSeq = newISFS;
	
		switch ( seq[i] ) {
		case CMD_SEP1:
			if ( i && seq[i-1] == CMD_SEP1 ) { i--;  continue; }	 		// skip undos
			m1exp[0] = -m1exp[0];  m1exp[1] = -m1exp[1];
			m2exp[0] += m1exp[0];  m2exp[1] += m1exp[1];
//			bipoly_xchg(&map1n, &map1d);						// simpler approach
//			bipoly_mult(&map2n, &map1n);
//			bipoly_mult(&map2d, &map1d);
			break;
		case CMD_SEP2:
			if ( i && seq[i-1] == CMD_SEP2 ) { i--; continue; }	 		// skip undos
			m2exp[0] = -m2exp[0];  m2exp[1] = -m2exp[1];
			m1exp[0] += m2exp[0];  m1exp[1] += m2exp[1];			
//			bipoly_xchg(&map2n, &map2d);
//			bipoly_mult(&map1n, &map2n);
//			bipoly_mult(&map1d, &map2d);
			break;
		case CMD_FLIP1:
			if ( i && seq[i-1] == CMD_FLIP1 ) { i--; continue; }	 		// skip undos
			m1exp[0] = -m1exp[0];  m1exp[1] = -m1exp[1];
			break;
		case CMD_FLIP2:
			if ( i && seq[i-1] == CMD_FLIP2 ) { i--; continue; }	 		// skip undos
			m2exp[0] = -m2exp[0];  m2exp[1] = -m2exp[1];
			break;
		case CMD_ADD1:
			if ( i && seq[i-1] == CMD_SUB1 ) { i--; continue; }	 		// skip undos
			bipoly_add(&map1n, &map1d);
			break;
		case CMD_ADD2:
			if ( i && seq[i-1] == CMD_SUB2 ) { i--; continue; }	 		// skip undos
			bipoly_add(&map2n, &map2d);
			break;
		case CMD_SUB1:
			if ( i && seq[i-1] == CMD_ADD1 ) { i--; continue; }	 		// skip undos
			bipoly_sub(&map1n, &map1d);
			break;
		case CMD_SUB2:
			if ( i && seq[i-1] == CMD_ADD2 ) { i--; continue; }	 		// skip undos
			bipoly_sub(&map2n, &map2d);
			break;
		case CMD_NEG1:
			bipoly_scale(&map1n, COEFF_NEG_1);
			break;
		case CMD_NEG2:
			bipoly_scale(&map2n, COEFF_NEG_1);
			break;
		case CMD_SWAP:
			bipoly_xchg(&map1n, &map2n);
			bipoly_xchg(&map1d, &map2d);
			break;
		case CMD_SKEW1P:
			bipoly_mult(&map1n, &map2d);
			bipoly_copy(&bp_tmp[0], &map1d);
			bipoly_mult(&bp_tmp[0], &map2n);
			bipoly_add(&map1n, &bp_tmp[0]);
			bipoly_mult(&map1d, &map2d);
			break;
		case CMD_SKEW2P:
			bipoly_mult(&map2n, &map1d);
			bipoly_copy(&bp_tmp[0], &map2d);
			bipoly_mult(&bp_tmp[0], &map1n);
			bipoly_add(&map2n, &bp_tmp[0]);
			bipoly_mult(&map2d, &map1d);
			break;
		case CMD_SKEW1M:
			bipoly_mult(&map1n, &map2d);
			bipoly_copy(&bp_tmp[0], &map1d);
			bipoly_mult(&bp_tmp[0], &map2n);
			bipoly_sub(&map1n, &bp_tmp[0]);
			bipoly_mult(&map1d, &map2d);
			break;
		case CMD_SKEW2M:
			bipoly_mult(&map2n, &map1d);
			bipoly_copy(&bp_tmp[0], &map2d);
			bipoly_mult(&bp_tmp[0], &map1n);
			bipoly_sub(&map2n, &bp_tmp[0]);
			bipoly_mult(&map2d, &map1d);
			break;
		default:
			fprintf(stderr, "Invalid cmd %d\n", seq[i]); return;
		}
	}
	
	if (inSFSeq)
		bipolyrat_apply_powers(&map1n, &map1d, &map2n, &map2d, &m1exp, &m2exp);

	bipolyrat_monomial_reduce(&map1n, &map1d);
	bipolyrat_monomial_reduce(&map2n, &map2d);

	bipoly_sort(&map1n, 1);
	bipoly_sort(&map1d, 1);
	bipoly_sort(&map2n, 1);
	bipoly_sort(&map2d, 1);
	
	if ((map1n.n == 0) || (map1d.n == 0)) {
		printf("%c = (could not evaluate, degree > MAX_DEGREE)\n", v1);
	} else {
		printf("%c = (", v1);
		bipoly_print_sorted2(&map1n, w1, w2, 1, 0, stdout);
		printf(") / (");
		bipoly_print_sorted2(&map1d, w1, w2, 1, 0, stdout);
		puts(")");
	}

	if ((map2n.n == 0) || (map2d.n == 0)) {
		printf("%c = (could not evaluate, degree > MAX_DEGREE)\n", v2);
	} else {	
		printf("%c = (", v2);
		bipoly_print_sorted2(&map2n, w1, w2, 1, 0, stdout);
		printf(") / (");
		bipoly_print_sorted2(&map2d, w1, w2, 1, 0, stdout);
		puts(")");
	}
}

#ifndef MODULAR
// Runs the point (iv1,iv2) through the map given by seq and n.
//  Puts the result in (*ov1, *ov2) and sets *undef to false if it is defined;
//  Otherwise sets *undef to true and sets ov1->n to the number of moves reached before we got an undefined value.
void bipoly_apply_map_to_point (int seq[], int n, const ctq *iv1, const ctq *iv2, ctq *ov1, ctq *ov2, int *undef)
{
	register int i;
	ct an, ad, bn, bd, tmp, g;
	
	an = iv1->n;  ad = iv1->d;
	bn = iv2->n;  bd = iv2->d;
	*undef=1;
	for ( i = n-1 ; i >= 0 ; i-- ) {
		ov1->n = from_i(n-i);
		switch ( seq[i] ) {
		case CMD_SEP1:
			if ( i && seq[i-1] == CMD_SEP1 ) { i--;  continue; }	 		// skip undos
			if (an==0) return;
			tmp=an;  an=ad;  ad=tmp;
			bn *= an;  bd *= ad;
			g = ct_gcd(bn, bd); bn /= g; bd /= g;
			break;
		case CMD_SEP2:
			if ( i && seq[i-1] == CMD_SEP2 ) { i--; continue; }	 		// skip undos
			if (bn==0) return;
			tmp=bn;  bn=bd;  bd=tmp;
			an *= bn;  ad *= bd;
			g = ct_gcd(an, ad); an /= g; ad /= g;
			break;
		case CMD_FLIP1:
			if ( i && seq[i-1] == CMD_FLIP1 ) { i--; continue; }	 		// skip undos
			if (an==0) return;
			tmp=an;  an=ad;  ad=tmp;
			break;
		case CMD_FLIP2:
			if ( i && seq[i-1] == CMD_FLIP2 ) { i--; continue; }	 		// skip undos
			if (bn==0) return;
			tmp=bn;  bn=bd;  bd=tmp;
			break;
		case CMD_ADD1:
			if ( i && seq[i-1] == CMD_SUB1 ) { i--; continue; }	 		// skip undos
			an += ad;
			break;
		case CMD_ADD2:
			if ( i && seq[i-1] == CMD_SUB2 ) { i--; continue; }	 		// skip undos
			bn += bd;
			break;
		case CMD_SUB1:
			if ( i && seq[i-1] == CMD_ADD1 ) { i--; continue; }	 		// skip undos
			an -= ad;
			break;
		case CMD_SUB2:
			if ( i && seq[i-1] == CMD_ADD2 ) { i--; continue; }	 		// skip undos
			bn -= bd;
			break;
		case CMD_NEG1:
			an = -an;
			break;
		case CMD_NEG2:
			bn = -bn;
			break;
		case CMD_SWAP:
			tmp=an;  an=bn;  bn=tmp;
			tmp=ad;  an=bd;  bd=tmp;
			break;
		default:
			fprintf(stderr, "Invalid cmd %d\n", seq[i]); return;
		}
	}
	
	*undef=0;
	ov1->n=an;  ov1->d=ad;
	ov2->n=bn;  ov2->d=bd;
}
#endif

// Exchange the contents of fs and gs.
void bipoly_xchg(bipoly *fs, bipoly *gs)
{
	register int i, n1, n2;
	biterm tmp;
	
	n1 = fs->n; n2 = gs->n;
	fs->n = n2; gs->n = n1;

	if (n1 < n2) {
		for (i = 0; i < n1; i++) { tmp = fs->t[i]; fs->t[i] = gs->t[i]; gs->t[i] = tmp; }
		for (i = n1; i < n2; i++) { fs->t[i] = gs->t[i]; } 
	} else {
		for (i = 0; i < n2; i++) { tmp = fs->t[i]; fs->t[i] = gs->t[i]; gs->t[i] = tmp; }
		for (i = n2; i < n1; i++) { gs->t[i] = fs->t[i]; } 
	}
}

// Compute fs+gs and store the result in fs.
void bipoly_add(bipoly *fs, bipoly *gs)
{
	register int k;
	register biterm *f = &fs->t[0], *g = &gs->t[0];
	int n1 = fs->n, n2 = gs->n, m;
	ct c;
	
	for ( k = 0 ; k < n1 ; k++ ) coeff[f[k].e[0]][f[k].e[1]] = from_i(0);
	for ( k = 0 ; k < n2 ; k++ ) coeff[g[k].e[0]][g[k].e[1]] = g[k].c;
	m = 0;
	for ( k = 0 ; k < n1 ; k++ ) {
		c = A(f[k].c, coeff[f[k].e[0]][f[k].e[1]]);  coeff[f[k].e[0]][f[k].e[1]] = from_i(0);
		if (c == 0) continue;
		f[m] = f[k];  f[m].c = c;
		m++;
	}
	
	for ( k = 0; k < n2; k++ ) {
		if ( coeff[g[k].e[0]][g[k].e[1]] != 0 ) f[m++] = g[k];
	}
	
	fs->n = m;
}

// Compute fs-gs and store the result in fs.
void bipoly_sub(bipoly *fs, bipoly *gs)
{
	register int i, j, k;
	register biterm *f = &fs->t[0], *g = &gs->t[0];
	int n1 = fs->n, n2 = gs->n, m;
	ct c;
	
	for ( k = 0 ; k < n1 ; k++ ) coeff[f[k].e[0]][f[k].e[1]] = from_i(0);
	for ( k = 0 ; k < n2 ; k++ ) coeff[g[k].e[0]][g[k].e[1]] = g[k].c;
	m = 0;
	for ( k = 0 ; k < n1 ; k++ ) {
		c = S(f[k].c, coeff[f[k].e[0]][f[k].e[1]]);  coeff[f[k].e[0]][f[k].e[1]] = from_i(0);
		if (c == 0) continue;
		f[m] = f[k];  f[m].c = c;
		m++;
	}
		
	for ( k = 0; k < n2; k++ ) {
		i = g[k].e[0];  j = g[k].e[1];
		if ( coeff[i][j] != 0 ) {
			f[m].e[0] = i;  f[m].e[1] = j;
			f[m].c = MS(coeff[i][j], -1);
			m++;
		}
	}
	
	fs->n = m;
}

// Compute fs*gs and store the result in fs.
void bipoly_mult(bipoly *fs, bipoly *gs)
{
	register int i, j, k, fm1, fm2, gm1, gm2, fe0, fe1;
	register biterm *f = &fs->t[0], *g = &gs->t[0];
	register ct fc;
	int n1 = fs->n, n2 = gs->n;
	
	// Compute maximum degrees and prepare coefficient table
	bipoly_degs(fs, fm1, fm2);
	bipoly_degs(gs, gm1, gm2);
	fm1 += gm1;
	fm2 += gm2;
	
	// If the product exceeds MAX_DEGREE, print warnings everywhere and return 0.
	if ((fm1 > MAX_DEGREE) || (fm2 > MAX_DEGREE)) { printf ("Exceeded MAX_DEGREE=%d in bipoly_mult\n", MAX_DEGREE);
				fprintf (stderr, "Exceeded MAX_DEGREE=%d in bipoly_mult\n", MAX_DEGREE); fs->n = 0; return; }


	for ( i = 0 ; i <= fm1 ; i++ ) for ( j = 0 ; j <= fm2 ; j++ ) coeff[i][j] = from_i(0);
	
	// Compute coefficients for fs*gs
	for (i = 0; i < n1; i++) {
		fe0 = f[i].e[0]; fe1 = f[i].e[1];
		fc = f[i].c;
		for (j = 0; j < n2; j++) AE(coeff[fe0+g[j].e[0]][fe1+g[j].e[1]], M(fc, g[j].c));
	}
	
	// Store result
	k = 0;
	for ( i = fm1 ; i >= 0 ; i-- ) {
		for ( j = fm2 ; j >= 0 ; j-- ) {
			if ( coeff[i][j] != 0 ) {
				f[k].e[0] = i;  f[k].e[1] = j;
				f[k].c = coeff[i][j];
				k++;
			}
		}
	}
	
	fs->n = k;
}

// Compute fs+c and store the result in fs.
void bipoly_add_const(bipoly *fs, ct c)
{
	int k, i, n = fs->n;
	biterm *f = &fs->t[0];
	
	if (c == 0)
		return;
	
	for (k = n-1; k >= 0; k--) {
		if ((f[k].e[0] == 0) && (f[k].e[1] == 0)) {
			AE(f[k].c, c);
			if (f[k].c == 0) {
				fs->n = n-1;
				for (i = k+1; i < n; i++) f[i-1] = f[i];	// could be improved...
			}
			return;
		}
	}
	
	f[n].c = c;
	f[n].e[0] = 0;
	f[n].e[1] = 0;
	
	fs->n = n + 1;
}

// Compute fs * s, where s is a scalar, and store the result in fs.
void bipoly_scale(bipoly *fs, ct s)
{
	int k, n = fs->n;
	biterm *f = &fs->t[0];
	
	for (k = 0; k < n; k++)
		ME(f[k].c, s);
}

// Compute fs * var1^e0 * var2^e1 and store the result in fs.
void bipoly_scale_exp(bipoly *fs, int e0, int e1)
{
	int k, n=fs->n;
	biterm *f = &fs->t[0];

	for (k = 0; k < n; k++) { f[k].e[0] += e0;  f[k].e[1] += e1; }
}

// Divide fs by the monomial var1^e1*var2^e2 for the largest possible e1 and e2.
//  (as long as fs has more than one term)
void bipoly_reduce_exps(bipoly *fs)
{
	int d0, d1;
	
	if (fs->n <= 1) return;
	bipoly_mon_divider(fs, d0, d1);
	if ((d0 != 0) || (d1 != 0))
		bipoly_scale_exp(fs, -d0, -d1);
}

// Decide whether fs and gs are equal, without assumptions on the orders of their terms.
int bipoly_equal(bipoly *fs, bipoly *gs)
{
	int n, k;
	biterm *f, *g;
	
	n = fs->n;  if (n != gs->n) return 0;
	f = &fs->t[0];  g = &gs->t[0];
		
	for ( k = 0 ; k < n ; k++ ) coeff[f[k].e[0]][f[k].e[1]] = f[k].c;
	for ( k = 0 ; k < n ; k++ ) {
		if (coeff[g[k].e[0]][g[k].e[1]] != g[k].c) return 0;
		coeff[g[k].e[0]][g[k].e[1]] = from_i(0);
	}
	for ( k = 0 ; k < n ; k++ ) if (coeff[f[k].e[0]][f[k].e[1]] != 0) return 0;
	
	return 1;
}

// ********* SA *********


// Return random float in the range [0, 1)
float randomfloat()
{
	return ((float) random()) / ((unsigned long) RAND_MAX+1);
}

#define SA_UPDATE_CLOCKS	(CLOCKS_PER_SEC * 240)		// Amount of time between updates of the SA algorithm
#define SA_DONE_TEMP		(0.000005)			// Temperature at which SA finishes (for sa_maxrad)
#define SA_TWO_OPS_PER_MOVE	(1)				// 1:  do one or two basic operations per move in SA
								// 0:  always do one basic operation per move in SA

#if RECORD_RAW_SA
int raw_sa_moves;
int raw_sa_move_seq[256000];
#endif

void bipoly_sa()
{
	clock_t start, last, cur;
	int whichop, whichop2, z, accepts, iters, same;
	float temp = 0;
	#define cand (&h)

	last = start = clock();
	temp = init_temp; z = 0; accepts = 0, iters = 0, same = 0;
#if RECORD_RAW_SA
	raw_sa_moves = 0;
#endif
	bipoly_cost_diff_config(&f);

	while (1) {
		// Every 500 loops, update temperature and check to see if it's time to
		//  output progress information or to stop.
		if (++z == 500) {
			z = 0; cur = clock();
//			temp = (init_temp / exp((cur - start) / time_parm));		// temperature based on time elapsed, with exponential decay
			temp = (init_temp * exp(-(float) iters / start_rate));		// temperature based on iterations, with exponential decay
			if ((cur - last) >= SA_UPDATE_CLOCKS) {
				last = cur;
				printf(VT_BOLD "--- accepts=%d (of %d), temp = %f, time = %lds\n" VT_RESET, accepts, iters, temp, (long) (cur - start) / CLOCKS_PER_SEC);
				bipoly_print(&f, var1, var2, 0, 1, stdout);
			}
			if (temp < SA_DONE_TEMP) break;
		}
		
		whichop = 1 + (((unsigned long) random()) % (CMD_TOP-1));
		whichop2 = 1 + (((unsigned long) random()) % (CMD_TOP-1));
		
		bipoly_copy(cand, &f);
		bipoly_exec_cmd(cand, cand, whichop);
		if (whichop2 != inverse_cmds_2[whichop])
			bipoly_exec_cmd(cand, cand, whichop2);
		else
			whichop2 = 0;

		if (exp(-(bias+bipoly_cost_diff(cand, &f, &same))/temp) > (randomfloat())) {
			bipoly_copy(&f, cand), accepts++, same=0;
#if RECORD_RAW_SA
			raw_sa_move_seq[raw_sa_moves++] = whichop;
			if (whichop2 != 0) raw_sa_move_seq[raw_sa_moves++] = whichop2;
#endif
		} else
			same=1;
		iters++;
		
		if (f.n == 1 && (f.t[0].e[0] == 0 || f.t[0].e[1] == 0)) break;	// If we just have x^n or y^n, we're done.
	}
	
	puts("SA done.");
	printf(VT_BOLD "*** accepts=%d (of %d), temp = %f, total time = %lds\n" VT_RESET, accepts, iters, temp, (long) (clock() - start) / CLOCKS_PER_SEC);
	bipoly_print(&f, var1, var2, 0, 1, stdout);

//	bipoly_interactive_try_map(raw_sa_move_seq, raw_sa_moves);
#undef cand
}

#ifndef MODULAR		// not yet converted to modular arithmetic
void bipoly_interactive_try_map(int move_seq[], int moves)
{
	ctq x0, y0, x, y;
	int undef;
	for (;;) {
		puts("Enter a point to run through the map:");
		x0.n = from_i(0); x0.d = from_i(1); y0.n = from_i(0); y0.d = from_i(1);
		printf ("x> ");
		(void) scanf("" CTFS "/" CTFS "", I(x0.n), I(x0.d));
		printf ("y> ");
		(void) scanf("" CTFS "/" CTFS "", I(y0.n), I(y0.d));

		printf("Map(" CTF "/" CTF "," CTF "/" CTF ") = ", O(x0.n), O(x0.d), O(y0.n), O(y0.d));
		bipoly_apply_map_to_point(move_seq, moves, &x0, &y0, &x, &y, &undef);
		
		if (undef) printf("undef -- after " CTF " moves\n", O(x.n));
		else printf("(" CTF "/" CTF ", " CTF "/" CTF ")\n", O(x.n), O(x.d), O(y.n), O(y.d));
	}
}
#endif

void bipoly_sa_loop()
{
	int whichop, i, same;
	float temp = 0, p;
	char buf[256];
	
	bipoly_cost_diff_config(&f);
	temp = init_temp;
	i = 0;
	while (1) {
		printf(VT_BOLD "Current poly: ");
		bipoly_print(&f, var1, var2, 0, 1, stdout);
		printf(VT_RESET);
		temp = (2*init_temp) / (++i);
		if (temp < 0.000001) temp = 0.000001;
		printf("Temperature=%f\n", (double) temp);

		printf ("cmd> ");  (void) fgets(buf,sizeof(buf),stdin); buf[strlen(buf)-1]='\0';
		if ( strcmp(buf, "end") == 0)
			break;
		
		whichop = atoi(buf);
		if ( whichop == 0 )
			whichop = 1 + (((unsigned long) random()) & 7);
		
		printf("Chose operation %d\n", whichop);
		bipoly_copy(&h, &f);
		bipoly_exec_cmd(&h, &h, whichop);
		same = 0;
		p = exp(-bipoly_cost_diff(&h, &f, &same)/temp);
		printf("Acceptance probability=%f    ", (double) p);
		if (p > (randomfloat())) {
			bipoly_copy(&f, &h);
			puts("(accepted)");
		} else
			puts("(rejected)");
	}
}

void bipoly_thresh()
{
	clock_t start, last, cur;
	int whichop, whichop2, z, accepts, iters, same;
	float temp = 0;
	#define cand (&h)

	last = start = clock();
	temp = init_temp; z = 0; accepts = 0, iters = 0, same = 0;
//	raw_sa_moves = 0;

	bipoly_cost_diff_config(&f);

	while (1) {
		// Every 500 loops, update temperature and check to see if it's time to
		//  output progress information or to stop.
		if (++z == 500) {
			z = 0; cur = clock();
//			temp = (init_temp / exp((cur - start) / time_parm));
			temp = (init_temp * exp(-(float) iters / start_rate));
			if ((cur - last) >= SA_UPDATE_CLOCKS) {
				last = cur;
				printf(VT_BOLD "--- accepts=%d (of %d), temp = %f, time = %lds\n" VT_RESET, accepts, iters, (double) temp, (long) (cur - start) / CLOCKS_PER_SEC);
				bipoly_print(&f, var1, var2, 0, 1, stdout);
				
			}
			if (temp < SA_DONE_TEMP) break;
		}
		
		whichop = 1 + (((unsigned long) random()) & 7);
		whichop2 = 1 + (((unsigned long) random()) & 7);
		
		bipoly_copy(cand, &f);
		bipoly_exec_cmd(cand, cand, whichop);
//		raw_sa_move_seq[raw_sa_moves++] = whichop;
		if (whichop2 != inverse_cmds_2[whichop])
			bipoly_exec_cmd(cand, cand, whichop2); // raw_sa_move_seq[raw_sa_moves++] = whichop2;

		if (bipoly_cost_diff(cand, &f, &same) < temp)
			bipoly_copy(&f, cand), accepts++, same=0;
		else
			same=1;
		iters++;
	}
	
	printf("TS done.");
	
#undef cand
}

// Run multiple SA runs with radius limit until an improved polynomial cannot be
//  found after a certain number of tries.
// Input and output poly is the global variable f.
// The output poly just before the call to bipoly_normalized is stored in h (used for the avoidlist).
void bipoly_sa_iter()
{
	clock_t start, end;
	int trial, i, same, old_sa_moves, path_len;
	int rad_expand=0, rate_expand=0;
	int fxc;
	
	if ( sait_no_expand == 0 ) {
		rad_expand = start_radius / (sait_max_trials - 1);
		if (rad_expand == 0) rad_expand = 1;
		rate_expand = start_rate / 2;
	}
		
	sa_moves = 0;
	bipoly_copy(&f0, &f);			// Save original poly in f0
	bipoly_copy(&f1, &f);			// f1 holds the best poly found so far
	
	start = clock();
	trial = 0; i = 1;
	for (; trial<sait_max_trials; ) {
		if (sa_moves + start_radius + trial*rad_expand > (MAX_SA_MOVES - END_MOVES)) {	// allow for bipoly_normalize
			printf("SA: number of moves near MAX_SA_MOVES, halting\n");
			break;
		}
		printf("Iteration %d, trial %d\n", i, trial+1);
		old_sa_moves = sa_moves;
		bipoly_sa_maxrad(&f, start_radius + trial*rad_expand, start_rate + trial*rate_expand);
		(void) bph_search(&bph_tab, bipoly_hash(&f1), NULL, &fxc);
		same = 0;
		// Accept new polynomial if lexicographically better than original
		//  or if the old one is listed as a previously-found minimum
		if ((fxc >= HIGH_PENALTY) || (bipoly_cmp(&f, &f1, &same) < 0))		// if (bipoly_cost_diff(&f1, &f, &same) > 0)
			trial=0, i++, bipoly_copy(&f1, &f);
		else
			trial++, sa_moves=old_sa_moves, bipoly_copy(&f, &f1);
	}
	end = clock();
	
	
	putchar('\n'); putchar('\n');
	printf("%d moves:\n", sa_moves);
	for (i = 0; i < sa_moves; i++) printf("%d ", sa_move_seq[i]);
	putchar('\n');

	sa_moves = bipoly_simplify_map(&f0, sa_move_seq, sa_moves);
	printf("After simplifying, ");
	printf("%d moves:\n", sa_moves);
	for (i = 0; i < sa_moves; i++) printf("%d ", sa_move_seq[i]);
	putchar('\n');
	
	bipoly_copy(&h, &f0);
	for (i = 0; i < sa_moves; i++) {
		bipoly_exec_cmd(&h, &h, sa_move_seq[i]);
	}
	bipoly_sort(&h, 1);
	bipoly_print_sorted2(&h, var1, var2, 0, 1, stdout);
	
	path_len = sa_moves;
	sa_moves = bipoly_normalize(&f, sa_move_seq, sa_moves);
	if (sa_moves < 50 && compute_maps)
		bipoly_print_map_reduced(sa_move_seq, sa_moves,var1,var2,'x','y');

	printf("Final polynomial:\n");
	bipoly_sort(&f, 1);
	if (output_sorted)
		bipoly_print_sorted2(&f, var1, var2, 0, 1, stdout);
	else
		bipoly_print(&f, var1, var2, 0, 1, stdout);
	printf("Time = %lds\n", (long) (end - start) / CLOCKS_PER_SEC);
}


// Probability of trying a backtrack, multiplied by RAND_MAX
#define btprob (RAND_MAX / 2)			// 0.5 probability

void bipoly_sa_maxrad(bipoly *fs, int max_radius, int rate)
{
	clock_t start, last, cur=0;
	int whichop=0, whichop2=0, y, z, accepts, iters, same, bts, dobacktrack=0;
	float temp = 0, penalty;
	int moves;
	long b=0, r=0;
	bipoly *ca;
#define cand1 (&h)
#define cand2 (&h2)

	bipoly_cost_diff_config(fs);

	last = start = clock();
	temp = init_temp; y = z = 0; accepts = 0, iters = 0;
	moves = 0, bts = 0, same = 0;
	
	if (max_radius >= MAX_SEARCH_DEPTH) {
		fprintf(stderr, "SAMR: radius %d too large for MAX_SEARCH_DEPTH=%d!\n", max_radius, MAX_SEARCH_DEPTH);
		return;
	}
	
	bipoly_copy(g[0], fs);
	
	printf("SAMR: radius %d, rate %d\n", max_radius, rate);
	while (1) {
		if (++y == sa_temp_chg_iters) {
			y = 0;
			// We use an exponential decay annealing schedule
			temp = (init_temp * exp(-(float) iters / rate));
//			#undef SA_DONE_TEMP
//			#define SA_DONE_TEMP (0.005)
//			temp = (init_temp * (rate * rate) / ((float) (rate + iters) * (float) (rate + iters)));
			if (temp < SA_DONE_TEMP) break;
		}
		// Every 500 loops, check to see if it's time to
		//  output progress information.
		if (++z == 500) {
			z = 0; cur = clock();
			if ((cur - last) >= SA_UPDATE_CLOCKS) {
				last = cur;
				printf("--- dist = %d, temp = %f, time elapsed = %lds, accepts=%d (of %d)\n", moves, (double) temp, (long) (cur - start) / CLOCKS_PER_SEC, accepts, iters);
				for (r = 0; r < moves; r++) printf("%d ", sa_move_seq[sa_moves + r]);
				putchar('\n');
				bipoly_print(g[moves], var1, var2, 0, 1, stdout);
			}
		}

		b = 0;
		dobacktrack = (moves >= (max_radius - 1)) || ((moves != 0) && (random() < btprob));
		if (dobacktrack) {
			bts++;						// Count backtracks
			b = 2 + (random() & 3);		// Decide how far to backtrack
//			b = 2 + (random() & 1);
			if (b > moves) b = moves;	
		} 
		
		// Do one or two basic operations.
		#if SA_TWO_OPS_PER_MOVE
		whichop = 1 + (((unsigned long) random()) % (CMD_TOP-1));
		whichop2 = 1 + (((unsigned long) random()) % (CMD_TOP-1));
		#else
		whichop = 1 + (((unsigned long) random()) % (CMD_TOP-1));
		whichop2 = inverse_cmds_2[whichop];
		#endif
		
		bipoly_exec_cmd(cand1, g[moves - b], whichop);
		if (whichop2 != inverse_cmds_2[whichop]) {
			bipoly_exec_cmd(cand2, cand1, whichop2);
			ca = cand2;
		} else {
			whichop2 = 0;
			ca = cand1;
		}

		penalty = 0.;
#if 0
		if (moves >= 8) {
			static int tmp_move_seq[MAX_SEARCH_DEPTH+1];
			int tmp_moves;
			
			for (r = 0; r < moves; r++) tmp_move_seq[r] = sa_move_seq[sa_moves+r];
			tmp_moves = moves;
			
			if (dobacktrack) tmp_moves -= b;
			tmp_move_seq[tmp_moves++] = whichop; 
			if (whichop2 != 0) tmp_move_seq[tmp_moves++] = whichop2;
			if (bipoly_count_nt(tmp_moves, tmp_move_seq) < (2*tmp_moves/5)) penalty = 100.;
		}
		//if (r > (sa_moves / 2)) { penalty = (float) (r - (sa_moves / 2)); penalty *= 100. / sa_moves; }
#endif
//		for (r = 0; r < sa_moves+moves; r++) fprintf(stderr, "%d ", sa_move_seq[r]); fprintf(stderr, "\n");

		// Decide whether we should accept this candidate.
		if (exp(-(penalty+bipoly_cost_diff(ca, g[moves], &same))/temp) > (randomfloat())) {
			if (g_mark_accepts) fprintf(stderr, " (accepted)\n"); g_mark_accepts = 0;
			accepts++;  same=0;
			if (dobacktrack) {
				moves -= b;
			}
			sa_move_seq[sa_moves + moves] = whichop; 
			moves++;
			bipoly_copy(g[moves], cand1);
			if (whichop2 != 0) {
				sa_move_seq[sa_moves + moves] = whichop2;
				moves++;
				bipoly_copy(g[moves], cand2);
			}
		} else {
			if (g_mark_accepts) fprintf(stderr, " (rejected)\n"); g_mark_accepts = 0;
			same=1;
		}
		iters++;
	}

	cur = clock();
	bipoly_copy(fs, g[moves]);
	printf("--- dist = %d, temp = %f, time elapsed = %lds, accepts=%d (of %d)\n", moves, (double) temp, (long) (cur - start) / CLOCKS_PER_SEC, accepts, iters);
	for (r = 0; r < moves; r++) printf("%d ", sa_move_seq[sa_moves + r]); putchar('\n');
	bipoly_print(fs, var1, var2, 0, 1, stdout);

	sa_moves += bipoly_sm_remove_redundant_ranges_pc(sa_move_seq + sa_moves, moves, g);
#undef cand
}

void bipoly_samr_demo(bipoly *fs, int max_radius, int rate)
{
	clock_t start, last, cur=0;
	int whichop=0, whichop2=0, y, z, accepts, iters, same, bts, dobacktrack=0;
	float temp = 0, penalty;
	int moves;
	long b=0, r=0;
	bipoly *ca;
#define cand1 (&h)
#define cand2 (&h2)
#define TEMP_SCALE 1000000000.
#define SAMR_DEMO_TWO_OPS_PER_MOVE 0

	bipoly_cost_diff_config(fs);

	last = start = clock();
	temp = init_temp; y = z = 0; accepts = 0, iters = 0;
	moves = 0, bts = 0, same = 0;
	
	if (max_radius >= MAX_SEARCH_DEPTH) {
		fprintf(stderr, "SAMR: radius %d too large for MAX_SEARCH_DEPTH=%d!\n", max_radius, MAX_SEARCH_DEPTH);
		return;
	}
	
	bipoly_copy(g[0], fs);
	
	bph_begin_insertions(log(init_temp / SA_DONE_TEMP) * rate + 10000);
	bph_append(&bph_itab, bipoly_hash(fs), (int) (-temp * TEMP_SCALE));
	
	printf("SAMR demo: radius %d, rate %d\n", max_radius, rate);
	while (1) {
		if (++y == sa_temp_chg_iters) {
			y = 0;
			// We use an exponential decay annealing schedule
			temp = (init_temp * exp(-(float) iters / rate));
			if (temp < SA_DONE_TEMP) break;
		}
		// Every 500 loops, check to see if it's time to
		//  output progress information.
		if (++z == 500) {
			z = 0; cur = clock();
			if ((cur - last) >= SA_UPDATE_CLOCKS) {
				last = cur;
				printf("--- dist = %d, temp = %f, time elapsed = %lds, accepts=%d (of %d)\n", moves, (double) temp, (long) (cur - start) / CLOCKS_PER_SEC, accepts, iters);
				for (r = 0; r < moves; r++) printf("%d ", sa_move_seq[sa_moves + r]);
				putchar('\n');
				bipoly_print(g[moves], var1, var2, 0, 1, stdout);
			}
		}

		b = 0;
		dobacktrack = (moves >= (max_radius - 1)) || ((moves != 0) && (random() < btprob));
		if (dobacktrack) {
			bts++;
			b = 2 + (random() & 3);
			if (b > moves) b = moves;	
		} 
		
		// Do one or two basic operations.
		#if SAMR_DEMO_TWO_OPS_PER_MOVE
		whichop = 1 + (((unsigned long) random()) % (CMD_TOP-1));
		whichop2 = 1 + (((unsigned long) random()) % (CMD_TOP-1));
		#else
		whichop = 1 + (((unsigned long) random()) % (CMD_TOP-1));
		whichop2 = inverse_cmds_2[whichop];
		#endif
		
		bipoly_exec_cmd(cand1, g[moves - b], whichop);
		bph_append(&bph_itab, bipoly_hash(cand1), (int) (-temp * TEMP_SCALE));
		
		if (whichop2 != inverse_cmds_2[whichop]) {
			bipoly_exec_cmd(cand2, cand1, whichop2);
			bph_append(&bph_itab, bipoly_hash(cand2), (int) (temp * TEMP_SCALE));
			ca = cand2;
		} else {
			whichop2 = 0;
			ca = cand1;
		}

		penalty = 0.;

		// Decide whether we should accept this candidate.
		if (exp(-(penalty+bipoly_cost_diff(ca, g[moves], &same))/temp) > (randomfloat())) {
			if (g_mark_accepts) fprintf(stderr, " (accepted)\n"); g_mark_accepts = 0;
			accepts++;  same=0;
			if (dobacktrack) {
				moves -= b;
			}
			sa_move_seq[sa_moves + moves] = whichop; 
			moves++;
			bipoly_copy(g[moves], cand1);
			if (whichop2 != 0) {
				sa_move_seq[sa_moves + moves] = whichop2;
				moves++;
				bipoly_copy(g[moves], cand2);
			}
		} else {
			if (g_mark_accepts) fprintf(stderr, " (rejected)\n"); g_mark_accepts = 0;
			same=1;
		}
		iters++;
	}

	cur = clock();
	bipoly_copy(fs, g[moves]);

	printf("--- dist = %d, temp = %f, time elapsed = %lds, accepts=%d (of %d)\n", moves, (double) temp, (long) (cur - start) / CLOCKS_PER_SEC, accepts, iters);	
	
//	sa_moves += bipoly_sm_remove_redundant_ranges_pc(sa_move_seq + sa_moves, moves, g);
	
	printf("\nPolynomials:\n");
	for (y = 0; y <= moves; y++) {
		printf("%2d: ", y);
		bipoly_print(g[y], var1, var2, 1, 0, stdout);
		printf("  (" HTF ")\n", bipoly_hash(g[y]));
	}

	printf("\nCosts:\n");
	for (y = 0; y <= moves; y++) {
		printf("%2d: %f\n", y, (double) bipoly_cost(g[y]));
	}

	bph_finish_insertions();		// do this after computing costs

	printf("\nFirst Temps:\n");
	for (y = 0; y <= moves; y++) {
		int t;
		
		if (!bph_search_first(&bph_tab, bipoly_hash(g[y]), NULL, &t))
			printf(" could not find\n");
		else
			printf("%2d: %f\n", y, (double) (-t) / TEMP_SCALE);
	}
	
	for (r = 0; r < moves; r++) printf("%d ", sa_move_seq[sa_moves + r]); putchar('\n');
	bipoly_print(fs, var1, var2, 0, 1, stdout);	
#undef cand
}

// Count the number of non-translations in move_seq. Pairs of flips that cancel
//  one another (as in flip1 sub2 flip1, for example) are treated as one non-translation
//  rather than 2.
int bipoly_count_nt (int moves, int move_seq[])
{
	int i, nt = 0;
	int checkF1 = 0, checkF2 = 0;
	
	for (i = 0; i < moves; i++) {
		if (move_seq[i] > CMD_SUB2) nt++;
		if (move_seq[i] == CMD_FLIP1) {
			if (checkF1) { checkF1 = 0; nt--; }
			else { checkF1 = 1; }
		}
		else if (move_seq[i] != CMD_ADD2 && move_seq[i] != CMD_SUB2) checkF1 = 0;
		if (move_seq[i] == CMD_FLIP2) {
			if (checkF2) { checkF2 = 0; nt--; }
			else { checkF2 = 1; }
		}
		else if (move_seq[i] != CMD_ADD1 && move_seq[i] != CMD_SUB1) checkF2 = 0;
	}

//	for (i = 0; i < moves; i++) fprintf(stderr, "%d ", move_seq[i]);
//	fprintf(stderr, " #--> (%d)\n", nt);
	
	return nt;
}

// Compute cost(fs). Just informational: the SA algorithms use bipoly_cost_diff instead.
float bipoly_cost(bipoly *fs)
{
	ct fevec[6];
	float cost = 0.;
#ifndef MODULAR
	ct fmaxc;
	int i;
#endif
	int func;
	int fxc;
	
	if (fs->n == 0) return INFINITY;	//  bipoly_sep returns zero when an error occurs (i.e. degrees too large)
	
	bipoly_evaluate(fs, fevec);
	(void) bph_search(&bph_tab, bipoly_hash(fs), NULL, &fxc);
	
	func = CD_Function;
	if (func == 12) func = random() & 1;
	
	switch (func) {
	// Cost function 0: minimum degree, then monic. If monic matches up, then maxd, totd, and terms.
//	case 0:
//		break;
	// Cost function 1: standard cost function, linear combination of differences in evaluation vector entries
	case 1:
		cost = to_i(fevec[MIND]);
		cost += to_i(fevec[MONIC]) * 0.3;
		cost += to_i(fevec[MAXD]) * 0.05 + to_i(fevec[TOTD]) * 0.047;
		cost += to_i(fevec[TERMS]) * 0.0046;
#ifndef MODULAR
		cost += to_i(fevec[TOTC]) * 0.00045;
#endif
		break;
	// Cost function 2: pseudo-lexicographic weighting. Only the most significant difference is evaluated.
//	case 2:
//		break;
	// Cost function 3: only attempts to minimize total degree
	case 3:
		cost = to_i(fevec[TOTD]);
#ifndef MODULAR
		cost += to_i(fevec[TOTC]) * 0.00045;
#endif
		break;
	// Cost function 4: same as #1 except with a greater weight for monic polynomials -- better for larger polynomials
	case 4:
		cost = to_i(fevec[MIND]);
		cost += to_i(fevec[MONIC]) * 0.9;
		cost += to_i(fevec[MAXD]) * 0.05 + to_i(fevec[TOTD]) * 0.047;
		cost += to_i(fevec[TERMS]) * 0.0046;
#ifndef MODULAR
		cost += to_i(fevec[TOTC]) * 0.00045;
#endif
		break;
	// Cost function 5: similar to #1 except that the degree weights are scaled to be (roughly) independent of the degrees themselves. The global weight factor (-w) scales the weight given to the degrees.
//	case 5:
//		break;
	// Cost function 6: same as #5 except that the -w switch scales the weight given to monic
//	case 6:
//		break;
	case 7:
		cost = to_i(fevec[MIND]);
		cost += to_i(fevec[MAXD]) * 0.05 + to_i(fevec[TOTD]) * 0.047;
		cost += to_i(fevec[TERMS]) * 0.0046;
		cost *= global_weightf;
		cost += to_i(fevec[MONIC]) * 0.3;
#ifndef MODULAR
		cost += to_i(fevec[TOTC]) * 0.00045;
#endif
		break;

	default:
		if (func <= 6)
			fprintf(stderr, "bipoly_cost: cost function %d not permitted\n", func);
		else {
			fprintf(stderr, "bipoly: unknown cost function (%d) specified\n", func);
			exit(0);
		}
		break;
	}
	cost += fxc * 0.1;
#if 1	
	if (to_i(fevec[TOTD]) > CD_totd_ceiling)
		cost += (to_i(fevec[TOTD]) - CD_totd_ceiling) * (to_i(fevec[TOTD]) - CD_totd_ceiling);
#endif


#ifndef MODULAR
	fmaxc = 0;
	for (i = 0; i < fs->n; i++) { if (ct_abs(fs->t[i].c) > fmaxc) fmaxc = ct_abs(fs->t[i].c); }

	// Strongly discourage large constants
	if (fmaxc > CD_maxc_discourage) cost += 100;
	
	// Disallow if constants even larger
	if (fmaxc > CD_maxc_disallow) return INFINITY;
#endif

//	printf("cost: fmind=%ld, fmonic=%ld, ftotd=%ld, fmaxd=%ld, fterms=%ld, ftotc=%ld\n",
//	  fevec[MIND], fevec[MONIC], fevec[TOTD], fevec[MAXD], fevec[TERMS], fevec[TOTC]);
//	printf("  ==> %f\n", (double) cost);

	return cost;
}

// Compute "cost(fs)-cost(gs)"
//  Not the same as bipoly_cost(fs) - bipoly_cost(gs). These cost functions
//  ensure that relatively unimportant comparisons (e.g. #terms, coeff size) don't
//  override the more essential comparisons.
float bipoly_cost_diff (bipoly *fs, bipoly *gs, int *same)
{
	ct fevec[6];
	static ct gevec[6];
	float diff, weightf;
#ifndef MODULAR
	ct a, fmaxc, gmaxc;
	int i;
#endif
	int func, q;
	
	int fxc;
	static int gxc;
	
	if (gs->n == 0) return -INFINITY;	// Zero is never a valid model for a curve. The function
	if (fs->n == 0) return INFINITY;	//  bipoly_sep returns zero when an error occurs (i.e. degrees too large)
	
	bipoly_evaluate(fs, fevec);
	(void) bph_search(&bph_tab, bipoly_hash(fs), NULL, &fxc);

	if ( ! *same ) { bipoly_evaluate(gs, gevec); (void) bph_search(&bph_tab, bipoly_hash(gs), NULL, &gxc); }
	
	
//	if ((fxc > 100000) && (gxc <= 100000)) return INFINITY;
//	if ((fxc <= 100000) && (gxc > 100000)) return -INFINITY;
	
	func = CD_Function;
	if (func == 12) func = random() & 1;
	
	switch (func) {
	// Cost function 0: minimum degree, then monic. If monic matches up, then maxd, totd, and terms.
	case 0:
		diff = (to_i(fevec[MIND]) - to_i(gevec[MIND]));
		if (to_i(gevec[MONIC]) != to_i(fevec[MONIC])) {
			diff += (to_i(fevec[MONIC]) - to_i(gevec[MONIC])) * 0.3;
		} else {
			diff += (to_i(fevec[MAXD]) - to_i(gevec[MAXD])) * 0.25 + (to_i(fevec[TOTD]) - to_i(gevec[TOTD])) * 0.25;
			
			q = (to_i(fevec[TERMS]) - to_i(gevec[TERMS]));
			q = _min(q, 10);
			q = _max(q, -10);
			diff += q * 0.005;
		}
		break;
	// Cost function 1: standard cost function, linear combination of differences in evaluation vector entries
	case 1:
		diff = (to_i(fevec[MIND]) - to_i(gevec[MIND]));
		diff += (to_i(fevec[MONIC]) - to_i(gevec[MONIC])) * 0.3;
		diff += (to_i(fevec[MAXD]) - to_i(gevec[MAXD])) * 0.05 + (to_i(fevec[TOTD]) - to_i(gevec[TOTD])) * 0.047;
			
		q = (to_i(fevec[TERMS]) - to_i(gevec[TERMS]));
		q = _min(q, 10);
		q = _max(q, -10);
		diff += q * 0.0046;
		break;
	// Cost function 2: pseudo-lexicographic weighting. Only the most significant difference is evaluated.
	case 2:
		diff = (to_i(fevec[MIND]) - to_i(gevec[MIND]));
		if (!diff) {
			diff = (to_i(fevec[MONIC]) - to_i(gevec[MONIC])) * 0.5;
			if (!diff) {
				diff = (to_i(fevec[MAXD]) - to_i(gevec[MAXD])) * 0.2;
				if (!diff) {
					diff = (to_i(fevec[TOTD]) - to_i(gevec[TOTD])) * 0.1;
					if (!diff) {
						q = (to_i(fevec[TERMS]) - to_i(gevec[TERMS]));
						q = _min(q, 10);
						q = _max(q, -10);
						diff += q * 0.008;
					}
				}
			}
		}
		break;
	// Cost function 3: only attempts to minimize total degree
	case 3:
		diff = to_i(fevec[TOTD]) - to_i(gevec[TOTD]);
		break;
	// Cost function 4: same as #1 except with a greater weight for monic polynomials -- better for larger polynomials
	case 4:
		diff = (to_i(fevec[MIND]) - to_i(gevec[MIND]));
		diff += (to_i(fevec[MONIC]) - to_i(gevec[MONIC])) * 0.9;
		diff += (to_i(fevec[MAXD]) - to_i(gevec[MAXD])) * 0.05 + (to_i(fevec[TOTD]) - to_i(gevec[TOTD])) * 0.047;
			
		q = (to_i(fevec[TERMS]) - to_i(gevec[TERMS]));
		q = _min(q, 10);
		q = _max(q, -10);
		diff += q * 0.0046;
		break;
	// Cost function 5: similar to #1 except that the degree weights are scaled to be (roughly) independent of the degrees themselves. The global weight factor (-w) scales the weight given to the degrees.
	case 5:
		// Scaling factor for the degrees
		weightf = global_weightf / (to_i(fevec[MIND]) + to_i(gevec[MIND]));
		
		diff = (to_i(fevec[MIND]) - to_i(gevec[MIND]));
		diff += (to_i(fevec[MAXD]) - to_i(gevec[MAXD])) * 0.05 + (to_i(fevec[TOTD]) - to_i(gevec[TOTD])) * 0.047;
		diff *= weightf;
		diff += (to_i(fevec[MONIC]) - to_i(gevec[MONIC])) * 0.3;
			
		q = (to_i(fevec[TERMS]) - to_i(gevec[TERMS]));
		q = _min(q, 10);
		q = _max(q, -10);
		diff += q * 0.0046;
		break;
	// Cost function 6: same as #5 except that the -w switch scales the weight given to monic
	case 6:
		// Scaling factor for the degrees
		weightf = 30. / (to_i(fevec[MIND]) + to_i(gevec[MIND]));
		
		diff = (to_i(fevec[MIND]) - to_i(gevec[MIND]));
		diff += (to_i(fevec[MAXD]) - to_i(gevec[MAXD])) * 0.05 + (to_i(fevec[TOTD]) - to_i(gevec[TOTD])) * 0.047;
		diff *= weightf;
		diff += (to_i(fevec[MONIC]) - to_i(gevec[MONIC])) * 0.3 * global_weightf;
			
		q = (to_i(fevec[TERMS]) - to_i(gevec[TERMS]));
		q = _min(q, 10);
		q = _max(q, -10);
		diff += q * 0.0046;
		break;
	// for demo purposes
	case 7:
		diff = (to_i(fevec[MIND]) - to_i(gevec[MIND]));
		diff += (to_i(fevec[MAXD]) - to_i(gevec[MAXD])) * 0.05 + (to_i(fevec[TOTD]) - to_i(gevec[TOTD])) * 0.047;
		diff *= global_weightf;
		diff += (to_i(fevec[MONIC]) - to_i(gevec[MONIC])) * 0.3;
			
		q = (to_i(fevec[TERMS]) - to_i(gevec[TERMS]));
		q = _min(q, 10);
		q = _max(q, -10);
		diff += q * 0.0046;
		break;
	default:
		fprintf(stderr, "bipoly: unknown cost function specified\n");
		exit(0);
		break;
	}
#ifndef MODULAR
	a = fevec[TOTC] - gevec[TOTC];
	if (a > 10) q = 10;
	else if (a < -10) q = -10;
	else q = to_i(a);
	
	diff += q * 0.00045;
#endif

	diff += (fxc - gxc) * 0.1;
#if 1	
	// Discourage total degrees that are too high
	if (to_i(fevec[TOTD]) > CD_totd_ceiling) diff += (to_i(fevec[TOTD]) - CD_totd_ceiling) * (to_i(fevec[TOTD]) - CD_totd_ceiling);
	if (to_i(gevec[TOTD]) > CD_totd_ceiling) diff -= (to_i(gevec[TOTD]) - CD_totd_ceiling) * (to_i(gevec[TOTD]) - CD_totd_ceiling);
#endif


#ifndef MODULAR
#if 1
	fmaxc = gmaxc = 0;
	for (i = 0; i < fs->n; i++) { if (ct_abs(fs->t[i].c) > fmaxc) fmaxc = ct_abs(fs->t[i].c); }
	for (i = 0; i < gs->n; i++) { if (ct_abs(gs->t[i].c) > gmaxc) gmaxc = ct_abs(gs->t[i].c); }

	// Strongly discourage large constants
	if (fmaxc > CD_maxc_discourage) diff += 100;
	if (gmaxc > CD_maxc_discourage) diff -= 100;
	
	// Disallow if constants even larger
	if ((fmaxc > CD_maxc_disallow) && (gmaxc <= CD_maxc_disallow)) return INFINITY;
	if ((gmaxc > CD_maxc_disallow) && (fmaxc <= CD_maxc_disallow)) return -INFINITY;
#else
	// Strongly discourage large values of totc
	if (fevec[TOTC] > CD_maxc_discourage) diff += 100;
	if (gevec[TOTC] > CD_maxc_discourage) diff -= 100;
	
	// Disallow if totc is even larger
	if ((fevec[TOTC] > CD_maxc_disallow) && (gevec[TOTC] <= CD_maxc_disallow)) return INFINITY;
	if ((gevec[TOTC] > CD_maxc_disallow) && (fevec[TOTC] <= CD_maxc_disallow)) return -INFINITY;
#endif
#endif

//	printf("diff: fmind=%ld, gmind=%ld, fmonic=%ld, gmonic=%ld, ftotd=%ld, gtotd=%ld, fmaxd=%ld, gmaxd=%ld, fterms=%ld, gterms=%ld, ftotc=%ld, gtotc=%ld\n",
//	  fevec[MIND], gevec[MIND], fevec[MONIC], gevec[MONIC], fevec[TOTD], gevec[TOTD], fevec[MAXD], gevec[MAXD], fevec[TERMS], gevec[TERMS], fevec[TOTC], gevec[TOTC]);
//	printf("  ==> %f\n", (double) diff);

	return diff;
}

void bipoly_cost_diff_config(bipoly *fs)
{
	ct evec[6];
#ifndef MODULAR
	ct maxc;
	int i;
#endif
	
	bipoly_evaluate(fs, evec);
#ifndef MODULAR
	maxc = from_i(0);
	for (i = 0; i < fs->n; i++) { if (ct_abs(fs->t[i].c) > maxc) maxc = ct_abs(fs->t[i].c); }

	// Raising the coefficient ceilings helps with some polynomials.
	//  -- perhaps make these ceilings configurable on the command line?
#if HAVE_COEFF_MAX
//	CD_maxc_discourage = (maxc > (COEFF_MAX / 10000)) ? (COEFF_MAX / 100) : (100 * maxc);
//	CD_maxc_disallow = (maxc > (COEFF_MAX / 1000)) ? COEFF_MAX : (1000 * maxc);
	CD_maxc_discourage = (maxc > (COEFF_MAX / 100000)) ? (COEFF_MAX / 100) : (1000 * maxc);
	CD_maxc_disallow = (maxc > (COEFF_MAX / 10000)) ? COEFF_MAX : (10000 * maxc);
#else
//	CD_maxc_discourage = (100 * maxc);
//	CD_maxc_disallow = (1000 * maxc);
	CD_maxc_discourage = (1000 * maxc);
	CD_maxc_disallow = (10000 * maxc);
#endif
#endif

	CD_totd_ceiling = to_i(evec[TOTD]) * 2;

//	fprintf(stderr, "Cost configuration:\n totd ceiling " CTF "; maxc discourage at " CTF ", disallow at " CTF "\n", CD_totd_ceiling, CD_maxc_discourage, CD_maxc_disallow);
}

// Run through the chain of moves. At each point along the way, check
//  whether the current polynomial is equal to any later
//  polynomials. If so, zero out the series of moves that lead to the
//  latest matching polynomial.
// Could be vastly improved using hashing.
void bipoly_sm_zero_redundant_ranges(bipoly *init, int move_seq[], int moves)
{
	int i, j, match;

//	bipoly_copy(&bp_tmp[0], init); printf("00: "), bipoly_print(&bp_tmp[0], var1, var2, 0, 0, stdout); 
//	for ( i = 0; i < moves ; i++) bipoly_exec_cmd(&bp_tmp[0], move_seq[i]), printf("\n%02d: ", i+1), bipoly_print(&bp_tmp[0], var1, var2, 0, 0, stdout);
//	putchar('\n');
	
	bipoly_copy(&bp_tmp[0], init);
	
	for ( i = 0 ; i < moves ; i++ ) {
		if (move_seq[i] == 0) continue;
		bipoly_copy(&bp_tmp[1], &bp_tmp[0]);

		for ( match = -1, j = i ; j < moves ; j++ ) {
			if (move_seq[j] != 0) bipoly_exec_cmd(&bp_tmp[1], &bp_tmp[1], move_seq[j]);
			if (bipoly_equal(&bp_tmp[0], &bp_tmp[1])) match = j;
		}
		// If match is nonzero, then poly after i moves is equal to poly after match+1 moves.
		//  Zero out all the unnecessary moves.
		if ( match != -1 ) {
			printf("Removing moves %d-%d\n", i, match);
			for ( j = i ; j <= match ; j++ ) move_seq[j] = 0;
			if ( match == (moves - 1)) break;
			i = match + 1;
		}
	
		bipoly_exec_cmd(&bp_tmp[0], &bp_tmp[0], move_seq[i]);
	}
}

// Same as bipoly_sm_zero_redundant_ranges, except that all the intermediate polynomials
//  have already been computed and the redundant moves are removed rather than
//  zeroed out. It returns the number of moves after removal.
int bipoly_sm_remove_redundant_ranges_pc(int move_seq[], int moves, bipoly *polys[])
{
	int i, j, off, match;
	
	for ( i = 0 ; i < moves ; i++ ) {
		if (move_seq[i] == 0) continue;

		for ( match = 0, j = i + 1 ; j <= moves ; j++ ) {
			if (bipoly_equal(polys[i], polys[j])) match = j;
		}
		// If match is nonzero, then poly after i moves is equal to poly after match moves.
		//  Zero out all the unnecessary moves.
		if ( match != 0 ) {
			for ( j = i ; j < match ; j++ ) move_seq[j] = 0;
		}
	}
	
	off = 0;
	for ( i = 0 ; i < moves ; i++ ) {
		if ( move_seq[i] != 0)
			move_seq[off] = move_seq[i], off++;
	}
	
	return off;
}

// Remove moves from the given move_seq that have no effect
//   and remove pairs of moves that cancel one another.
// *init is the polynomial that the sequence of moves is supposed to
//   have started with.
int bipoly_simplify_map(bipoly *init, int move_seq[], int moves)
{
	int i, off, j, k;
	
	bipoly_copy(&bp_tmp[0], init);
	
//	bipoly_sm_zero_redundant_ranges(init, move_seq, moves);

	for ( i=0; i<(moves-1); i++) {
		k = move_seq[i];
		if (k == CMD_ADD1 || k == CMD_SUB1) {
			for ( j=i+1; j<moves; j++ ) {
				if (move_seq[j] == inverse_cmds[k]) {
					move_seq[i] = move_seq[j] = 0;
					break;
				} else if (move_seq[j] == CMD_FLIP1 || move_seq[j] == CMD_SEP1 || move_seq[j] == CMD_SEP2 || move_seq[j] >= CMD_TOP)
					break;
			}
		} else if (k == CMD_ADD2 || k == CMD_SUB2) {
			for ( j=i+1; j<moves; j++ ) {
				if (move_seq[j] == inverse_cmds[k]) {
					move_seq[i] = move_seq[j] = 0; break;
				} else if (move_seq[j] == CMD_FLIP2 || move_seq[j] == CMD_SEP1 || move_seq[j] == CMD_SEP2 || move_seq[j] >= CMD_TOP)
					break;
			}
		} else if (k == CMD_FLIP1) {
			for ( j=i+1; j<moves; j++ ) {
				if (move_seq[j] == CMD_FLIP1) {
					move_seq[i] = move_seq[j] = 0;
					break;
				} else if (move_seq[j] != CMD_ADD2 && move_seq[j] != CMD_SUB2 && move_seq[i] != 0)
					break;
			}
		} else if (k == CMD_FLIP2) {
			for ( j=i+1; j<moves; j++ ) {
				if (move_seq[j] == CMD_FLIP2) {
					move_seq[i] = move_seq[j] = 0;
					break;
				} else if (move_seq[j] != CMD_ADD1 && move_seq[j] != CMD_SUB1 && move_seq[i] != 0)
					break;
			}
		} else if (k == CMD_SEP1 || k == CMD_SEP2) {
			for ( j=i+1; j<moves; j++ ) {
				if (move_seq[j] == k) {
					move_seq[i] = move_seq[j] = 0;
					break;
				} else if (move_seq[i] != 0)
					break;
			}
		}
	}

	bipoly_sm_zero_redundant_ranges(init, move_seq, moves);

	off=0;
	for ( i=0 ; i<moves ; i++ ) {
		if ( move_seq[i] != 0)
			move_seq[off] = move_seq[i], off++;
	}
	
	return off;
}

#ifndef MODULAR
// Generic code for computing the gcd of two coefficients
ct ct_gcd(ct a, ct b)
{
	ct t;
	
	if (a < 0) a = -a;
	if (b < 0) b = -b;
	if (a < b) { t = a; a = b; b = t; }			// Let a be the largest.
	while (b != 0) {
		t = b;
		b = a % b;
		a = t;
	}
	return a;
}
#endif

// Initialize the binomial coefficient table
void bc_init ()
{
	register int i, j;
	
	bc[0][0] = from_i(1);
	for ( i = 1 ; i <= MAX_DEGREE ; i++ ) {
		for ( j = 0 ; j <= i/2 ; j++ ) {
			bc[i][j] = A(bc[i-1][j], ( j ? bc[i-1][j-1] : from_i(0) ));
			bc[i][i-j] = bc[i][j];
		}
	}
}

// Print out an evaluation vector.
void fprint_evec(FILE *fp, ct evec[6])
{
	fprintf(fp, "(" CTF "," CTF "," CTF "," CTF "," CTF "," CTF ")", O(evec[0]), O(evec[1]), O(evec[2]), O(evec[3]), O(evec[4]), O(evec[5]));
}

// General string conversion function.
ct atoct_gen(const char *str)
{
	ct z = COEFF_0;
	int n;
	char cc;
	const char *p = str;
	int sign = 0;
    
	for (;;) {
		cc = *p++; if (cc == '\0') break;
		if (isspace(cc)) continue;
		if (cc == '-') { sign = -1; continue; }
		if (cc == '+') continue;
		n = cc - '0';
		if (n >= 10) break;
		ME(z, 10);
		AE(z, n);
	}

	if (sign) ME(z, COEFF_NEG_1);

	return z;
}

#ifdef MODULAR
// Modular version of string conversion.
// Idea:
//  Let k=floor(log_10 P)
//  Collect k digits at a time; multiply n by 10^k (mod P) then add these k digits (mod P)
ct atomodp(const char *str)
{
	ct z = COEFF_0;
	unsigned int n, collect = 0, numc = 0;
	char cc;
	const char *p = str;
	int sign = 0;
	static const unsigned int pow10[] = {1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000};
    
	for (;;) {
		cc = *p++; if (cc == '\0') break;
		if (isspace(cc)) continue;
		if (cc == '-') { sign = -1; continue; }
		if (cc == '+') continue;
		n = cc - '0';
		if (n >= 10) break;
		numc++;
		collect = (collect * 10) + n;
		if (numc >= (P_digits - 1)) {
			ME(z, P_power10); AE(z, collect);
			numc = 0; collect = 0;
		}
	}
	if (numc > 0) { ME(z, pow10[numc]);  AE(z, collect); }
	if (sign) ME(z, COEFF_NEG_1);

	return z;
}

#endif

// Hash a bipoly.
ht bipoly_hash(const bipoly *fs)
{
	ht hash = 0;
	int i;
	const biterm *f = &fs->t[0];
	
	for (i = 0; i < fs->n; i++) {
		hash += (to_ht(f[i].c) + f[i].e[0]) * (((ht) f[i].e[0] * 139 + (ht) f[i].e[1] * 65521) * 100003 + 1);
	}
//		hash += f[i].c * (((ht) f[i].e[0] * 139 + (ht) f[i].e[1] * 6) + 1);				// Things of this kind make
//		hash += f[i].c * (((ht) f[i].e[0] * 139 + (ht) f[i].e[1] * 65521) * 104729 + 1);		//  the raw X1(N) polys hash to zero!!
	
	return hash;
}


inline void bph_update_exp(bpht *t)
{
	int n, e = -1;
	for ( n = t->size ; n != 0 ; n >>= 1 ) e++;
	t->exp = e;
}

#if 0
// Linear search. Slow!
int bph_search(bpht *t, ht n, int *index)
{
	int i;
	
//	if (bph_tab != NULL) {
		for (i = 0; i < t->size; i++) {
			if (t->e[i].hash == n) {
				printf("Got one: %d\n", i);
				if (index != NULL) *index = i;
				return t->e[i].data;
			}
		}
//	}
	
	if (index != NULL) *index = -1;
	return NULL;
}

#else
// Binary search.
int bph_search(bpht *t, ht n, int *index, int *data)
{
	int i, k;
	ht h = 0;
	
//	fprintf(stderr, "search: exp %d\n", t->exp);
//	for (i = 0 ; i < t->size; i++) fprintf(stderr, "        " HTF "\n", t->e[i].hash);
	
	if (t->size == 0) {
		i = 0;
	} else {
		k = (1 << t->exp);
		i = k - 1;
		do {	
			k >>= 1;
			h = t->e[i].hash;
			if (n == h) { if (index != NULL) *index = i; if (data != NULL) *data = t->e[i].data; return 1; }
			if (n < h)
				i -= k;
			else {
				i += k;
				while (i >= t->size && (k != 0)) { k >>= 1; i -= k;  }		// could be done better...
			}
		} while (k != 0);
		if (i == t->size) i--;
	}
	if (index != NULL) *index = i;
	if (data != NULL) *data = 0;
	return 0;
}
#endif

int bph_search_first(bpht *t, ht n, int *index, int *data)
{
	int i;

	if (!bph_search(t, n, &i, NULL))
		return 0;
	for (; i >= 0; i--)
		if (t->e[i].hash != n) { i++; break; }
	
	if (i < 0) i = 0;
	if (index != NULL) *index = i;
	if (data != NULL) *data = t->e[i].data;
	
	return 1;
}

void bph_alloc(bpht *t, int entries)
{
	al_entry *ntab;
	
	t->alloc += entries;
	ntab = (al_entry *) realloc(t->e, t->alloc * sizeof(al_entry));
	
	if (!ntab) { fprintf(stderr, "bipoly: could not allocate memory\n");  exit(0); }
	t->e = ntab;
	
	bph_update_exp(t);
}

// Compare function for qsort()'ing hash tables
int al_entry_compare(const void *a, const void *b)
{
	{
		ht ka, kb;
		
		ka = ((al_entry *) a)->hash;
		kb = ((al_entry *) b)->hash;
		if (ka < kb) return -1;
		if (ka > kb) return 1;	
	}
	{
		int da, db;
		
		da = ((al_entry *) a)->data;
		db = ((al_entry *) b)->data;
		if (da < db) return -1;
		if (da > db) return 1;	
	}
	
	return 0;
}

void bph_sort(bpht *t)
{
	qsort(&t->e[0], t->size, sizeof(al_entry), al_entry_compare);
}

void bph_append(bpht *t, ht hash, int data)
{	
	if (t->size >= t->alloc) { fprintf(stderr, "bipoly: no more room in avoidlist\n"); exit(0); }
	t->e[t->size].hash = hash;
	t->e[t->size].data = data;
	t->size++;

	bph_update_exp(t);
}

void bph_begin_insertions(int num_inserts)
{
	if (bph_itab.e != NULL) { fprintf(stderr, "bipoly: already in insertion phase\n"); exit(0); }
	
	bph_itab.e = (al_entry *) malloc(num_inserts * sizeof(al_entry));
	if (bph_itab.e == NULL) { fprintf(stderr, "bipoly: could not allocate memory\n"); exit(0); }
	
	bph_itab.alloc = num_inserts;
	bph_itab.size = 0;
	bph_update_exp(&bph_itab);
}

void bph_insert(ht hash, int data)
{
	int d, nd, res;
	int od;
	
//	fprintf(stderr, "insert: " HTF ", %d\n", hash, data);
	res = bph_search(&bph_tab, hash, &d, &od);
	if (res == 0) {
		res = bph_search(&bph_itab, hash, &d, &od);
		if (res == 0) {
			if (bph_itab.size >= bph_itab.alloc) { fprintf(stderr, "bipoly: no more room in avoidlist\n"); exit(0); }
			bph_itab.size++;
			nd = d;
			if (bph_itab.size > 1) {
				if (bph_itab.e[d].hash < hash) nd++;
				memcpy(&bph_itab.e[nd+1], &bph_itab.e[nd], (bph_itab.size - nd) * sizeof(al_entry));
			}
			bph_itab.e[nd].hash = hash;
			bph_itab.e[nd].data = data;
			
			bph_update_exp(&bph_itab);
		} else {
			bph_itab.e[d].data = _max(od, data);
		}
	} else {
		bph_tab.e[d].data = od + data;
	}
}

void bph_finish_insertions()
{
	int nal;

	if (bph_itab.e == NULL) { fprintf(stderr, "bipoly: not in insertion phase\n"); exit(0); }

	nal = bph_itab.size + bph_tab.size - bph_tab.alloc;
	if (nal > 0) bph_alloc(&bph_tab, nal);
	
	memcpy(&bph_tab.e[bph_tab.size], &bph_itab.e[0], bph_itab.size * sizeof(al_entry));
	bph_tab.size += bph_itab.size;
	bph_update_exp(&bph_tab);
	bph_sort(&bph_tab);	// xxx Just for testing
				// replace with code that adds the tables together. Then don't search through bph_tab in bph_search.
				//  This is not only faster but is also more correct.
	free(bph_itab.e);
	bph_itab.e = NULL;
	bph_itab.size = bph_itab.alloc = 0;


//	fprintf(stderr, "%x %d %d\n", (int) bph_tab.e, bph_tab.size, bph_tab.alloc);
//	fprintf(stderr, "%x %d %d\n", (int) bph_itab.e, bph_itab.size, bph_itab.alloc);
//	{
//		int i;
//		for (i = 0; i < bph_tab.size; i++) printf("" HTF " %d\n", bph_tab.e[i].hash, bph_tab.e[i].data);
//	}
}

void bph_load(const char *name)
{
	FILE *infp;
	int n, rets;
	ht h;
	
	infp = fopen (name,"r");
	if ( ! infp ) { fprintf (stderr, "bipoly: warning: avoidlist file %s does not exist\n", name); return; }

	for (;;) {
		rets = fscanf(infp, HTFS " %d\n", &h, &n);
		if (rets == EOF) break;
		if (rets != 2) { fprintf (stderr, "bipoly: parsing error in avoidlist file\n"); exit(0); }
		if (bph_tab.size >= bph_tab.alloc)
			bph_alloc(&bph_tab, 4096);
		
		bph_tab.e[bph_tab.size].hash = h;
		bph_tab.e[bph_tab.size].data = n;
		bph_tab.size++;
	}
	bph_update_exp(&bph_tab);

	fclose (infp);
}

void bph_save(const char *name, int append)
{
	FILE *outfp;
	int i;
	
	outfp = fopen(name, append ? "a" : "w");
	if (!outfp) { fprintf (stderr, "bipoly: error opening file %s\n", name); exit(0); }
	
	for (i = 0; i < bph_tab.size; i++) {
		fprintf(outfp, HTF " %d\n", bph_tab.e[i].hash, bph_tab.e[i].data);
	}
	
	fclose(outfp);
}

// Add a single poly to the avoidlist.
void bph_add_poly(const bipoly *fs, int penalty)
{
	bph_begin_insertions(1);
	bph_insert(bipoly_hash(fs), penalty);
	bph_finish_insertions();
}

void bph_add_poly_assocs(bipoly *fs, int penalty)
{
// This is OK for now. All these hashes could be easily computed at once and using much less memory inside bipoly_hash itself.

	bipoly_neg(&bp_tmp[1], fs, 0);
	bipoly_neg(&bp_tmp[2], fs, 1);
	bipoly_neg(&bp_tmp[3], &bp_tmp[1], 1);

	bph_insert(bipoly_hash(fs), penalty); 
	bph_insert(bipoly_hash(&bp_tmp[1]), penalty); 
	bph_insert(bipoly_hash(&bp_tmp[2]), penalty); 
	bph_insert(bipoly_hash(&bp_tmp[3]), penalty); 

	bipoly_copy(&bp_tmp[0], fs);
	bipoly_scale(&bp_tmp[0], COEFF_NEG_1);
	bipoly_scale(&bp_tmp[1], COEFF_NEG_1);
	bipoly_scale(&bp_tmp[2], COEFF_NEG_1);
	bipoly_scale(&bp_tmp[3], COEFF_NEG_1);

	bph_insert(bipoly_hash(&bp_tmp[0]), penalty); 
	bph_insert(bipoly_hash(&bp_tmp[1]), penalty); 
	bph_insert(bipoly_hash(&bp_tmp[2]), penalty); 
	bph_insert(bipoly_hash(&bp_tmp[3]), penalty); 

	bipoly_swap(&bp_tmp[0], &bp_tmp[0]);
	bipoly_swap(&bp_tmp[1], &bp_tmp[1]);
	bipoly_swap(&bp_tmp[2], &bp_tmp[2]);
	bipoly_swap(&bp_tmp[3], &bp_tmp[3]);

	bph_insert(bipoly_hash(&bp_tmp[0]), penalty); 
	bph_insert(bipoly_hash(&bp_tmp[1]), penalty); 
	bph_insert(bipoly_hash(&bp_tmp[2]), penalty); 
	bph_insert(bipoly_hash(&bp_tmp[3]), penalty); 

	bipoly_scale(&bp_tmp[0], COEFF_NEG_1);
	bipoly_scale(&bp_tmp[1], COEFF_NEG_1);
	bipoly_scale(&bp_tmp[2], COEFF_NEG_1);
	bipoly_scale(&bp_tmp[3], COEFF_NEG_1);

	bph_insert(bipoly_hash(&bp_tmp[0]), penalty); 
	bph_insert(bipoly_hash(&bp_tmp[1]), penalty); 
	bph_insert(bipoly_hash(&bp_tmp[2]), penalty); 
	bph_insert(bipoly_hash(&bp_tmp[3]), penalty); 
}

// Insert fs, and every polynomial within three moves of it, into the avoidlist.
void bph_add_poly_nbhd3(bipoly *fs)
{
	#define HLEVEL 20
	int j[4];
	
	bph_begin_insertions((1 + (8) + (8*8) + (8*8*8)) * 16);
	
	bph_add_poly_assocs(fs, HIGH_PENALTY);
	
	for ( j[0] = 1 ; j[0] < CMD_TOP ; j[0]++ ) {
		bipoly_exec_cmd(g[1],fs,j[0]);
		if ( !g[1]->n ) continue;
		bph_add_poly_assocs(g[1], HLEVEL * 3 * 3);
		for ( j[1] = 1 ; j[1] < CMD_TOP ; j[1]++ ) {
			if ( j[1] == inverse_cmds[j[0]] ) continue;
			bipoly_exec_cmd(g[2],g[1],j[1]);
			if ( !g[2]->n ) continue;
			bph_add_poly_assocs(g[2], HLEVEL * 2 * 2);
			for ( j[2] = 1 ; j[2] < CMD_TOP ;  j[2]++ ) {
				if ( j[2] == inverse_cmds[j[1]] ) continue;
				bipoly_exec_cmd(g[3],g[2],j[2]);
				if ( !g[3]->n ) continue;
				bph_add_poly_assocs(g[3], HLEVEL * 1 * 1);
			}
		}
	}
	
	bph_finish_insertions();
}

