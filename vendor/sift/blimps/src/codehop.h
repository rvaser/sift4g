/* COPYRIGHT 1997 Fred Hutchinson Cancer Research Center
	codehop.h
   Header file for codehop.c
 8/18/01 Updated list of genetic codes
10/14/04 Updated list of genetic codes (see ../include/gcode.h)
==========================================================================*/
#include <sys/time.h>

#define MAXNAME 80	/* Maximum file name length */

#define BLOCKS_ALLOC 10 /* Initial number of block data structures 
                           to be allocated to blocks family array. */
#define REGIONS_ALLOC 20 /* Initial number of region structures 
                           to be allocated to core degenerate regions array. */
#define OK    0

#define ERROR 1

#define BLOCK_AC_LEN 8

#define COMMENTS "Comments: \"-\" instead of a input/output file name => stdin/stdout\n          The order of the arguments is not important but the\n          first name is for the input file and the second for the output.\npssm_type = 2 odds ratios normalized to 100\n\t3 pseudos, log odds, nats\n\t5 pseudos, log odds, half bits\n\t6 pseudos, log odds, third bits\n\t10 pseudos, odds ratios, nats\n\t20 just counts+pseudo counts\n\t30 average score\nGenetic code type = 0 Standard\n\t1 Vertebrate Mitochondrial\n\t2 Yeast Mitochondrial\n\t3 Mold Mitochondrial and Mycoplasma\n\t4 Invertebrate Mitochondrial\n\t5 Ciliate Nuclear\n\t6 Echinoderm Mitochondrial\n\t7 Euplotid Nuclear\n\t8 Bacterial and Plant Plastid\n\t9 Alternative Yeast Nuclear\n\t10 Ascidian Mitochondrial\n\t11 Flatworm Mitochondrial\n\t12 Blepharisma Macronuclear\n\t13 Chlorophycean Mitochondrial\n\t14 Trematode Mitochondrial\n\t15 Scenedesmus obliquus mitochondrial\n\t16 Thraustochytrium mitochondrial\nDegeneracy parameters=0.0 all nucleotides that actually appear are counted\n\t1.0 only nucleotide with highest value counted\n\tBetween 0 and 1 nucleotides with high values counted,\n\t if value/highest-value >= degeneracy parameter\nClamp temperature     Target melting temperature for clamp in degC \nConcentration         Probe concentration in nM\nRose restrictions     If set, uses Tim Rose's restrictions on boundaries\n                      of core degenerate region.\nMost common codons    If set, uses most common codons in clamp.\nBegin oligo           If set, oligo must start on a conserved column,\n                       otherwise core strictness is applied.\nApoly-x       Maxiumum number of consecutive nucleotides of same type.\n\n"

#define AA_FREQUENCY_FNAME "default.amino.frq"
#define ID_FREQUENCY_FNAME "identity.frq"

#define ALPHABET 14 /* position in structure array nt_adegen (ntbet.h) with
                       number (4) and chars (A,C,G,T) of non-degenerate nucleotides. */

#define PSSM_TYPE_DFLT 2

#define PSSM_TYPES " 2 3 5 6 10 20 30 " /* each type must be flanked by blanks ! */

#define PSSM_DEFS {"","","odds ratios normalized to 100","pseudos, log odds, nats","","pseudos, log odds, half bits","pseudos, log odds, third bits","","","","pseudos, odds ratios, nats","","","","","","","","","","just counts+pseudo counts","","","","","","","","","","average score"}

#define GCODE_TYPE_DFLT 0

#define CODON_USAGE_FILE_DFLT "default.codon.use"

#define MIN_STRICTNESS 0.0  /* these 2 values are by definition and */
#define MAX_STRICTNESS 1.0  /* should not be changed.               */

#define CORE_MIN_LEN 11   /* minimal length of core degenerate regions */
#define CORE_MAX_LEN 12   /* maximal length of core degenerate regions */

#define CLAMP_MIN_LEN 17    /* min length of clamp (non-degenerate) region*/
#define CLAMP_MAX_LEN 165    /* max length of clamp (non-degenerate) region*/

#define DEFAULT_CORE_DEG 128    /* default degeneracy of core region */
#define DEFAULT_TEMP 60.0	/* default clamp temperature in degC */
#define DEFAULT_CONC 50.0	/* default concentration in nM units */
#define DEFAULT_KONC 50		/* default salt concentration in mM units */
#define DEFAULT_POLYX 5		/* default max run of any nuc in clamp */
#define DEFAULT_PRODLEN 200	/* default product length */

#define INT_TO_LOG2 {DBL_MIN, 0., 1., 1.5850, 2.} /* a table to find log2(X) */

/* # of positions (4 chars each) per line */
#define WWW_MATRIX_PRINT_WIDTH 240
#define SHL_MATRIX_PRINT_WIDTH 18

#define HELP_REQUEST " ? help HELP usage USAGE "   /* flanking blanks required */


/***** Structure type definitions *****/

struct position_degeneracy_struct {
  int degeneracy;                    /* the number of residues deemed to 
                                        appear in the position */
  char *residues;                    /* a pointer to a cahr array with the 
                                        residues appearing in the position */
};
typedef struct position_degeneracy_struct PosDegen;


struct degeneracy_struct {
  double strictness;                 /* the strictness value used to
                                        calculate the degeneracy 
                                        (see procedure J_degeneracy) */
  int  length;                       /* number of positions in this array */
  PosDegen **positions;              /* pointer to array of position degeneracy
                                        structure pointers */
};
typedef struct degeneracy_struct Dgnrcy;

struct region_struct {
  int start;                         /* region start position */

  int len;                           /* region length */

  int degen;                         /* degeneracy of region */
};
typedef struct region_struct Region;

/*	Modified matrix structure for DNA PSSMs  */
struct dna_matrix {
   Block *block;
   int strand;
   int length;
   MatType *weights[4];		/* ACGT weights for each column */
   MatType *max_weight;		/* Max weight for each column */
   int *nres;			/* Number of actual different residues in col */
   int *nsres;		/* Num of diff residues in col using strictness */
   int *nt_core;		/* Pointer into nt_btoa for degenerate res */
   int *nt_clamp;		/* Pointer into nt_btoa for nondegenerate res */
   int *nt_common;		/* Pointer into nt_btoa for most common codon */
   int *nt_max;		/* Pointer into nt_btoa for max PSSM weight */
};

/*		Temporary sort structure  */
struct sort_array {
   int clump, number;
   double core_score, clamp_score;
   struct oligo_list *oligo;
};

/*		Sorted list of oligos  */
struct oligo_list {
   struct dna_matrix *pssm;	/*  Pointer to the DNA PSSM */
   int strand;			/*  Negative is complemented PSSM */
   int core_first;		/*  1st position of core region */
   int core_last;		/*  Last position of core region */
   double core_score;		/*  Degen. score of core region (log base 4) */
   int clamp_first;		/*  1st position of clamp region */
   int clamp_last;		/*  Last position of clamp region */
   double clamp_score;		/*  Degen. score of clamp region (log base 4) */
   double clamp_prob;		/*  Clamp probability score */
   double clamp_temp;		/*  Clamp temperature */
   double clamp_temp_nn;	/*  Nearest neighbor clamp temperature */
   int clump;			/*  Left to right clump number */
   int number;			/*  Left to right display number */
   int out_order;		/*  Printing order  */
   struct oligo_list *next_oligo;
};

