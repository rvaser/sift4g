struct algnmnt_data_struct {
  double score;                      /* the score between the matrices of the
                                        aligned blocks */
  int positionE0;                    /* the alignment end position in the first 
                                        matrix */
  int positionE1;                    /* the alignment end position in the second 
                                        matrix */
  int complength;                    /* the alignment width */
};
typedef struct algnmnt_data_struct algnmnt;


struct blocks_algnmnt_struct {
  int  reported_alignments;          /* how many reported alignments in the 
                                        structure */
  char matrix0_number[NUMBER_WIDTH]; /* the first matrix/block number(accession)*/
  char matrix1_number[NUMBER_WIDTH]; /* the second matrix/block number(accession)*/
  algnmnt *alignment;                /* pointer to array of alignment structures - 
                                        each pair of blocks might have a number of
					significant non-overlapping alignments. */
};
typedef struct blocks_algnmnt_struct BlocksAlgnmnt;


struct blocks_score_struct1 {
  double score;                      /* the score between the matrices of the*/
                                     /* aligned blocks */
  int complength;                    /* the alignment width */
  char matrix0_number[NUMBER_WIDTH]; /* the first matrix/block number(accession) */
  int positionE0;                    /* the alignment end position in the first 
                                        matrix */
  char matrix1_number[NUMBER_WIDTH]; /* the second matrix/block number(accession)*/
  int positionE1;                    /* the alignment end position in the second 
                                        matrix */
};
typedef struct blocks_score_struct1 BlocksScore1;


struct blocks_score_struct {
  int score;                         /* the score between the matrices of the*/
                                     /* aligned blocks */
  int complength;                    /* the alignment width */
  char matrix0_number[NUMBER_WIDTH]; /* the scanning matrix/block number (accession) */
  int positionS0;                    /* the alignment start position in the scanning matrix */
  int strength0;                     /* the scanning block/matrix strength */
  char matrix1_number[NUMBER_WIDTH]; /* the scanned matrix/block number (accession) */
  int positionS1;                    /* the alignment start position in the scanned matrix */
  int strength1;                     /* the scanned block/matrix strength */
};
typedef struct blocks_score_struct BlocksScore;


struct confidence_limits_struct {
  double low_score ;                 /* the low limit score */
  double high_score ;                /* the high limit score */
  double P ;                         /* the confidence */
};
typedef struct confidence_limits_struct Conf_limits ;


struct stat_struct {
  double mean ;                      /* the mean score */
  double variance ;                  /* the variance of the scores */
  double cutoff ;                    /* the score reporting cutoff */
};
typedef struct stat_struct Stat ;

struct ZvsPRCNTL_struct {
  double Z ;                         /* the Z score */
  double prcntl ;                    /* the percentile in scores from shuffled 
                                        unbiassed blocks*/
};
typedef struct ZvsPRCNTL_struct ZvsPRCNTL ;

struct line_struct {
  double intrcpt ;                /* the intercept (a, y=a+bx) */
  double slope ;                  /* the slope     (b, y=a+bx) */
  double std ;                    /* the standard deviation of y values  */
  double cutoff ;                 /* the score reporting cutoff */
};
typedef struct line_struct Line ;

#define AA_FREQUENCY_FNAME "default.amino.frq"
                                     /* file name of aa frequencies for 
					load_frequencies procedure */

#define MAXNAME 100 	             /* Maximum file name length */

#define MAXLINELEN 5000 	     /* Maximum line length */

#define ALGNMNTS_ALLOC 10            /* Initial number of alignment
                                        data structures allocated to
                                        alignment structure. */
#define OK    0

#define ERROR 1

#define WIDTHS 56                    /* block widths range. 0 can be used 
					for data from all the scores. 
					Currently the narrowest BLOCKS block 
					is 4 and the widest is 55. */ 

#define NARROWEST_WIDTH 4            /* default narrowest alignment to calculate */

#define Z_CUTOFF 5.6                   /* default Z score cutoff to report */

#define LOW_SCORE  5.6   /* scores with lower Z values are probably 
                            not significant in searching Blocks DB vs itself. */
#define HIGH_SCORE 8.3   /* scores with equal or higher Z values are probably 
                            significant in searching Blocks DB vs itself. */

#define PERSEARCHES 3179             /* number of searches the expected number 
                                        of scores is calculated for.
                                        BLOCKS 9.0 has 3179 blocks.*/

#define MATRIX_CMPRD_WIDTH 21        /* compared "width" (possible chars in
                                        each position) of matrices - 
					20 aa and gap. */

#define SMALL_VAL 0.000001

#define ENTRYSIZE 5000

#define BLOCK_DB 'B'

#define MATRIX_DB 'M'

#define ProDom_mul_DB 'P'

#define Conservative 'C'

#define Liberal 'L'

#define BLOCK_SEQUENCE_INCREASE_SIZE 20 /* default number of sequences to */
                                   /* allocate if the number of sequences */
                                   /* in the block is not specified */

#define POS_SUM 10000.0            /* sum of all counts in a matrix position */
                                   /* for matrix with probability values. */
