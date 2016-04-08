/* COPYRIGHT 1997 Fred Hutchinson Cancer Research Center, Seattle, WA USA */
/*	blockmap.h	Used by block_vis.c               */

struct block_pos_struct {
  char code ;                        /* the block's code */
  int  start ;                       /* the block's start position */
  int  end ;                         /* the block's end position */
};
typedef struct block_pos_struct block_pos ;


struct seq_map_struct {
  char seq_name[SMALL_BUFF_LENGTH] ; /* sequence name */
  int seq_len ;                      /* sequence length */
  int num_blocks ;                   /* number of blocks in this sequence */
  block_pos *blocks ;                /* pointer to array of block positions */
};
typedef struct seq_map_struct sequence_map ;


struct blocks_map_struct {
  char block_family[SMALL_BUFF_LENGTH] ; /* the blocks family code */
  char description[SMALL_BUFF_LENGTH] ; /* the family description */
  char id[SMALL_BUFF_LENGTH] ;       /* the short family description */
  int num_seqs ;                     /* number of sequences */
  int tot_num_blocks ;               /* number of blocks */
  int max_seq_len ;                  /* longest sequence length */
  sequence_map *seq_map ;            /* pointer to array of sequence maps */
};
typedef struct blocks_map_struct blocks_map ;


#define MAXNAME 80	/* Maximum file name length */

#define BLOCKS_ALLOC 10 /* Initial number of block data structures 
                           to be allocated to blocks family array. */
#define OK    0

#define ERROR 1

