/* (C) Copyright 1993, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* ntbet.h: Definitions for nucleotides */
/* Modified by: Bill Alford */
/* Change log information is at the end of the file. */


#ifndef __NTBET_H__
#define __NTBET_H__

#define NUCID_MIN	0
#define NUCID_MAX	15
#define NUCID_CNT	16
#define NUCID_NAR	(NUCID_MAX+1)
#define NUCID_IGNORE	(NUCID_MAX+3)

#define NA NUCID_NAR	/* Not A Residue */
#define EL (NA+1)		/* End of Line */
#define ES (NA+2)		/* End of Sequence */
#define IC NUCID_IGNORE	/* Ignore this Character */
#define UNKNOWN_NT_CHR	'N'
#define GAP_NT_CHR '-'

/* Alpha chars that are not codes are treated as wildcard place holders */
/* (unknown nucleotides) */
#define UN 14			/* 14 = 'N' an unknown nt */


/*********************************************
*
*  ASCII-to-binary and binary-to-ASCII
*  translation tables
*
*********************************************/

EXTERN int nt_atob[1<<CHAR_BIT]
#ifdef INIT
	= {
	EL,NA,NA,NA,NA,NA,NA,NA,NA,IC,EL,NA,NA,EL,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	IC,NA,NA,NA,NA,NA,NA,NA,NA,NA,ES,NA,NA,15,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	ES, 0,10, 1,11,UN,UN, 2,12,UN,UN, 7,UN, 6,14,UN,
	UN,UN, 4, 9, 3, 3,13, 8,UN, 5,UN,NA,NA,NA,NA,NA,
	ES, 0,10, 1,11,UN,UN, 2,12,UN,UN, 7,UN, 6,14,UN,
	UN,UN, 4, 9, 3, 3,13, 8,UN, 5,UN,NA,NA,NA,NA,NA,

	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA
	}
#endif
	;

EXTERN char	nt_btoa[] INITIAL("ACGTRYMKWSBDHVN-");


/* Not used.

EXTERN char nt_n4toa[]
*/
/*                 111111 */
/*       0123456789012345 */
/*
INITIAL("-ACMGRSVTWYHKDBN");

EXTERN unsigned char nt_aton4[1<<CHAR_BIT]
#ifdef INIT
= {
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 1,14, 2,13, 0, 0, 4,11, 0, 0,12, 0, 3,15, 0,
     5,10, 5, 6, 8, 8, 7, 9, 0,10, 0, 0, 0, 0, 0, 0,
     0, 1,14, 2,13, 0, 0, 4,11, 0, 0,12, 0, 3,15, 0,
     5,10, 5, 6, 8, 8, 7, 9, 0,10, 0, 0, 0, 0, 0, 0,

     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
}
#endif
    ;
*/

/*********************************************
*
*  Reverse complement tables
*
*********************************************/

EXTERN char nt_brevcomp[128] /* binary-to-binary reverse complement */
#ifdef INIT
	= {
		 3, 2, 1, 0, 5, 4, 7, 6, 8, 9,13,12,11,10,14,15,
		16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,
		16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,
		16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,
		16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,
		16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,
		16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,
		16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16
	}
#endif
	;

EXTERN char nt_arevcomp[128] /* ASCII-to-ASCII reverse complement */
#ifdef INIT
#define NANT	'?' /* Not A Nucleotide, used when rev. complement is undefined */
	= {
        NANT,NANT,NANT,NANT,NANT,NANT,NANT,NANT,
        NANT,NANT,NANT,NANT,NANT,NANT,NANT,NANT,
        NANT,NANT,NANT,NANT,NANT,NANT,NANT,NANT,
        NANT,NANT,NANT,NANT,NANT,NANT,NANT,NANT,
        NANT,NANT,NANT,NANT,NANT,NANT,NANT,NANT,
        NANT,NANT,NANT,NANT,NANT, '-',NANT,NANT,
        NANT,NANT,NANT,NANT,NANT,NANT,NANT,NANT,
        NANT,NANT,NANT,NANT,NANT,NANT,NANT,NANT,
    /*   @    A    B    C    D    E    F    G  */
        NANT,'T', 'V', 'G', 'H', NANT,NANT,'C',
    /*   H    I    J    K    L    M    N    O  */
        'D', NANT,NANT,'M', NANT,'K', 'N', NANT,
    /*   P    Q    R    S    T    U    V    W  */
        NANT,NANT,'Y', 'S', 'A', 'A', 'B', 'W',
    /*   X    Y    Z                           */
        'N', 'R', NANT,NANT,NANT,NANT,NANT,NANT,
    /*   `    a    b    c    d    e    f    g  */
        NANT,'T', 'V', 'G', 'H', NANT,NANT,'C',
    /*   h    i    j    k    l    m    n    o  */
        'D', NANT,NANT,'M', NANT,'K', 'N', NANT,
    /*   p    q    r    s    t    u    v    w  */
        NANT,NANT,'Y', 'S', 'A', 'A', 'B', 'W',
    /*   x    y    z                           */
        'N', 'R', NANT,NANT,NANT,NANT,NANT,NANT,
        }
#undef NANT
#endif
        ;

EXTERN Degen	nt_adegen[16] /* ASCII degeneracy table */
#ifdef INIT
		= {
	{'A',	1,	"A"},
	{'C',	1,	"C"},
	{'G',	1,	"G"},
	{'T',	1,	"T"},
	{'R',	2,	"AG"},
	{'Y',	2,	"CT"},
	{'M',	2,	"AC"},
	{'K',	2,	"GT"},
	{'W',	2,	"AT"},
	{'S',	2,	"CG"},
	{'B',	3,	"CGT"},
	{'D',	3,	"AGT"},
	{'H',	3,	"ACT"},
	{'V',	3,	"ACG"},
	{'N',	4,	"ACGT"},
	{'-',	1,	"-"}
	}
#endif
	;

EXTERN Degen	nt_bdegen[16] /* binary degeneracy table */
#ifdef INIT
		= {
	{'\000',	1,	"\000"},
	{'\001',	1,	"\001"},
	{'\002',	1,	"\002"},
	{'\003',	1,	"\003"},
	{'\004',	2,	"\000\002"},
	{'\005',	2,	"\001\003"},
	{'\006',	2,	"\000\001"},
	{'\007',	2,	"\002\003"},
	{'\010',	2,	"\000\003"},
	{'\011',	2,	"\001\002"},
	{'\012',	3,	"\001\002\003"},
	{'\013',	3,	"\000\002\003"},
	{'\014',	3,	"\000\001\003"},
	{'\015',	3,	"\000\001\002"},
	{'\016',	4,	"\000\001\002\003"},
	{'\017',	1,	"\017"}
	}
#endif
	;


/* Generate an appropriate random BINARY nucleotide from the list ACGT */
#define RANDOM_BNUC(c)	(c < nt_atob['-'] ? \
	nt_bdegen[c].list[(rnd_gen()>>8)%nt_bdegen[c].ndegen] \
	: ((rnd_gen()>>8) % 4) )

/* Generate an appropriate random ASCII nucleotide from the list ACGT */
#define RANDOM_ANUC(c)	(nt_atob[c] < nt_atob['-'] ? \
	nt_adegen[nt_atob[c]].list[(rnd_gen()>>8)%nt_adegen[nt_atob[c]].ndegen] \
	: nt_btoa[((rnd_gen()>>8) % 4)] )

#undef NA
#undef EL
#undef ES
#undef IC

#undef UN

EXTERN double	ntfq[4]
#ifdef INIT
	= { 0.25, 0.25, 0.25, 0.25 }
#endif
	;

#endif /* __NTBET_H__ */




/* Change log information follows. */
/* 
 * Revision 0.52  1993/07/27  20:46:50  billa
 * Modified nt_atob so that alphabetic characters that are not codes are
 * treated as unknown nucleotides for placeholdings.  Also made it closer
 * to IUB codes by removing 'P' = 'R' and 'Q' = 'Y', 'P' and 'Q' are
 * treated as just mentioned above.  To take into account GCG codes, made
 * the representation of 'X' be equivalent to 'N' (this false under the
 * first sentence too).  The changes to nt_atob are:
 *      'P' gets the same binary representation as 'N', not 'R' as before
 *      'Q' gets the same binary representation as 'N', not 'Y' as before
 *      'E', 'F', 'I', 'J', 'L', 'O', 'Z' get the same binary
 *              representation as 'N'
 *      'X' gets the same binary representation as 'N'
 * Commented out the n4 information (nt_n4toa, nt_aton4).  This was
 * not used.
 * Modified nt_arevcomp so that 'P' and 'Q' are not treated as possible
 * codes and 'X' is treated as a any nucleotide ('N') as in GCG.  The
 * other alphabetic characters that are not codes are considered to be
 * NANT, as before.  The changes to nt_arevcomp are:
 *      revcomp of 'P' changed to NANT ('?') from 'Y'
 *      revcomp of 'Q' changed to NANT from 'R'
 *      revcomp of 'X' changed to 'N' from NANT
 *
 */
