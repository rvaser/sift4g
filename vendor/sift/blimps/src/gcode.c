/* (C) Copyright 1993-9, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* gcode.c: Functions for the translation of nucleotides to amino acids */
/* Modified by: Bill Alford */
/* Change log information is at the end of the file. */

#include <global.h>
#include <residues.h>
#include <sequences.h>


/* -----------------------   initialize a genetic code   -------------------*/
void
init_gcode(gp, xltab, rcxltab)
    GeneticCodePtr  gp;
    register unsigned char   xltab[64], rcxltab[64];
{
    register char   *code;
    register int    i, j, k, tot;
    /* gctrans -- used to translate from the binary alphabet used in ntbet.h
    into a binary alphabet appropriate for the GeneticCode strings */
    static unsigned char    gctrans[] = { '\003', '\001', '\000', '\002' };

    code = gp->code;

    for (i=0; i<4; ++i) {
        for (j=0; j<4; ++j) {
            for (k=0; k<4; ++k, ++code) {
                tot = gctrans[i]*4*4 + gctrans[j]*4 + gctrans[k];
                xltab[tot] = aa_atob[(int)*code];
                tot = (3-gctrans[i])*4*4 + (3-gctrans[j])*4 + (3-gctrans[k]);
                rcxltab[tot] = aa_atob[(int)*code];
            }
        }
    }    

	return;
}



/*
codon2aa

Translate 3 binary-encoded nucleotides (n1, n2, and n3) into a binary
amino acid in the specified genetic code.
*/
unsigned char
codon2aa(gcode, n1, n2, n3)
    unsigned char   *gcode;
    unsigned    n1, n2, n3;
{
    if (n1 < 4 && n2 < 4 && n3 < 4)
        return gcode[n1*4*4 + n2*4 + n3];

    if (n1 >= NUCID_MAX || n2 >= NUCID_MAX || n3 >= NUCID_MAX)
        return (unsigned char)aa_atob[UNKNOWN_AA_CHR];

    {
        unsigned char   aa;
        unsigned char   b1, b2, b3;
        int i1, i2, i3;

        b1 = nt_bdegen[n1].list[0];
        b2 = nt_bdegen[n2].list[0];
        b3 = nt_bdegen[n3].list[0];
        aa = gcode[b1*4*4 + b2*4 + b3];

        for (i1=0; i1 < nt_bdegen[n1].ndegen; ++i1) {
            b1 = nt_bdegen[n1].list[i1]*4*4;
            for (i2=0; i2 < nt_bdegen[n2].ndegen; ++i2) {
                b2 = b1 + nt_bdegen[n2].list[i2]*4;
                for (i3=0; i3 < nt_bdegen[n3].ndegen; ++i3) {
                    b3 = nt_bdegen[n3].list[i3];
                    if (gcode[b2 + b3] != aa)
                        return (unsigned char)aa_atob[UNKNOWN_AA_CHR];
                }
            }    
        }
        return aa;
    }
/*NOTREACHED*/
}   /* end of codon2aa */


/*
aa2codon

Translate a single binary-encoded amino acid into 3 binary-encoded
nucleotides (n1, n2, and n3)
Note:  Could add translation tables to gcode.h
	The translations here are simplified to one degenerate codon per aa
	and therefore not quite correct for:
	Arg (R) = CGT, CGC, CGA, CGG; AGA, AGG = CGN or AGR
	Leu (L) = CTT, CTC, CTA, CTG; TTA, TTG = CTN or TTR
	Ser (S) = TCT, TCC, TCA, TCG; AGT, AGC = TCN or AGY
*/
void aa2codon(aa, n1, n2, n3)
Residue aa, *n1, *n2, *n3;
{
	switch (aa)
	{
	case 0: /*  -  => --- */
		{ *n1 = *n2 = *n3 = nt_atob['-']; }
	break;
	case 1: /*  A  => GCN */
		{ *n1 = nt_atob['G']; *n2 = nt_atob['C']; *n3 = nt_atob['N']; }
	break;
	case 2: /*  R  => MGN */
		{ *n1 = nt_atob['M']; *n2 = nt_atob['G']; *n3 = nt_atob['N']; }
	break;
	case 3: /*  N  => AAY */
		{ *n1 = nt_atob['A']; *n2 = nt_atob['A']; *n3 = nt_atob['Y']; }
	break;
	case 4: /*  D  => GAY */
		{ *n1 = nt_atob['G']; *n2 = nt_atob['A']; *n3 = nt_atob['Y']; }
	break;
	case 5: /*  C  => UGY */
		{ *n1 = nt_atob['U']; *n2 = nt_atob['G']; *n3 = nt_atob['Y']; }
	break;
	case 6: /*  Q  => CAR */
		{ *n1 = nt_atob['C']; *n2 = nt_atob['A']; *n3 = nt_atob['R']; }
	break;
	case 7: /*  E  => GAR */
		{ *n1 = nt_atob['G']; *n2 = nt_atob['A']; *n3 = nt_atob['R']; }
	break;
	case 8: /*  G  => GGN */
		{ *n1 = nt_atob['G']; *n2 = nt_atob['G']; *n3 = nt_atob['N']; }
	break;
	case 9: /*  H  => CAY */
		{ *n1 = nt_atob['C']; *n2 = nt_atob['A']; *n3 = nt_atob['Y']; }
	break;
	case 10: /*  I  => ATH */
		{ *n1 = nt_atob['A']; *n2 = nt_atob['T']; *n3 = nt_atob['H']; }
	break;
	case 11: /*  L  => YTN */
		{ *n1 = nt_atob['Y']; *n2 = nt_atob['T']; *n3 = nt_atob['N']; }
	break;
	case 12: /*  K  => AAR */
		{ *n1 = nt_atob['A']; *n2 = nt_atob['A']; *n3 = nt_atob['R']; }
	break;
	case 13: /*  M  => ATG */
		{ *n1 = nt_atob['A']; *n2 = nt_atob['T']; *n3 = nt_atob['G']; }
	break;
	case 14: /*  F  => TTY */
		{ *n1 = nt_atob['T']; *n2 = nt_atob['T']; *n3 = nt_atob['Y']; }
	break;
	case 15: /*  P  => CCN */
		{ *n1 = nt_atob['C']; *n2 = nt_atob['C']; *n3 = nt_atob['N']; }
	break;
	case 16: /*  S  => WSN */
		{ *n1 = nt_atob['W']; *n2 = nt_atob['S']; *n3 = nt_atob['N']; }
	break;
	case 17: /*  T  => ACN */
		{ *n1 = nt_atob['A']; *n2 = nt_atob['C']; *n3 = nt_atob['N']; }
	break;
	case 18: /*  W  => TGG */
		{ *n1 = nt_atob['T']; *n2 = nt_atob['G']; *n3 = nt_atob['G']; }
	break;
	case 19: /*  Y  => TAY */
		{ *n1 = nt_atob['T']; *n2 = nt_atob['A']; *n3 = nt_atob['Y']; }
	break;
	case 20: /*  V  => GTN */
		{ *n1 = nt_atob['G']; *n2 = nt_atob['T']; *n3 = nt_atob['N']; }
	break;
	case 21: /*  B (D or N)  => RAY */
		{ *n1 = nt_atob['R']; *n2 = nt_atob['A']; *n3 = nt_atob['Y']; }
	break;
	case 22: /*  Z (E or Q)  => SAR */
		{ *n1 = nt_atob['S']; *n2 = nt_atob['A']; *n3 = nt_atob['R']; }
	break;
	case 23: /*  X => NNN  */
		{ *n1 = nt_atob['N']; *n2 = nt_atob['N']; *n3 = nt_atob['N']; }
	break;
	case 24: /*  * => TRR  (stop) */
		{ *n1 = nt_atob['T']; *n2 = nt_atob['R']; *n3 = nt_atob['R']; }
	break;
	default: /*  Unknown */
		{ *n1 = nt_atob['N']; *n2 = nt_atob['N']; *n3 = nt_atob['N']; }
	break;
	}
}   /*  end of aa2codon  */


/* Change log information follows. */
/* Changes since Blimps Version 3.1:
 * 1/19/1997 Added aa2codon() routine.
 * 
 */

