/* (C) Copyright 1993-9, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* global.h: basic global defines and variable declarations used in *
 *           all blimps modules                                     */
/* Put system headers here */
/* Written by: Bill Alford */
/* Change log information is at the end of the file. */

#ifndef GLOBAL_H_
#define GLOBAL_H_

/*	system headers used by all modules */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <ctype.h>
#include <malloc.h>

/*
 * Exported variables and data structures
 */


/*
 * Boolean type and values
 */
typedef int Boolean;

#ifndef TRUE
#define TRUE  1
#endif
#ifndef FALSE
#define FALSE 0
#endif



/*
 * min and max and round
 */
#define max(A, B) ((A) > (B) ? (A) : (B))
#define min(A, B) ((A) < (B) ? (A) : (B))
#define round(x) ((x >= 0.0) ? (int) (x+0.5) : (int) (x-0.5))


/*	blimps library headers used by all modules */
#include "license.h"
#include "strutil.h"
#include "errors.h"
#include "memory.h"

#endif /*  GLOBAL_H_ */

/* Change log information follows. */
/* Changes since version 3.0.0:
12/ 6/99 Put common system & blimps headers here.
*/
