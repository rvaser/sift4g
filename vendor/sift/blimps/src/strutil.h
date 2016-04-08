/* (C) Copyright 1993-7, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* strutil.h: general string utilities */
/* Written by: Bill Alford */
/* Change log information is at the end of the file. */

#ifndef STRUTIL_H_
#define STRUTIL_H_

/*
 * Exported variables and data structures
 */

#define EXTRA_LARGE_BUFF 1000
#define LARGE_BUFF_LENGTH 5000	/* The length of large buffers */
#define SMALL_BUFF_LENGTH 100	/* The length of small buffers and */
				/* predeclared strings */
extern char Buffer[LARGE_BUFF_LENGTH];


/*
 * Exported functions
 */

/* 
 * blank_line
 *  returns TRUE if the passed string is only whitespace
 *  Parameters:
 *    char* s: the string to check
 *  Return codes:
 *    TRUE if the line only has whitespace characters, FALSE otherwise.
 *  Error Codes:
 */

extern Boolean blank_line();


/*
 * eat_whitespace
 *   Eats the whitespace at the beginning of the string s1 and returns a 
 *   pointer to the first non whitespace character.  The whitespace characters
 *   are: space, tab, and newline.
 *   Parameters:
 *     char *s : input string
 *   Return Value: the return value is the pointer to the first non-whitespace
 *                 character
 *   Error codes:
 */

extern char *eat_whitespace();


/*
 * remove_trailing_whitespace
 *   Removes trailing whitespace by scanning back from the end of the string
 *   and replacing the leftmost trailing whitespace character with a NULL.
 *   Whitespace characters are: space, tab, newline and carriage return.
 *   NOTE: this is destructive
 *   Parameters:
 *     char *s : input string
 *   Return Value: returns its first argument
 *   Error codes:
 */

extern char *remove_trailing_whitespace();


/*
 * get_token
 *   get_token returns a pointer to a token delimited by whitespace.  The
 *   first call with a string returns the first token.  Subsequent calls (with
 *   the string as NULL) returns pointers to the following tokens.  A NULL
 *   is placed following the token.
 *   Note: get_token must be called with a non-null string pointer the first
 *         time or the behavior is ambiguous.
 *   Whitespace is: space, tab, newline, carriage return
 *   Parameters: 
 *     char *s: the string of tokens
 *   Return Value: the pointer to the token, NULL when there are no more tokens
 *   Error codes:
 *   Notes: get_token() wraps strtok() because the separators will always be 
 *          whitespace and strtok() is cryptic and easy to forget. (for me).
 */

extern char *get_token();




#endif  /* STRUTIL_H_ */

/* Change log information follows. */
/* 
 Changes since version 3.4:
12/23/00 Added EXTRA_LARGE_BUFF
 *
 */
