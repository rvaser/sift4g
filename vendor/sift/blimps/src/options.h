/* (C) Copyright 1993, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* options.h: */
/* Written by: Bill Alford */
/* Change log information is at the end of the file. */

#ifndef OPTIONS_H_
#define OPTIONS_H_

extern char **OptionsARGV;	/* Global options/arguments for */
extern int OptionsARGC;		/* sub-algorithms*/


/* 
 * get_option_args
 *  This retrieves the desired option arguments from the OptionsARGV
 *  array.  This is the function to use to retrieve options for any of
 *  the sub-algorithms (algorithms that are optional and not used
 *  generally).
 *  The options/arguments are accessed from the returned pointer.
 *  NOTE: there is no argument quoting yet.  All arguments are
 *    whitespace delimited.
 *  To identify which algorithm a set of arguments belongs to the
 *  arguments need to be prefixed and postfixed with a key.  For
 *  example (the colon is needed):
 *    OP Alt: arg0 arg1 arg2 arg3 :Alt
 *
 *  NOTE: The key is case INsensitive.  Also, if there is no
 *    terminating key the rest of the arguments are included in the
 *    returned argc.
 *  In the above example the returned value for argc would be 4 and
 *  the argv vector would start at arg0.
 *  WARNING: all of the following arguments after the key region are
 *    STILL accessible.  Only the address of the first arg is set.  No
 *    memory is duplicated.
 *  IMPORTANT NOTE: This does NOT follow the standard argc/argv
 *    conventions in main().  The argv[0] is the first ARGUMENT.
 *  Parameters:
 *    char *key:        the key (case INsensitive)
 *    int *argc_ptr:    a pointer to where to place argc.
 *    char ***argv_ptr: a pointer to where to place the argv address.
 *                      (remember argv is accessed as an array of char
 *                      ptrs.  An array of char ptrs is of type char **.  
 *                      so this is just a pointer to the place to put
 *                      the "array")
 *  Return Codes: returns TRUE if it can find the key region, FALSE
 *                otherwise
 *  Error Codes:  See above and argc is set to -1 and argv is set to NULL
 *
 *  Calling example:
 *    int argc;
 *    char **argv;  -- this can NOT be an array
 *    int result;
 *
 *    result = get_option_args("mykey", & argc, & argv); -- pass the addresses!
 *    if (result == TRUE) {
 *      printf("Got them\n");
 *    }
 *    else {
 *      printf("Could not find the args\n");
 *      -- error handling stuff here -- 
 *    }
 */

extern Boolean get_option_args();


/* 
 * insert_into_options
 *   Inserts the tokens in string into the options array for later
 *   retreival. 
 *   NOTE, get_token() is the function called to parse the tokens.  If
 *     a parse is in progress just pass in NULL, this will get passed
 *     to get_token() the first time, retreiving the next token.
 *   Parameters:
 *     char  *string: the string of tokens whitespace delimited
 *   Return codes: none
 *   Error codes: none
 */

extern void insert_into_options();


#endif /*  OPTIONS_H_ */

/* Change log information follows. 
 * $Log: options.h,v $
 * Revision 1.2  2011-06-03 18:33:01  gsims
 * *** empty log message ***
 *
 * Revision 1.1  2011-05-24 16:32:31  gsims
 * *** empty log message ***
 *
 * Revision 1.1.1.1  2011-02-22 22:16:02  gsims
 * Initial Import of sift 4.0.4
 *
 * Revision 2.2010  1995/07/28  23:47:14  billa
 * Added new convert method and pattern matching and minor updates merged.
 * */
