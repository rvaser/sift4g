/* (C) Copyright 1993, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* options.c: */
/* Written by: Bill Alford */
/* Change log information is at the end of the file. */

#define EXTERN

/*	system headers not in global.h */
/*	blimps library headers */
#include <global.h>
#include <options.h>	/* requires global.h ! */
#include <strings.h>

char **OptionsARGV = NULL;
int OptionsARGC = 0;

/* Eliminate compiler warning */
char *strdup(const char *s);

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

#define KEYCHAR ':'

Boolean 
get_option_args(key, argc_ptr, argv_ptr)
     char *key;
     int *argc_ptr;
     char ***argv_ptr;
{
  int x, start, stop;
  int keylen;

  keylen = strlen(key);
  x = 0;
  while ((x < OptionsARGC) &&
	 ! ((strncasecmp(OptionsARGV[x], key, keylen) == 0) &&
	    (OptionsARGV[x][keylen] == KEYCHAR))) {
    x++;
  }

  if (x >= OptionsARGC) {
    *argc_ptr = -1;
    *argv_ptr = NULL;
    return FALSE;
  } 

  start = x;
	 
  while ((x < OptionsARGC) &&
	 ! ((strncasecmp(&(OptionsARGV[x][1]), key, keylen) == 0) &&
	    (OptionsARGV[x][0] == KEYCHAR))) {
    x++;
  }

  stop = x;

  *argc_ptr = stop - start - 1;

/* debug stuff
printf("Start= %d, stop = %d\n", start, stop);
printf("OptionsARGV[start]= %s\n", OptionsARGV[start]);
*/

  *argv_ptr = &(OptionsARGV[start+1]);

  return TRUE;

}



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

void
insert_into_options(string)
     char *string;
{
  int x;
  char ** argv;
  int argc;
  int max_arg;
  char *tmp;
  char **tmp_argv;

  int cgt; /* whether a call to get token has been done with the string yet */
  cgt = 0;

  /* allocate a temporary argv */
  CheckMem(
	   argv = (char **) malloc(sizeof(char *) * 5)	/* allocate 5 ptrs */
	   );

  max_arg=5;
  argc = 0;

  /* fill it */
  while ( (tmp = ((cgt) ?
		  (get_token(NULL)) :
		  ((cgt=1), (get_token(string))) )) != NULL) {
    /* increase if needed */
    if (argc >= max_arg) {
      CheckMem(
	       tmp_argv = (char **) realloc(argv, sizeof(char *) * (max_arg+5))
	       );
      argv = tmp_argv;
      max_arg = max_arg+5;
    }

    /* allocate space for the duplicate string and put it in the array */
    /* strdup is a non-ansi function, therefore it needs it a prototype,
     * otherwise a compiler warning is generated. */
    CheckMem(
	     argv[argc] = strdup(tmp)
	     );
    argc++;
  }

  /* do NOT need check if the OptionsARG[VC] are already used since */
  /* realloc acts like malloc if the pointer was NULL */
  /* WRONG!!! Some systems do not have realloc like this. :P */
  
  /* (re)allocate enough space for the new argv */
  if (OptionsARGV == NULL) {
    CheckMem(
	     tmp_argv = (char **) malloc(sizeof(char *) * (OptionsARGC+argc))
	     );
  }
  else {
    CheckMem(
	     tmp_argv = (char **) realloc(OptionsARGV, 
					  sizeof(char *) * (OptionsARGC+argc))
	     );
  }
  OptionsARGV = tmp_argv;
  
  /* add the new args to OptionsARGV */
  for (x=0; x<argc; x++) {
    OptionsARGV[OptionsARGC+x] = argv[x];
  }
  OptionsARGC += argc;

  free(argv);

  /* DEBUG: checking if valid */
  /*
  {
    int argc;
    char **argv;
    
    if (get_option_args("BLT", & argc, & argv)) {
      for (x=0; x<argc; x++) {
	printf("\"%s\" ", argv[x]);
      }
      printf("\nargc= %d\n", argc);
    }
    else {
      printf("request failed\n");
    }
  }
  */
}




/* Change log information follows. 
* $Log: options.c,v $
* Revision 1.3  2011-09-28 22:10:55  gsims
* *** empty log message ***
*
* Revision 1.2  2011-06-03 18:33:01  gsims
* *** empty log message ***
*
* Revision 1.1  2011-05-24 16:32:28  gsims
* *** empty log message ***
*
* Revision 1.1.1.1  2011-02-22 22:16:00  gsims
* Initial Import of sift 4.0.4
*
 * Revision 2.2010  1995/07/28  23:47:14  billa
 * Added new convert method and pattern matching and minor updates merged.
 * */
