/* (C) Copyright 1993-9, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* strutil.c: general string utilities */
/* Written by: Bill Alford */
/* Change log information is at the end of the file. */


/*	system headers not in global.h */
/*	blimps library headers not in global.h */
#include <global.h>

/*
 * Exported variables and data structures definitions
 */

char Buffer[LARGE_BUFF_LENGTH];

/*
 * Local variables and data structures
 */

/*
 * Function definitions
 */

/* 
 *  Boolean blank_line(s)
 *  char *eat_whitespace(s)
 *  void remove_trailing_whitespace(s)
 *  char *get_token(s)
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

Boolean blank_line(s)
     char *s;
{
  char *tmp;
  int length;

  tmp = eat_whitespace(s);
  length = strlen(tmp);

  if (length == 0) {
    return TRUE;
  }
  else {
    return FALSE;
  }
}


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

char *eat_whitespace(s)
     char *s;
{
/*  
  int white_length; 
*/

  while ((*s == ' ') ||
	 (*s == '\t') ||
	 (*s == '\n') ||
	 (*s == '\r')) {
    s++;
  }
  
  return (s);

/*
  white_length = strspn(s, WHITE_SPACE_CHARS);
  
  return (s+white_length);
*/
}


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

char *remove_trailing_whitespace(s)
     char *s;
{
  int length;

  length = strlen(s);

  /* if the string is empty, exit */
  if (length == 0) {
    return(NULL);
  }

  length--;				/* indexing starts at zero */

  while ((s[length] == ' ') ||
	 (s[length] == '\t') ||
	 (s[length] == '\n') ||
	 (s[length] == '\r')) {
    length--;
  }
  
  s[length+1] = '\0';

  return(s);
}

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

char *get_token(s)
     char *s;
{
#ifdef NO_STRTOK
  char c;
  char *return_string;
  static char *saved_string_pointer; /* pointer right after NULL char */



  if (s != NULL) {
    /* save the new string pointer */
    saved_string_pointer = s;
    remove_trailing_whitespace(saved_string_pointer);
  }

  /* if finished with this string last time, return NULL right away */
  if (saved_string_pointer == NULL) {
    return NULL;
  }

  /* eat up the beginning whitespace, if any */
  saved_string_pointer = eat_whitespace(saved_string_pointer);

  /* read ahead from saved pointer */
  return_string = saved_string_pointer;
  
  c = *saved_string_pointer;
  while ((c != ' ')  &&
    	 (c != '\t') &&
	 (c != '\n') &&
	 (c != '\r') &&
	 (c != '\0')) {
    saved_string_pointer++;
    c = *saved_string_pointer;
  }

  /* if saw '\0', finished the string */
  if (*saved_string_pointer == '\0') {
    saved_string_pointer = NULL;
    return return_string;
  }
  else {
    /* write in the null delimiting character */
    *saved_string_pointer = '\0';
    saved_string_pointer++;
    return return_string;
  }

#else
  return (char *)strtok(s, " \t\n\r\0");  
                                             /* typecast shouldn't be needed */
					     /* but compiler says it is a */
					     /* type mismatch of pointer and */
					     /* integer */
#endif /* NO_STRTOK */

}



/* Change log information follows. */
/* 
 * Revision 2.2002  1994/07/07  00:45:58  billa
 * Made remove_trailing_whitespace() return its argument.
 *
 * Revision 2.2001  1994/05/18  19:16:26  billa
 * Removed global.c and moved Buffer variable to strutil.[ch].
 *
 * Revision 2.1001  1994/04/21  01:03:18  billa
 * Creation.  Moved these from files.[ch].
 *
 */
