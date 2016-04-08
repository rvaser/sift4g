/* (C) Copyright 1993-2006, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* version.h: Header file to version.c. keeps track of the system */
/*            version level. */
/* Written by: Bill Alford */
/* Change log information is at the end of the file. */

#ifndef VERSION_H_
#define VERSION_H_

/*  Remember to put the DEVELOPMENT VERSION string in the VERSION define */
#define VERSION "3.8"
#define VERSION_DATE "2006/04"

#ifndef VERSION_FROM_MAKE
#define VERSION_FROM_MAKE ""
#endif


/****
Versions are now as follows:
major_release.minor_release.updates

Major releases, minor releases, and updates are described as below.

*****/
/**************************************************************************
A NOTE ABOUT VERSIONS:

The plan for the difference between versions:
  Major release:  
    A major release would be something like a reconstruction of the
    interface.  New file formats or configuration keys that are
    necessary in order to run.  Also changes to the program structure
    should be limited to a major release.  An example would be
    changing the way the frequency file is laid out.

    ex. system: 1.x c    --> 2.0 A
        dev.    1.xxxyzz --> 2.0

  Minor release:  
    A minor release would be something like an important enhancement
    to a file format or a configuration that could change the
    performance (but not the output) of the program, enhancements to
    the program that do not change the data, but make the program
    better.  An example would be memory management code or new
    names to the configuration file or reading of a matrix database/file.

    ex. system: 1.2 c    --> 1.3 A
        dev.    1.2yzz   --> 1.3000

  Updates:     
    An update would be something that would not get in the way of
    earlier file formats and that could be used by earlier versions of
    this minor release.  Nothing that changes the way the program is
    called.  Examples of updates would be new conversion and scoring
    routines added or a new configuration key that will change the way
    data is output (eg. outputs a matrix) or the ability to read a new
    sequence database format.  Note that reading a new block format
    would be a minor release because there is only one block database
    and it is key to the way the program behaves.  Also note that
    something like being able to read a new sequence database format
    is more suited to an update since there are many different formats
    and external to the program they can be converted to one that
    blimps can already read.  Updates are also for bug fixes.
   
    ex. system: 1.2 B    --> 1.2 C
        dev.    1.22yy   --> 1.2300

**************************************************************************/



/*
 * Exported variables and data structures
 */

/*
 * Exported functions
 */


/*
 * version_strings
 *   Enters the strings to be used by the version output function 
 *   (print_version).
 *   Parameters:
 *     char *program_string:   the name of the program
 *     char *version_string:   the RCS revision number 
 *     char *version_info:     additional version info
 *     char *date_string:      the RCS date
 *     char *copyright_string: the copyright info
 *   Error codes: none
 */

extern void version_strings();


/*
 * print_version
 *   Prints the version of the current release.  
 *   NOTE: version.c must be checked in each time a change is made so that
 *         the version info will be correct
 *   Parameters:
 *     FILE *fp: the file to print the verson info, typically stdout and the
 *               output file
 *   Error codes: none
 */

extern void print_version();


#endif /*  VERSION_H_ */

/* Change log information follows. */
