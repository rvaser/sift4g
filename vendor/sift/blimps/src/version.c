/* (C) Copyright 1993-9, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* version.c: Keeps track of the system version level. */
/* Written by: Bill Alford */
/* Change log information is at the end of the file. */

/*	system headers not in global.h  */
/*	blimps library headers; global.h must be first  */
#include <global.h>
#include <version.h>

/****
The note about versions is out of date.

Versions are now as follows:
major_release.minor_release.updates

Major releases, minor releases, and updates are described as below.

*****/
/**************************************************************************
A NOTE ABOUT VERSIONS:

The way I am computing system versions is this:
  Take the RCS version number of version.c (the development version)
  The whole number part of the system version is the same as the whole
    number part of the version.c version.
  The fractional part of the system version is the integer value of
    the fractional part of version.c modulus 1000.
  The letter part of the version is the 100's place of the fractional
    part (taken as an integer).  eg. A = 0, B = 1, C = 2, ...

  version.c version: a.xxxxyzz  --> system version: a.xxxx ('A' + y)

With this method of versions every system version update can be done
by checking out and checking in everything at that level.  This makes
the system version level consistient and allows 99 incremental
updates.  Also to retrieve an earlier system version, all that is
needed to be done is to check out all the files of that development
version.

A caveat: ALWAYS check in a new version of version.c everytime another
          system file is checked in (use ci -f).  And ALWAYS say what
	  files and their version number in the log of the version.c
	  checkin.  If this is always done, we should be able to get
	  to any of the incremental versions between the system
	  versions.

Another caveat: Only the first system version of the development
                version is released.  This way to get any particular
		version, all the files only have to be checked out at
		one development version level.  Also the minor changes
		in the incremental versions will not be in different
		copies of the same system version. eg.:
		  development    system
		  1.0          = 1.0 A
		  1.1-1.99     = incremental increase from 1.0 A to 1.0 B
		  1.100        = 1.0 B
		  1.101-1.199  = incremental increase from 1.0 B to 1.0 C
		  ...
		  1.1000       = 1.1 A
		  ...

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
 * Local variables and data structures
 */

static char ProgramString[SMALL_BUFF_LENGTH];
static char VersionString[SMALL_BUFF_LENGTH];
static char VersionInfo[SMALL_BUFF_LENGTH];
static char DateString[SMALL_BUFF_LENGTH];
static char CopyrightString[SMALL_BUFF_LENGTH];

static char TitleString[SMALL_BUFF_LENGTH] = { '\0' };
				/* make sure it is initialized to null */

/*
 * Function definitions
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

void version_strings(program_string, version_string, version_info,
		     date_string, copyright_string)
     char *program_string, *version_string, *version_info,
       *date_string, *copyright_string;
{
  strncpy(ProgramString,   program_string,   SMALL_BUFF_LENGTH);
  strncpy(VersionString,   version_string,   SMALL_BUFF_LENGTH);
  strncpy(VersionInfo,     version_info,     SMALL_BUFF_LENGTH);
  strncpy(DateString,      date_string,      SMALL_BUFF_LENGTH);
  strncpy(CopyrightString, copyright_string, SMALL_BUFF_LENGTH);
}


/*
 * print_version
 *   Prints the version of the current release.
 *   Parameters:
 *     FILE *fp: the file to print the verson info, typically stdout and the
 *               output file
 *   Error codes: none
 */

void print_version(fp)
     FILE *fp;
{

  if (TitleString[0] == '\0') {
    /* if the development version is different from the system version */
    /* then make the development_version string */
    /* if something in ones or tens place, this is a develoment version */
    if (VersionInfo[0] != '\0') {
      sprintf(TitleString,
	      "%s   Version %s (%s)  %s%s\n",
	      ProgramString, VersionString, VersionInfo,
	      DateString, VERSION_FROM_MAKE);
    }
    else {
      sprintf(TitleString,
	      "%s   Version %s  %s%s\n",
	      ProgramString, VersionString, DateString, VERSION_FROM_MAKE);
    }
  }

  fprintf(fp, "%s\n", TitleString);
  fprintf(fp, "%s\n", CopyrightString);

}

/* Change log information follows. */
/*
 *
 */
