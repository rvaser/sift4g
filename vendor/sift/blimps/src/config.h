/* (C) Copyright 1993-7, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* config.h: reads the configuration file. */
/* Written by: Bill Alford */
/* Change log information is at the end of the file. */

#ifndef CONFIG_H_
#define CONFIG_H_

/*
 * Exported variables and data structures
 */
 
/*
 * Exported functions
 */
 
/*
 * read_config_file
 *   reads the configuration file given by the filename
 *   Parameters:
 *     char * filename  : the configuration filename
 *   Error codes:
 */

extern int read_config_file ();


#endif /*  CONFIG_H_ */

/* Change log information follows. */
/* 
 * System version 2.1 A [2.1000]
 * Added reading and scoring of precomputed site specific scoring matricies.
 * Reconfigured files.[ch] into files.[ch] and config.[ch].
 *
 * Revision 2.1  1994/02/23  20:42:46  billa
 * Creation.  Moved the configuration file stuff from files.[ch] to
 * config.[ch].
 *
 */
