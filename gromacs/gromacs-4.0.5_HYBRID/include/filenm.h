/*
 * $Id: filenm.h,v 1.19 2008/09/19 08:09:33 lindahl Exp $
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _filenm_h
#define _filenm_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#ifdef CPLUSPLUS
extern "C" {
#endif
  
#include "futil.h"

void set_default_file_name(const char *name);
/* Set the default file name for all file types to name */

extern const char *ftp2ext(int ftp);
/* Return extension for filetype */

  extern const char *ftp2ext_generic(int ftp);
/* Return extension for filetype, and a generic name for generic types (e.g. trx)*/

extern const char *ftp2desc(int ftp);
/* Return description for file type */

extern const char *ftp2defnm(int ftp);
/* Return default file name for file type */

extern const char *ftp2ftype(int ftp);
/* Return Binary or ASCII depending on file type */

extern void pr_def(FILE *fp,int ftp);
/* Print definitions for filename ftp */

extern void pr_defs(FILE *fp);
/* Print definitions for all filename */

extern void pr_fns(FILE *fp,int nf,t_filenm tfn[]);
/* Print nf file names and types */

extern void pr_fopts(FILE *fp,int nf,t_filenm tfn[], int shell);
/* prints file options in tcsh 'complete' format */

extern void parse_file_args(int *argc,char *argv[],int nf,t_filenm fnm[],
			    bool bKeep);
/* Parse command line for file names. When bKeep is set args are 
 * not removed from argv. */

extern char *opt2fn(char *opt,int nfile,t_filenm fnm[]);
/* Return the filename belonging to cmd-line option opt, or NULL when 
 * no such option. */

extern int opt2fns(char **fns[], char *opt,int nfile,t_filenm fnm[]);
/* Return the filenames belonging to cmd-line option opt, or NULL when 
 * no such option. */

#define opt2FILE(opt,nfile,fnm,mode) ffopen(opt2fn(opt,nfile,fnm),mode)
/* Return a file pointer from the filename (see above) */

extern int fn2ftp(const char *fn);
/* Return the filetype corrsponding to filename */

extern char *ftp2fn(int ftp,int nfile,t_filenm fnm[]);
/* Return the first file name with type ftp, or NULL when none found. */

extern int ftp2fns(char **fns[], int ftp,int nfile,t_filenm fnm[]);
/* Return the number of files for the first option with type ftp
   and the files in **fns[] (will be allocated), or NULL when none found. */
  
extern char *ftp2filter(int ftp);
/* Return a file extension filter for file type */

#define ftp2FILE(ftp,nfile,fnm,mode) ffopen(ftp2fn(ftp,nfile,fnm),mode)
/* Return a file pointer from the filename (see above) */

extern bool ftp2bSet(int ftp,int nfile,t_filenm fnm[]);
/* Return TRUE when this file type has been found on the cmd-line */

extern bool opt2bSet(char *opt,int nfile,t_filenm fnm[]);
/* Return TRUE when this option has been found on the cmd-line */

extern char *opt2fn_null(char *opt,int nfile,t_filenm fnm[]);
/* Return the filenm belonging top cmd-line option opt, or NULL when 
 * no such option. 
 * Also return NULL when opt is optional and option is not set. 
 */

extern char *ftp2fn_null(int ftp,int nfile,t_filenm fnm[]);
/* Return the first file name with type ftp, or NULL when none found.
 * Also return NULL when ftp is optional and option is not set.
 */

extern bool is_optional(t_filenm *fnm);
/* Return whether or not this filenm is optional */

extern bool is_output(t_filenm *fnm);
/* Return whether or not this filenm is output */

extern bool is_set(t_filenm *fnm);
/* Return whether or not this filenm is set */

/* When we do checkpointing, this routine is called to check for previous
 * output files and append a 'partNNNN' suffix before the (output) file extensions.
 */
int
add_suffix_to_output_names(t_filenm *fnm, int nfile, char *suffix);
	
#ifdef CPLUSPLUS
}
#endif

#endif	/* _filenm_h */
