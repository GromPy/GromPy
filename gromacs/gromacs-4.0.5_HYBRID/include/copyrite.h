/*
 * $Id: copyrite.h,v 1.22 2008/02/28 14:07:18 spoel Exp $
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

#ifndef _copyrite_h
#define _copyrite_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef CPLUSPLUS
extern "C" {
#endif
  
#include <stdio.h>
  
/* Has to be a function, so we can get version number from autoconf */   
const char *GromacsVersion(void);
  
  
static const char *CopyrightText[] = {
  "Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.",
  "Copyright (c) 1991-2000, University of Groningen, The Netherlands.",
  "Copyright (c) 2001-2008, The GROMACS development team,",
  "check out http://www.gromacs.org for more information.\n"
};

static const char *GPLText[] = {
  "This program is free software; you can redistribute it and/or",
  "modify it under the terms of the GNU General Public License",
  "as published by the Free Software Foundation; either version 2",
  "of the License, or (at your option) any later version."
};


extern void pr_difftime(FILE *out,double dt);

void CopyRight(FILE *out,char *szProgram);
  
/* For both bromacs() and cool_quote() you have to provide a pointer to
 * a string of reasonable length (say 256) and the string length. This
 * is necessary to make the routines threadsafe and avoid allocating
 * a new string each time. The retstring pointer will be the return value.
 */
extern void bromacs(char *retstring, int retsize);
  
/* For cool_quote, the number of the quote used will be returned in cqnum 
 * if it is non-NULL. 
 */
extern void cool_quote(char *retstring, int retsize, int *cqnum);

extern int be_cool(void);
/* Return TRUE when the user is COOL, FALSE otherwise */

extern void thanx(FILE *fp);

enum { eCITEGMX, eCITEBATH, eCITESHAKE, eCITESETTLE, eCITESOR, 
       eCITEDISRE, eCITERF, eCITELINCS, eCITENR };

extern void please_cite(FILE *fp,char *key);
/* Print a message asking to cite something... */

#ifdef CPLUSPLUS
}
#endif

#endif	/* _copyright_h */
