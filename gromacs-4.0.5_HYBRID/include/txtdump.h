/*
 * $Id: txtdump.h,v 1.26 2008/09/19 08:09:33 lindahl Exp $
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

#ifndef _txtdump_h
#define _txtdump_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <stdio.h>
#include "typedefs.h"
#include "tpxio.h"

#define LINE_WIDTH	80
#define RMARGIN		10
#define	USE_WIDTH	((LINE_WIDTH)-(RMARGIN))
#define INDENT		3

extern char *atomname(t_atoms *a,int i);
/* Return pointer to a buffer which holds the atomname in the
 * form resname resnr atomname. Pointer can be freed afterwards.
 */
int pr_indent(FILE *fp,int n);
int available(FILE *fp,void *p,int indent,const char *title);
int pr_title(FILE *fp,int indent,const char *title);
int pr_title_n(FILE *fp,int indent,const char *title,int n);
int pr_title_nxn(FILE *fp,int indent,const char *title,int n1,int n2);
void pr_ivec(FILE *fp,int indent,const char *title,int vec[],int n, bool bShowNumbers);
void pr_ivecs(FILE *fp,int indent,const char *title,ivec vec[],int n, bool bShowNumbers);
void pr_rvec(FILE *fp,int indent,const char *title,real vec[],int n, bool bShowNumbers);
void pr_dvec(FILE *fp,int indent,const char *title,double vec[],int n, bool bShowNumbers);
void pr_rvecs(FILE *fp,int indent,const char *title,rvec vec[],int n);
void pr_rvecs_len(FILE *fp,int indent,const char *title,rvec vec[],int n);
void pr_reals(FILE *fp,int indent,const char *title,real vec[],int n);
void pr_doubles(FILE *fp,int indent,const char *title,double *vec,int n);
void pr_block(FILE *fp,int indent,const char *title,t_block *block,bool bShowNumbers);
void pr_blocka(FILE *fp,int indent,const char *title,t_blocka *block,bool bShowNumbers);
void pr_ilist(FILE *fp,int indent,const char *title,
	       t_functype *functype,t_ilist *ilist, bool bShowNumbers);
void pr_iparams(FILE *fp,t_functype ftype,t_iparams *iparams);
void pr_idef(FILE *fp,int indent,const char *title,t_idef *idef, bool bShowNumbers);
void pr_inputrec(FILE *fp,int indent,const char *title,t_inputrec *ir,
		 bool bMDPformat);
void pr_atoms(FILE *fp,int indent,const char *title,t_atoms *atoms, 
	      bool bShownumbers);
void pr_atomtypes(FILE *fp,int indent,const char *title,
		  t_atomtypes *atomtypes,bool bShowNumbers);
void pr_mtop(FILE *fp,int indent,const char *title,gmx_mtop_t *mtop,
	     bool bShowNumbers);
void pr_top(FILE *fp,int indent,const char *title,t_topology *top, bool bShowNumbers);
/*
 * This routine prints out a (human) readable representation of 
 * the topology to the file fp. Ident specifies the number of 
 * spaces the text should be indented. Title is used to print a 
 * header text.
 */
void pr_header(FILE *fp,int indent,const char *title,t_tpxheader *sh);
/*
 * This routine prints out a (human) readable representation of
 * a header to the file fp. Ident specifies the number of spaces
 * the text should be indented. Title is used to print a header text.
 */

void pr_commrec(FILE *fp,int indent,t_commrec *cr);

#endif	/* _txtdump_h */