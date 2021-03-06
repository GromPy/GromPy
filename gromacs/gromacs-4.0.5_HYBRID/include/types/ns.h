/*
 * $Id: ns.h,v 1.1.2.1 2008/10/07 22:49:40 lindahl Exp $
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "nsgrid.h"

enum { eNL_VDWQQ, eNL_VDW, eNL_QQ, 
       eNL_VDWQQ_FREE, eNL_VDW_FREE, eNL_QQ_FREE, 
       eNL_VDWQQ_WATER, eNL_QQ_WATER, 
       eNL_VDWQQ_WATERWATER, eNL_QQ_WATERWATER, 
       eNL_NR };

#define MAX_CG 1024

typedef struct {
  int     ncg;
  int     nj;
  atom_id jcg[MAX_CG];
} t_ns_buf;

typedef unsigned long t_excl;

typedef struct {
  atom_id  *simple_aaj;
  t_grid   *grid;
  t_excl   *bexcl;
  bool     *bHaveVdW;
  t_ns_buf **ns_buf;
  bool     *bExcludeAlleg;
  int      nra_alloc;
  int      cg_alloc;
  atom_id  **nl_sr;
  int      *nsr;
  atom_id  **nl_lr_ljc;
  atom_id  **nl_lr_one;
  int      *nlr_ljc;
  int      *nlr_one;
} gmx_ns_t;
