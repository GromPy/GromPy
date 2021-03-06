/*
 * $Id: typedefs.c,v 1.51.2.1 2008/11/06 10:03:32 hess Exp $
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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
/* This file is completely threadsafe - keep it that way! */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "smalloc.h"
#include "symtab.h"
#include "vec.h"
#include "pbc.h"
#include <string.h>

static bool bOverAllocDD=FALSE;

void set_over_alloc_dd(bool set)
{
  bOverAllocDD = set;
}

int over_alloc_dd(int n)
{
  if (bOverAllocDD)
    return OVER_ALLOC_FAC*n + 100;
  else
    return n;
}

void init_block(t_block *block)
{
  int i;

  block->nr           = 0;
  block->nalloc_index = 1;
  snew(block->index,block->nalloc_index);
  block->index[0]     = 0;
}

void init_blocka(t_blocka *block)
{
  int i;

  block->nr           = 0;
  block->nra          = 0;
  block->nalloc_index = 1;
  snew(block->index,block->nalloc_index);
  block->index[0]     = 0;
  block->nalloc_a     = 0;
  block->a            = NULL;
}

void init_atom(t_atoms *at)
{
  int i;

  at->nr       = 0;
  at->nres     = 0;
  at->atom     = NULL;
  at->resname  = NULL;
  at->atomname = NULL;
  at->atomtype = NULL;
  at->atomtypeB= NULL;
  at->pdbinfo  = NULL;
}

void init_atomtypes(t_atomtypes *at)
{
  at->nr = 0;
  at->radius = NULL;
  at->vol = NULL;
  at->atomnumber = NULL;
}

void init_groups(gmx_groups_t *groups)
{
  int g;

  groups->ngrpname = 0;
  groups->grpname  = NULL;
  for(g=0; (g<egcNR); g++) {
    groups->grps[g].nm_ind = NULL;
    groups->ngrpnr[g] = 0;
    groups->grpnr[g]  = NULL;
  }

}

void init_mtop(gmx_mtop_t *mtop)
{
  mtop->name = NULL;
  mtop->nmoltype = 0;
  mtop->moltype = NULL;
  mtop->nmolblock = 0;
  mtop->molblock = NULL;
  init_groups(&mtop->groups);
  init_block(&mtop->mols);
  open_symtab(&mtop->symtab);
}

void init_top (t_topology *top)
{
  int i;
  
  top->name = NULL;
  init_atom (&(top->atoms));
  init_atomtypes(&(top->atomtypes));
  init_block(&top->cgs);
  init_block(&top->mols);
  init_blocka(&top->excls);
  open_symtab(&top->symtab);
}

void init_inputrec(t_inputrec *ir)
{
  memset(ir,0,(size_t)sizeof(*ir));
}

void stupid_fill_block(t_block *grp,int natom,bool bOneIndexGroup)
{
  int i;

  if (bOneIndexGroup) {
    grp->nalloc_index = 2;
    snew(grp->index,grp->nalloc_index);
    grp->index[0]=0;
    grp->index[1]=natom;
    grp->nr=1;
  }
  else {
    grp->nalloc_index = natom+1;
    snew(grp->index,grp->nalloc_index);
    snew(grp->index,natom+1);
    for(i=0; (i<=natom); i++)
      grp->index[i]=i;
    grp->nr=natom;
  }
}

void stupid_fill_blocka(t_blocka *grp,int natom)
{
  int i;

  grp->nalloc_a = natom;
  snew(grp->a,grp->nalloc_a);
  for(i=0; (i<natom); i++)
    grp->a[i]=i;
  grp->nra=natom;
  
  grp->nalloc_index = natom + 1;
  snew(grp->index,grp->nalloc_index);
  for(i=0; (i<=natom); i++)
    grp->index[i]=i;
  grp->nr=natom;
}

void copy_blocka(const t_blocka *src,t_blocka *dest)
{
  int i;

  dest->nr = src->nr;
  dest->nalloc_index = dest->nr + 1;
  snew(dest->index,dest->nalloc_index);
  for(i=0; i<dest->nr+1; i++) {
    dest->index[i] = src->index[i];
  }
  dest->nra = src->nra;
  dest->nalloc_a = dest->nra + 1;
  snew(dest->a,dest->nalloc_a);
  for(i=0; i<dest->nra+1; i++) {
    dest->a[i] = src->a[i];
  }
}

void done_block(t_block *block)
{
  block->nr    = 0;
  sfree(block->index);
  block->nalloc_index = 0;
}

void done_blocka(t_blocka *block)
{
  block->nr    = 0;
  block->nra   = 0;
  sfree(block->index);
  if (block->a)
    sfree(block->a);
  block->nalloc_index = 0;
  block->nalloc_a = 0;
}

void done_atom (t_atoms *at)
{
  at->nr       = 0;
  at->nres     = 0;
  sfree(at->atom);
  sfree(at->resname);
  sfree(at->atomname);
}

void done_moltype(gmx_moltype_t *molt)
{
  int f;
  
  done_atom(&molt->atoms);
  done_block(&molt->cgs);
  done_blocka(&molt->excls);

  for(f=0; f<F_NRE; f++) {
    sfree(molt->ilist[f].iatoms);
    molt->ilist[f].nalloc = 0;
  }
}

void done_molblock(gmx_molblock_t *molb)
{
  if (molb->nposres_xA > 0) {
    molb->nposres_xA = 0;
    free(molb->posres_xA);
  }
  if (molb->nposres_xB > 0) {
    molb->nposres_xB = 0;
    free(molb->posres_xB);
  }
}

void done_mtop(gmx_mtop_t *mtop,bool bDoneSymtab)
{
  int i;

  if (bDoneSymtab) {
    done_symtab(&mtop->symtab);
  }

  sfree(mtop->ffparams.functype);
  sfree(mtop->ffparams.iparams);

  for(i=0; i<mtop->nmoltype; i++) {
    done_moltype(&mtop->moltype[i]);
  }
  sfree(mtop->moltype);
  for(i=0; i<mtop->nmolblock; i++) {
    done_molblock(&mtop->molblock[i]);
  }
  sfree(mtop->molblock);
  done_block(&mtop->mols);
}

void done_top(t_topology *top)
{
  int i;
  
  done_atom (&(top->atoms));
  done_symtab(&(top->symtab));
  done_block(&(top->cgs));
  done_block(&(top->mols));
  done_blocka(&(top->excls));
}

static void done_pullgrp(t_pullgrp *pgrp)
{
  sfree(pgrp->ind);
  sfree(pgrp->ind_loc);
  sfree(pgrp->weight);
  sfree(pgrp->weight_loc);
}

static void done_pull(t_pull *pull)
{
  int i;

  for(i=0; i<pull->ngrp+1; i++) {
    done_pullgrp(pull->grp);
    done_pullgrp(pull->dyna);
  }
}

void done_inputrec(t_inputrec *ir)
{
  int m;
  
  for(m=0; (m<DIM); m++) {
    if (ir->ex[m].a)   sfree(ir->ex[m].a);
    if (ir->ex[m].phi) sfree(ir->ex[m].phi);
    if (ir->et[m].a)   sfree(ir->et[m].a);
    if (ir->et[m].phi) sfree(ir->et[m].phi);
  }

  sfree(ir->opts.nrdf);
  sfree(ir->opts.ref_t);
  sfree(ir->opts.annealing); 
  sfree(ir->opts.anneal_npoints); 
  sfree(ir->opts.anneal_time); 
  sfree(ir->opts.anneal_temp); 
  sfree(ir->opts.tau_t);
  sfree(ir->opts.acc);
  sfree(ir->opts.nFreeze);
  sfree(ir->opts.QMmethod);
  sfree(ir->opts.QMbasis);
  sfree(ir->opts.QMcharge);
  sfree(ir->opts.QMmult);
  sfree(ir->opts.bSH);
  sfree(ir->opts.CASorbitals);
  sfree(ir->opts.CASelectrons);
  sfree(ir->opts.SAon);
  sfree(ir->opts.SAoff);
  sfree(ir->opts.SAsteps);
  sfree(ir->opts.bOPT);
  sfree(ir->opts.bTS);

  if (ir->pull) {
    done_pull(ir->pull);
    sfree(ir->pull);
  }
}

static void init_ekinstate(ekinstate_t *eks)
{
  eks->ekinh_n = 0;
  eks->ekinh   = NULL;
  eks->dekindl = 0;
  eks->mvcos   = 0;
}

static void init_energyhistory(energyhistory_t *enh)
{
  enh->ener_ave = NULL;
  enh->ener_sum = NULL;
  enh->nener    = 0;
}

void init_gtc_state(t_state *state,int ngtc)
{
  int i;

  state->ngtc = ngtc;
  if (state->ngtc > 0) {
    snew(state->nosehoover_xi, state->ngtc);
    snew(state->therm_integral,state->ngtc);
    for(i=0; i<state->ngtc; i++) {
      state->nosehoover_xi[i]  = 0.0;
      state->therm_integral[i] = 0.0;
    }
  } else {
    state->nosehoover_xi  = NULL;
    state->therm_integral = NULL;
  }
}


void init_state(t_state *state,int natoms,int ngtc)
{
  int i;

  state->natoms = natoms;
  state->nrng   = 0;
  state->flags  = 0;
  state->lambda = 0;
  clear_mat(state->box);
  clear_mat(state->box_rel);
  clear_mat(state->boxv);
  clear_mat(state->pres_prev);
  init_gtc_state(state,ngtc);
  state->nalloc = state->natoms;
  if (state->nalloc > 0) {
    snew(state->x,state->nalloc);
    snew(state->v,state->nalloc);
  } else {
    state->x = NULL;
    state->v = NULL;
  }
  state->sd_X = NULL;
  state->cg_p = NULL;

  init_ekinstate(&state->ekinstate);

  init_energyhistory(&state->enerhist);

  state->ddp_count = 0;
  state->ddp_count_cg_gl = 0;
  state->cg_gl = NULL;
  state->cg_gl_nalloc = 0;
}

void done_state(t_state *state)
{
  if (state->nosehoover_xi) sfree(state->nosehoover_xi);
  if (state->x) sfree(state->x);
  if (state->v) sfree(state->v);
  if (state->sd_X) sfree(state->sd_X);
  if (state->cg_p) sfree(state->cg_p);
  state->nalloc = 0;
  if (state->cg_gl) sfree(state->cg_gl);
  state->cg_gl_nalloc = 0;
}

static void do_box_rel(t_inputrec *ir,matrix box_rel,matrix b,bool bInit)
{
  int d,d2;

  for(d=YY; d<=ZZ; d++) {
    for(d2=XX; d2<=(ir->epct==epctSEMIISOTROPIC ? YY : ZZ); d2++) {
      /* We need to check if this box component is deformed
       * or if deformation of another component might cause
       * changes in this component due to box corrections.
       */
      if (ir->deform[d][d2] == 0 &&
	  !(d == ZZ && d2 == XX && ir->deform[d][YY] != 0 &&
	    (b[YY][d2] != 0 || ir->deform[YY][d2] != 0))) {
	if (bInit) {
	  box_rel[d][d2] = b[d][d2]/b[XX][XX];
	} else {
	  b[d][d2] = b[XX][XX]*box_rel[d][d2];
	}
      }
    }
  }
}

void set_box_rel(t_inputrec *ir,t_state *state)
{
  /* Make sure the box obeys the restrictions before we fix the ratios */
  correct_box(NULL,0,state->box,NULL);

  clear_mat(state->box_rel);

  if (PRESERVE_SHAPE(*ir))
    do_box_rel(ir,state->box_rel,state->box,TRUE);
}

void preserve_box_shape(t_inputrec *ir,matrix box_rel,matrix b)
{
  if (PRESERVE_SHAPE(*ir))
    do_box_rel(ir,box_rel,b,FALSE);
}

void init_t_atoms(t_atoms *atoms, int natoms, bool bPdbinfo)
{
  atoms->nr=natoms;
  atoms->nres=0;
  snew(atoms->atomname,natoms);
  atoms->atomtype=NULL;
  atoms->atomtypeB=NULL;
  snew(atoms->resname,natoms);
  snew(atoms->atom,natoms);
  if (bPdbinfo)
    snew(atoms->pdbinfo,natoms);
  else
    atoms->pdbinfo=NULL;
}

void free_t_atoms(t_atoms *atoms,bool bFreeNames)
{
  int i;

  if (bFreeNames) {
    for(i=0; i<atoms->nr; i++) {
      sfree(*atoms->atomname[i]);
      *atoms->atomname[i]=NULL;
    }
    for(i=0; i<atoms->nres; i++) {
      sfree(*atoms->resname[i]);
      *atoms->resname[i]=NULL;
    }
  }
  sfree(atoms->atomname);
  /* Do we need to free atomtype and atomtypeB as well ? */
  sfree(atoms->resname);
  sfree(atoms->atom);
  if (atoms->pdbinfo)
    sfree(atoms->pdbinfo);
  atoms->nr=0; 
  atoms->nres=0;
}     
