/*
 * $Id: tpbconv.c,v 1.56.2.3 2009/02/13 08:25:37 hess Exp $
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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include "index.h"
#include "gmx_fatal.h"
#include "string2.h"
#include "sysstuff.h"
#include "smalloc.h"
#include "macros.h"
#include "names.h"
#include "typedefs.h"
#include "tpxio.h"
#include "trnio.h"
#include "enxio.h"
#include "readir.h"
#include "statutil.h"
#include "copyrite.h"
#include "futil.h"
#include "vec.h"
#include "mtop_util.h"
#include "random.h"

#define RANGECHK(i,n) if ((i)>=(n)) gmx_fatal(FARGS,"Your index file contains atomnumbers (e.g. %d)\nthat are larger than the number of atoms in the tpr file (%d)",(i),(n))

static bool *bKeepIt(int gnx,int natoms,atom_id index[])
{
  bool *b;
  int  i;
  
  snew(b,natoms);
  for(i=0; (i<gnx); i++) {
    RANGECHK(index[i],natoms);
    b[index[i]] = TRUE;
  }
  
  return b;
}

static atom_id *invind(int gnx,int natoms,atom_id index[])
{
  atom_id *inv;
  int     i;
  
  snew(inv,natoms);
  for(i=0; (i<gnx); i++) {
    RANGECHK(index[i],natoms);
    inv[index[i]] = i;
  }
  
  return inv;
}

static void reduce_block(bool bKeep[],t_block *block,
			 const char *name)
{
  atom_id *index;
  int i,j,newi,newj;
  
  snew(index,block->nr);
  
  newi = newj = 0;
  for(i=0; (i<block->nr); i++) {
    for(j=block->index[i]; (j<block->index[i+1]); j++) {
      if (bKeep[j]) {
	newj++;
      }
    }
    if (newj > index[newi]) {
      newi++;
      index[newi] = newj;
    }
  }
  
  fprintf(stderr,"Reduced block %8s from %6d to %6d index-, %6d to %6d a-entries\n",
	  name,block->nr,newi,block->index[block->nr],newj);
  block->index = index;
  block->nr    = newi;
}

static void reduce_blocka(atom_id invindex[],bool bKeep[],t_blocka *block,
			  const char *name)
{
  atom_id *index,*a;
  int i,j,k,newi,newj;
  
  snew(index,block->nr);
  snew(a,block->nra);
  
  newi = newj = 0;
  for(i=0; (i<block->nr); i++) {
    for(j=block->index[i]; (j<block->index[i+1]); j++) {
      k = block->a[j];
      if (bKeep[k]) {
	a[newj] = invindex[k];
	newj++;
      }
    }
    if (newj > index[newi]) {
      newi++;
      index[newi] = newj;
    }
  }
  
  fprintf(stderr,"Reduced block %8s from %6d to %6d index-, %6d to %6d a-entries\n",
	  name,block->nr,newi,block->nra,newj);
  block->index = index;
  block->a     = a;
  block->nr    = newi;
  block->nra   = newj;
}

static void reduce_rvec(int gnx,atom_id index[],rvec vv[])
{
  rvec *ptr;
  int  i;
  
  snew(ptr,gnx);
  for(i=0; (i<gnx); i++)
    copy_rvec(vv[index[i]],ptr[i]);
  for(i=0; (i<gnx); i++)
    copy_rvec(ptr[i],vv[i]);
  sfree(ptr);
}

static void reduce_atom(int gnx,atom_id index[],t_atom atom[],char ***atomname,
			int *nres, char ***resname)
{
  t_atom *ptr;
  char   ***aname,***rname;
  int    i,nr;
  
  snew(ptr,gnx);
  snew(aname,gnx);
  snew(rname,atom[index[gnx-1]].resnr+1);
  for(i=0; (i<gnx); i++) {
    ptr[i]   = atom[index[i]];
    aname[i] = atomname[index[i]];
  }
  nr=-1;   
  for(i=0; (i<gnx); i++) {
    atom[i]     = ptr[i];
    atomname[i] = aname[i];
    if ((i==0) || (atom[i].resnr != atom[i-1].resnr)) {
      nr++;
      rname[nr]=resname[atom[i].resnr];
    }
    atom[i].resnr=nr;
  }
  nr++;
  for(i=0; (i<nr); i++)
    resname[i]=rname[i];
  *nres=nr;

  sfree(aname);
  sfree(ptr);
  sfree(rname);
}

static void reduce_ilist(atom_id invindex[],bool bKeep[],
			 t_ilist *il,int nratoms,char *name)
{
  t_iatom *ia;
  int i,j,newnr;
  bool bB;

  if (il->nr) {  
    snew(ia,il->nr);
    newnr=0;
    for(i=0; (i<il->nr); i+=nratoms+1) {
      bB = TRUE;
      for(j=1; (j<=nratoms); j++) {
	bB = bB && bKeep[il->iatoms[i+j]];
      }
      if (bB) {
	ia[newnr++] = il->iatoms[i];
	for(j=1; (j<=nratoms); j++)
	  ia[newnr++] = invindex[il->iatoms[i+j]];
      }
    }
    fprintf(stderr,"Reduced ilist %8s from %6d to %6d entries\n",
	    name,il->nr/(nratoms+1),
	  newnr/(nratoms+1));
    
    il->nr = newnr;
    for(i=0; (i<newnr); i++)
      il->iatoms[i] = ia[i];
    
    sfree(ia);
  }
}

static void reduce_topology_x(int gnx,atom_id index[],
			      gmx_mtop_t *mtop,rvec x[],rvec v[])
{
  t_topology top;
  bool    *bKeep;
  atom_id *invindex;
  int     i;
  
  top = gmx_mtop_t_to_t_topology(mtop); 
  bKeep    = bKeepIt(gnx,top.atoms.nr,index);
  invindex = invind(gnx,top.atoms.nr,index);
  
  reduce_block(bKeep,&(top.cgs),"cgs");
  reduce_block(bKeep,&(top.mols),"mols");
  reduce_blocka(invindex,bKeep,&(top.excls),"excls");
  reduce_rvec(gnx,index,x);
  reduce_rvec(gnx,index,v);
  reduce_atom(gnx,index,top.atoms.atom,top.atoms.atomname,
	      &(top.atoms.nres),top.atoms.resname);

  for(i=0; (i<F_NRE); i++) {
    reduce_ilist(invindex,bKeep,&(top.idef.il[i]),
		 interaction_function[i].nratoms,
		 interaction_function[i].name);
  }
    
  top.atoms.nr = gnx;

  mtop->nmoltype = 1;
  snew(mtop->moltype,mtop->nmoltype);
  mtop->moltype[0].name  = mtop->name;
  mtop->moltype[0].atoms = top.atoms;
  for(i=0; i<F_NRE; i++) {
    mtop->moltype[0].ilist[i] = top.idef.il[i];
  }
  mtop->moltype[0].atoms = top.atoms;
  mtop->moltype[0].cgs   = top.cgs;
  mtop->moltype[0].excls = top.excls;

  mtop->nmolblock = 1;
  snew(mtop->molblock,mtop->nmolblock);
  mtop->molblock[0].type       = 0;
  mtop->molblock[0].nmol       = 1;
  mtop->molblock[0].natoms_mol = top.atoms.nr;
  mtop->molblock[0].nposres_xA = 0;
  mtop->molblock[0].nposres_xB = 0;
  
  mtop->natoms                 = top.atoms.nr;
}

static void zeroq(int n,atom_id index[],gmx_mtop_t *mtop)
{
  int mt,i;
  
  for(mt=0; mt<mtop->nmoltype; mt++) {
    for(i=0; (i<mtop->moltype[mt].atoms.nr); i++) {
      mtop->moltype[mt].atoms.atom[index[i]].q = 0;
      mtop->moltype[mt].atoms.atom[index[i]].qB = 0;
    }
  }
}

int main (int argc, char *argv[])
{
  static char *desc[] = {
    "tpbconv can edit run input files in four ways.[PAR]"
    "[BB]1st.[bb] by modifying the number of steps in a run input file",
    "with option [TT]-nsteps[tt] or option [TT]-runtime[tt].[PAR]",
    "[BB]2st.[bb] (OBSOLETE) by creating a run input file",
    "for a continuation run when your simulation has crashed due to e.g.",
    "a full disk, or by making a continuation run input file.",
    "This option is obsolete, since mdrun now writes and reads",
    "checkpoint files.",
    "Note that a frame with coordinates and velocities is needed.",
    "When pressure and/or Nose-Hoover temperature coupling is used",
    "an energy file can be supplied to get an exact continuation",
    "of the original run.[PAR]",
    "[BB]3nd.[bb] by creating a tpx file for a subset of your original",
    "tpx file, which is useful when you want to remove the solvent from",
    "your tpx file, or when you want to make e.g. a pure Ca tpx file.",
    "[BB]WARNING: this tpx file is not fully functional[bb].",
    "[BB]4rd.[bb] by setting the charges of a specified group",
    "to zero. This is useful when doing free energy estimates",
    "using the LIE (Linear Interaction Energy) method."
  };

  char         *top_fn,*frame_fn;
  int          fp,fp_ener=-1;
  t_trnheader head;
  int          i,frame,run_step,nsteps_org;
  real         run_t,state_t;
  bool         bOK,bNsteps,bExtend,bUntil,bTime,bTraj;
  bool         bFrame,bUse,bSel,bNeedEner,bReadEner,bScanEner;
  gmx_mtop_t   mtop;
  t_atoms      atoms;
  t_inputrec   *ir,*irnew=NULL;
  t_gromppopts *gopts;
  t_state      state;
  rvec         *newx=NULL,*newv=NULL,*tmpx,*tmpv;
  matrix       newbox;
  int          gnx;
  char         *grpname;
  atom_id      *index=NULL;
  int          nre;
  char         **enm=NULL;
  t_enxframe   *fr_ener=NULL;
  t_filenm fnm[] = {
    { efTPX, NULL,  NULL,    ffREAD  },
    { efTRN, "-f",  NULL,    ffOPTRD },
    { efENX, "-e",  NULL,    ffOPTRD },
    { efNDX, NULL,  NULL,    ffOPTRD },
    { efTPX, "-o",  "tpxout",ffWRITE }
  };
#define NFILE asize(fnm)

  /* Command line options */
  static int  nsteps_req = -1;
  static real runtime_req = -1;
  static real start_t = -1.0, extend_t = 0.0, until_t = 0.0;
  static bool bContinuation = TRUE,bZeroQ = FALSE;
  static t_pargs pa[] = {
    { "-nsteps",        FALSE, etINT,  {&nsteps_req},
      "Change the number of steps" },
    { "-runtime",       FALSE, etREAL, {&runtime_req},
      "Set the run time (ps)" },
    { "-time",          FALSE, etREAL, {&start_t}, 
      "Continue from frame at this time (ps) instead of the last frame" },
    { "-extend",        FALSE, etREAL, {&extend_t}, 
      "Extend runtime by this amount (ps)" },
    { "-until",         FALSE, etREAL, {&until_t}, 
      "Extend runtime until this ending time (ps)" },
    { "-zeroq",         FALSE, etBOOL, {&bZeroQ},
      "Set the charges of a group (from the index) to zero" },
    { "-cont",          FALSE, etBOOL, {&bContinuation},
      "For exact continuation, the constraints should not be solved before the first step" }
  };
  int nerror = 0;
  
  CopyRight(stdout,argv[0]);
  
  /* Parse the command line */
  parse_common_args(&argc,argv,0,NFILE,fnm,asize(pa),pa,
		    asize(desc),desc,0,NULL);

  bNsteps = (nsteps_req >= 0 || runtime_req >= 0);
  bExtend = opt2parg_bSet("-extend",asize(pa),pa);
  bUntil  = opt2parg_bSet("-until",asize(pa),pa);
  bTime   = opt2parg_bSet("-time",asize(pa),pa);
  bTraj   = (opt2bSet("-f",NFILE,fnm) || bTime);

  top_fn = ftp2fn(efTPX,NFILE,fnm);
  fprintf(stderr,"Reading toplogy and shit from %s\n",top_fn);
  
  snew(ir,1);
  read_tpx_state(top_fn,&run_step,&run_t,ir,&state,NULL,&mtop);

  if (bTraj) {
    fprintf(stderr,"\n"
	    "NOTE: Reading the state from trajectory is an obsolete feaure of tpbconv.\n"
	    "      Continuation should be done by loading a checkpoint file with mdrun -cpi\n"
	    "      This guarantees that all state variables are transferred.\n"
	    "      tpbconv is now only useful for increasing nsteps,\n"
	    "      but even that can often be avoided by using mdrun -maxh\n"
	    "\n");

    if (ir->bContinuation != bContinuation)
      fprintf(stderr,"Modifying ir->bContinuation to %s\n",
	      bool_names[bContinuation]);
    ir->bContinuation = bContinuation;
    

    bNeedEner = (ir->epc == epcPARRINELLORAHMAN || ir->etc == etcNOSEHOOVER);
    bReadEner = (bNeedEner && ftp2bSet(efENX,NFILE,fnm));
    bScanEner = (bReadEner && !bTime);
    
    if (ir->epc != epcNO || EI_SD(ir->eI) || ir->eI == eiBD) {
      fprintf(stderr,"NOTE: The simulation uses pressure coupling and/or stochastic dynamics.\n"
	      "tpbconv can not provide binary identical continuation.\n"
	      "If you want that, supply a checkpoint file to mdrun\n\n");
    }

    if (EI_SD(ir->eI) || ir->eI == eiBD) {
      fprintf(stderr,"\nChanging ld-seed from %d ",ir->ld_seed);
      ir->ld_seed = make_seed();
      fprintf(stderr,"to %d\n\n",ir->ld_seed);
    }

    frame_fn = ftp2fn(efTRN,NFILE,fnm);
    fprintf(stderr,
	    "\nREADING COORDS, VELS AND BOX FROM TRAJECTORY %s...\n\n",
	    frame_fn);
    
    fp = open_trn(frame_fn,"r");
    if (bScanEner) {
      fp_ener = open_enx(ftp2fn(efENX,NFILE,fnm),"r");
      do_enxnms(fp_ener,&nre,&enm);
      snew(fr_ener,1);
      fr_ener->t = -1e-12;
    }

    /* Now scan until the last set of x and v (step == 0)
     * or the ones at step step.
     */
    bFrame = TRUE;
    frame  = 0;
    while (bFrame) {
      bFrame = fread_trnheader(fp,&head,&bOK);
      if (bOK && frame == 0) {
	if (mtop.natoms != head.natoms) 
	  gmx_fatal(FARGS,"Number of atoms in Topology (%d) "
		      "is not the same as in Trajectory (%d)\n",
		      mtop.natoms,head.natoms);
	snew(newx,head.natoms);
	snew(newv,head.natoms);
      }
      bFrame = bFrame && bOK;
      if (bFrame) {
	
	bOK = fread_htrn(fp,&head,newbox,newx,newv,NULL);
      }
      bFrame = bFrame && bOK;
      bUse = FALSE;
      if (bFrame &&
	  (head.x_size) && (head.v_size || !EI_STATE_VELOCITY(ir->eI))) {
	bUse = TRUE;
	if (bScanEner) {
	  /* Read until the energy time is >= the trajectory time */
	  while (fr_ener->t < head.t && do_enx(fp_ener,fr_ener));
	  bUse = (fr_ener->t == head.t);
	}
	if (bUse) {
	  tmpx    = newx;
	  newx    = state.x;
	  state.x = tmpx;
	  tmpv    = newv;
	  newv    = state.v;
	  state.v = tmpv;
	  run_t        = head.t;
	  run_step     = head.step;
	  state.lambda = head.lambda;
	  copy_mat(newbox,state.box);
	}
      }
      if (bFrame || !bOK) {
	fprintf(stderr,"\r%s %s frame %6d: step %6d time %8.3f",
		bUse ? "Read   " : "Skipped",ftp2ext(fn2ftp(frame_fn)),
		frame,head.step,head.t);
	frame++;
	if (bTime && (head.t >= start_t))
	  bFrame = FALSE;
      }
    }
    if (bScanEner) {
      close_enx(fp_ener);
      free_enxframe(fr_ener);
      for(i=0; i<nre; i++)
	sfree(enm[i]);
      sfree(enm);
    }
    close_trn(fp);
    fprintf(stderr,"\n");

    if (!bOK)
      fprintf(stderr,"%s frame %d (step %d, time %g) is incomplete\n",
	      ftp2ext(fn2ftp(frame_fn)),frame-1,head.step,head.t);
    fprintf(stderr,"\nUsing frame of step %d time %g\n",run_step,run_t);

    if (bNeedEner) {
      if (bReadEner) {
	get_enx_state(ftp2fn(efENX,NFILE,fnm),run_t,&mtop.groups,ir,&state);
      } else {
	fprintf(stderr,"\nWARNING: The simulation uses %s temperature and/or %s pressure coupling,\n"
		"         the continuation will only be exact when an energy file is supplied\n\n",
		ETCOUPLTYPE(etcNOSEHOOVER),
		EPCOUPLTYPE(epcPARRINELLORAHMAN));
      }
    }
  }

  if (bNsteps) {
    if (nsteps_req < 0) {
      if (!EI_DYNAMICS(ir->eI)) {
	gmx_fatal(FARGS,"Can not set the run time with integrator '%s'",
		  EI(ir->eI));
      }
      nsteps_req = (int)(runtime_req/ir->delta_t + 0.5);
    }
    fprintf(stderr,"Setting nsteps to %d\n",nsteps_req);
    ir->nsteps = nsteps_req;
  } else {
    /* Determine total number of steps remaining */
    if (bExtend) {
      ir->nsteps = ir->nsteps - (run_step - ir->init_step) + (int)(extend_t/ir->delta_t + 0.5);
      printf("Extending remaining runtime of by %g ps (now %d steps)\n",
	     extend_t,ir->nsteps);
    }
    else if (bUntil) {
      printf("nsteps = %d, run_step = %d, current_t = %g, until = %g\n",
	     ir->nsteps,run_step,run_t,until_t);
      ir->nsteps = (int)((until_t - run_t)/ir->delta_t + 0.5);
      printf("Extending remaining runtime until %g ps (now %d steps)\n",
	     until_t,ir->nsteps);
    }
    else {
      ir->nsteps -= run_step - ir->init_step; 
      /* Print message */
      printf("%d steps (%g ps) remaining from first run.\n",
	   ir->nsteps,ir->nsteps*ir->delta_t);
	  
    }
  }

  if (bNsteps || bZeroQ || (ir->nsteps > 0)) {
    ir->init_step = run_step;
    
    if (ftp2bSet(efNDX,NFILE,fnm) || 
	!(bNsteps || bExtend || bUntil || bTraj)) {
      atoms = gmx_mtop_global_atoms(&mtop);
      get_index(&atoms,ftp2fn_null(efNDX,NFILE,fnm),1,
		&gnx,&index,&grpname);
      if (!bZeroQ) {
	bSel = (gnx != state.natoms);
	for (i=0; ((i<gnx) && (!bSel)); i++)
	  bSel = (i!=index[i]);
      }
      else
	bSel = FALSE;
      if (bSel) {
	fprintf(stderr,"Will write subset %s of original tpx containing %d "
		"atoms\n",grpname,gnx);
	reduce_topology_x(gnx,index,&mtop,state.x,state.v);
	state.natoms = gnx;
      } 
      else if (bZeroQ) {
	zeroq(gnx,index,&mtop);
	fprintf(stderr,"Zero-ing charges for group %s\n",grpname);
      }
      else
	fprintf(stderr,"Will write full tpx file (no selection)\n");
    }    

    state_t = ir->init_t + ir->init_step*ir->delta_t;
    fprintf(stderr,"Writing statusfile with starting step %10d and length %10d steps...\n",
	    ir->init_step,ir->nsteps);
    fprintf(stderr,"                                 time %10.3f and length %10.3f ps\n",
	    state_t,ir->nsteps*ir->delta_t);
    write_tpx_state(opt2fn("-o",NFILE,fnm),0,state_t,ir,&state,&mtop);
  }
  else
    printf("You've simulated long enough. Not writing tpr file\n");
	      
  thanx(stderr);
  
  return 0;
}
