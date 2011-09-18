/*
 * $Id: enxio.c,v 1.43.2.5 2009/01/18 23:06:26 lindahl Exp $
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "futil.h"
#include "string2.h"
#include "gmx_fatal.h"
#include "smalloc.h"
#include "gmxfio.h"
#include "enxio.h"
#include "vec.h"

void free_enxframe(t_enxframe *fr)
{
  int b;

  if (fr->e_alloc)
    sfree(fr->ener);
  if (fr->d_alloc) {
    sfree(fr->disre_rm3tav);
    sfree(fr->disre_rt);
  }
  for(b=0; b<fr->nblock; b++)
    sfree(fr->block[b]);
  sfree(fr->block);
  sfree(fr->b_alloc);
  sfree(fr->nr);
}

static void wr_ener_nms(FILE *out,int nre,char *nms[])
{
	int i,rc;

	rc = 0;
	
	rc = fprintf(out,"%d\n",nre);
	for(i=0; (i<nre) && rc>=0; i++)
		rc=fprintf(out,"%s\n",nms[i]);
	
	if(rc<0)
	{
		gmx_file("Cannot write energy names to file; maybe you are out of quota?");
	}
}

static void rd_ener_nms(FILE *in,int *nre,char ***nm)
{
  char line[256];
  int  i;
  
  fgets2(line,255,in);
  if (sscanf(line,"%d",nre) == 0) {
    *nre=0;
    return;
  }
  snew((*nm),*nre);
  for(i=0; (i< (*nre)); i++) {
    fgets2(line,255,in);
    trim(line);
    (*nm)[i]=strdup(line);
  }
}

static void edr_nms(int fp,int *nre,char ***nms)
{
  XDR  *xdr;
  bool bRead = gmx_fio_getread(fp);
  int  i;
  char **NM;

  xdr = gmx_fio_getxdr(fp);

  NM=*nms;
  
  if (!xdr_int(xdr,nre)) {
	  if(!bRead)
		  gmx_file("Cannot write energy names to file; maybe you are out of quota?");
    *nre=0;
    return;
  }
  if (*nre == -55555) {
    /* In the 4.1 edr format the first entry is the magic number */
    gmx_fatal(FARGS,"You are trying to read an edr file of GROMACS version 4.1 or later with GROMACS version 4.0");
  }
  if (NM == NULL) {
    snew(NM,*nre);
  }
  for(i=0; i<*nre; i++) {
    if (bRead && NM[i]) {
      sfree(NM[i]);
      NM[i] = NULL;
    }
	  if(!xdr_string(xdr,&(NM[i]),STRLEN))
		  gmx_file("Cannot write energy names to file; maybe you are out of quota?");
  }
  *nms=NM;
}

static bool do_eheader(int fp,t_enxframe *fr,bool *bOK)
{
  int  block,i,dum=0;
  bool bRead = gmx_fio_getread(fp);
  int  tempfix_nr=0;

  *bOK=TRUE;
  if (!do_real(fr->t))
  {
      return FALSE;
  }
  if (!do_int (fr->step)    ||
      !do_int (fr->nre)     ||
      !do_int (fr->ndisre)  || 
      !do_int (fr->nblock))
  {
	  *bOK = FALSE;
  }

  if (*bOK && bRead && fr->nblock>70) 
  {
    /* Temporary fix for intermediate file version, only used by B. Hess */
    tempfix_nr = fr->nblock;
    fr->nblock = 1;
  }
	
  if (*bOK && bRead && fr->nblock>fr->nr_alloc)
  {
    srenew(fr->nr,fr->nblock);
    srenew(fr->b_alloc,fr->nblock);
    srenew(fr->block,fr->nblock);
    for(i=fr->nr_alloc; i<fr->nblock; i++) {
      fr->block[i]   = NULL;
      fr->b_alloc[i] = 0;
    }
    fr->nr_alloc = fr->nblock;
  }
  if (tempfix_nr)
    fr->nr[0] = tempfix_nr;
  else 
  {
    for(block=0; block<fr->nblock && *bOK; block++)
		if (!do_int (fr->nr[block])) 
			*bOK = FALSE;
  }
	
  if(*bOK)
  {
	  if (!do_int (fr->e_size)   ||
		  !do_int (fr->d_size)   ||
		  !do_int (dum))
	  {
		  *bOK = FALSE;
	  }
  }
	
  return *bOK;
}

void do_enxnms(int fp,int *nre,char ***nms)
{
  bool bRead;
  
  bRead = gmx_fio_getread(fp);
  if (gmx_fio_getftp(fp) == efEDR) {
    gmx_fio_select(fp);
    edr_nms(fp,nre,nms);
  }
  else if (bRead)
    rd_ener_nms(gmx_fio_getfp(fp),nre,nms);
  else
    wr_ener_nms(gmx_fio_getfp(fp),*nre,*nms);
}

void close_enx(int fp)
{
  if(gmx_fio_close(fp) != 0)
  {
	  gmx_file("Cannot close energy file; it might be corrupt, or maybe you are out of quota?");  
  }
}

static bool empty_file(char *fn)
{
  FILE *fp;
  char dum;
  int  ret;
  bool bEmpty;
  
  fp = gmx_fio_fopen(fn,"r");
  ret = fread(&dum,sizeof(dum),1,fp);
  bEmpty = feof(fp);
  gmx_fio_fclose(fp);

  return bEmpty;
}

static int  framenr;
static real frametime;

int open_enx(char *fn,char *mode)
{
  int        fp,nre,i;
  char       **nm=NULL;
  t_enxframe *fr;
  bool       bDum=TRUE;

  if (mode[0]=='r') {
    fp=gmx_fio_open(fn,mode);
    gmx_fio_select(fp);
    gmx_fio_setprecision(fp,FALSE);
    do_enxnms(fp,&nre,&nm);
    snew(fr,1);
    do_eheader(fp,fr,&bDum);
	if(!bDum)
	{
		gmx_file("Cannot read energy file header. Corrupt file?");
	}
	  
    /* Now check whether this file is in single precision */
    if (((fr->e_size && (fr->nre == nre) && 
	  (nre*4*sizeof(float) == fr->e_size)) ||
	 (fr->d_size && 
	  (fr->ndisre*sizeof(float)*2+sizeof(int) == fr->d_size)))){
      fprintf(stderr,"Opened %s as single precision energy file\n",fn);
      for(i=0; (i<nre); i++)
	sfree(nm[i]);
      sfree(nm);
    }
    else {
      gmx_fio_rewind(fp);
      gmx_fio_select(fp);
      gmx_fio_setprecision(fp,TRUE);
      do_enxnms(fp,&nre,&nm);
      do_eheader(fp,fr,&bDum);
  	  if(!bDum)
	  {
		  gmx_file("Cannot write energy file header; maybe you are out of quota?");
	  }
		
      if (((fr->e_size && (fr->nre == nre) && 
	    (nre*4*sizeof(double) == fr->e_size)) ||
	   (fr->d_size && 
	    (fr->ndisre*sizeof(double)*2+sizeof(int) == fr->d_size))))
	fprintf(stderr,"Opened %s as double precision energy file\n",fn);
      else {
	if (empty_file(fn))
	  gmx_fatal(FARGS,"File %s is empty",fn);
	else
	  gmx_fatal(FARGS,"Energy file %s not recognized, maybe different CPU?",
		      fn);
      }
      for(i=0; (i<nre); i++)
	  sfree(nm[i]);
      sfree(nm);
    }
    free_enxframe(fr);
    sfree(fr);
    gmx_fio_rewind(fp);
  }
  else 
    fp = gmx_fio_open(fn,mode);
    
  framenr=0;
  frametime=0;

  return fp;
}

bool do_enx(int fp,t_enxframe *fr)
{
  int       i,block;
  bool      bRead,bOK,bOK1,bSane;
  real      tmp1,tmp2;

  bOK = TRUE;
  bRead = gmx_fio_getread(fp);
  if (!bRead) {  
    fr->e_size = fr->nre*sizeof(fr->ener[0].e)*4;
    fr->d_size = fr->ndisre*(sizeof(fr->disre_rm3tav[0]) + 
			     sizeof(fr->disre_rt[0]));
  }
  gmx_fio_select(fp);

  if (!do_eheader(fp,fr,&bOK)) {
    if (bRead) {
      fprintf(stderr,"\rLast energy frame read %d time %8.3f           ",
	      framenr-1,frametime);
      if (!bOK)
	fprintf(stderr,
		"\nWARNING: Incomplete energy frame: nr %d time %8.3f\n",
		framenr,fr->t);
    }
    else
	{
		gmx_file("Cannot write energy file header; maybe you are out of quota?");
	}
	  return FALSE;
  }
  if (bRead) {
    if ((framenr <   20 || framenr %   10 == 0) &&
	(framenr <  200 || framenr %  100 == 0) &&
	(framenr < 2000 || framenr % 1000 == 0))
      fprintf(stderr,"\rReading energy frame %6d time %8.3f           ",
	      framenr,fr->t);
    framenr++;
    frametime = fr->t;
  }
  /* Check sanity of this header */
  bSane = (fr->nre > 0 || fr->ndisre > 0);
  for(block=0; block<fr->nblock; block++)
    bSane = bSane || (fr->nr[block] > 0);
  if (!((fr->step >= 0) && bSane)) {
    fprintf(stderr,"\nWARNING: there may be something wrong with energy file %s\n",
	    gmx_fio_getname(fp));
    fprintf(stderr,"Found: step=%d, nre=%d, ndisre=%d, nblock=%d, time=%g.\n"
	    "Trying to skip frame expect a crash though\n",
	    fr->step,fr->nre,fr->ndisre,fr->nblock,fr->t);
  }
  if (bRead && fr->nre>fr->e_alloc) {
    srenew(fr->ener,fr->nre);
    for(i=fr->e_alloc; (i<fr->nre); i++) {
      fr->ener[i].e = fr->ener[i].eav =
	fr->ener[i].esum = fr->ener[i].e2sum = 0;
    }
    fr->e_alloc = fr->nre;
  }
  for(i=0; i<fr->nre; i++) {
    bOK = bOK && do_real(fr->ener[i].e);
    
    tmp1 = fr->ener[i].eav;

	bOK = bOK && do_real(tmp1);
    if (bRead)
      fr->ener[i].eav = tmp1;
    
    /* This is to save only in single precision (unless compiled in DP) */
    tmp2 = fr->ener[i].esum;
    bOK = bOK && do_real(tmp2);
    if (bRead)
      fr->ener[i].esum = tmp2;
    
    bOK = bOK && do_real(fr->ener[i].e2sum);
  }
  if (fr->ndisre) {
    if (bRead && fr->ndisre>fr->d_alloc) {
      srenew(fr->disre_rm3tav,fr->ndisre);
      srenew(fr->disre_rt,fr->ndisre);
      fr->d_alloc = fr->ndisre;
    }
    ndo_real(fr->disre_rm3tav,fr->ndisre,bOK1);
    bOK = bOK && bOK1;
    ndo_real(fr->disre_rt,fr->ndisre,bOK1);
    bOK = bOK && bOK1;
  }
  for(block=0; block<fr->nblock; block++) {
    if (bRead && fr->nr[block]>fr->b_alloc[block]) {
      srenew(fr->block[block],fr->nr[block]);
      fr->b_alloc[block] = fr->nr[block];
    }
    ndo_real(fr->block[block],fr->nr[block],bOK1);
    bOK = bOK && bOK1;
  }

  if(!bRead)
  {
	  if( gmx_fio_flush(fp) != 0)
	  {
		  gmx_file("Cannot write energy file; maybe you are out of quota?");
	  }
  }

  if (!bOK) {
    if (bRead) {
      fprintf(stderr,"\nLast energy frame read %d",
	      framenr-1);
      fprintf(stderr,"\nWARNING: Incomplete energy frame: nr %d time %8.3f\n",
	      framenr,fr->t);
    } else 
      gmx_fatal(FARGS,"could not write energies");
    return FALSE; 
  }
  
  return TRUE;
}

static real find_energy(char *name, int nre, char **enm, t_enxframe *fr)
{
  int i;

  for(i=0; i<nre; i++)
    if (strcmp(enm[i],name) == 0)
      return  fr->ener[i].e;

  gmx_fatal(FARGS,"Could not find energy term named '%s'",name);

  return 0;
}


void get_enx_state(char *fn, real t, gmx_groups_t *groups, t_inputrec *ir,
		   t_state *state)
{
  /* Should match the names in mdebin.c */
  static char *boxvel_nm[] = {
  "Box-Vel-XX", "Box-Vel-YY", "Box-Vel-ZZ",
  "Box-Vel-YX", "Box-Vel-ZX", "Box-Vel-ZY"
  };
  
  static char *pcouplmu_nm[] = {
    "Pcoupl-Mu-XX", "Pcoupl-Mu-YY", "Pcoupl-Mu-ZZ",
    "Pcoupl-Mu-YX", "Pcoupl-Mu-ZX", "Pcoupl-Mu-ZY"
  };
  int ind0[] = { XX,YY,ZZ,YY,ZZ,ZZ };
  int ind1[] = { XX,YY,ZZ,XX,XX,YY };

  int in,nre,nfr,i,ni,npcoupl;
  char       **enm=NULL,buf[STRLEN];
  t_enxframe *fr;

  in = open_enx(fn,"r");
  do_enxnms(in,&nre,&enm);
  snew(fr,1);
  nfr = 0;
  while ((nfr==0 || fr->t != t) && do_enx(in,fr)) {
    nfr++;
  }
  close_enx(in);
  fprintf(stderr,"\n");

  if (nfr == 0 || fr->t != t)
    gmx_fatal(FARGS,"Could not find frame with time %f in '%s'",t,fn);
  
  npcoupl = TRICLINIC(ir->compress) ? 6 : 3;
  if (ir->epc == epcPARRINELLORAHMAN) {
    clear_mat(state->boxv);
    for(i=0; i<npcoupl; i++) {
      state->boxv[ind0[i]][ind1[i]] =
	find_energy(boxvel_nm[i],nre,enm,fr);
    }
    fprintf(stderr,"\nREAD %d BOX VELOCITIES FROM %s\n\n",npcoupl,fn);
  }

  if (ir->etc == etcNOSEHOOVER) {
    for(i=0; i<state->ngtc; i++) {
      ni = groups->grps[egcTC].nm_ind[i];
      sprintf(buf,"Xi-%s",*(groups->grpname[ni]));
      state->nosehoover_xi[i] = find_energy(buf,nre,enm,fr);
    }
    fprintf(stderr,"\nREAD %d NOSE-HOOVER Xi's FROM %s\n\n",state->ngtc,fn);
  }
  
  free_enxframe(fr);
  sfree(fr);
}

