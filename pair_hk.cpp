/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Michal Plucinski (Dalhousie University)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_hk.h"
#include "atom.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "fenv.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024
#define DELTA 4

/* ---------------------------------------------------------------------- */

PairHK::PairHK(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;

  nelements = 0;
  elements = NULL;
  nparams = maxparam = 0;
  params = NULL;
  elem2param = NULL;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairHK::~PairHK()
{
  if (elements)
    for (int i = 0; i < nelements; i++) delete [] elements[i];
  delete [] elements;
  memory->destroy(params);
  memory->destroy(elem2param);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    delete [] map;
  }
}

/* ---------------------------------------------------------------------- */

void PairHK::compute(int eflag, int vflag)
{
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  int i,j,k,ii,jj,kk,inum,jnum,jnumm1,itag,jtag;
  int itype,jtype,ktype,iparam,jparam,ijparam,jiparam,ikparam,ijkparam;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,rsq1,rsq2;
  double delr1[3],delr2[3],fj[3],fk[3];
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
//  int *tag = atom->tag;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over full neighbor list of my atoms
   char str[128];

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
//    itag = tag[i];
    itype = map[type[i]];
    if (itype == 1) continue;
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    iparam = elem2param[itype][itype][itype];
//   sprintf(str,"i %i",i);
//   error->warning(FLERR,str);
    calcBigZ(&params[iparam],i);

    // two-body interactions, skip half of them

    jlist = firstneigh[i];
    jnum = numneigh[i];

    jnumm1 = jnum - 1;

    for (jj = 0; jj < jnumm1; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = map[type[j]];
      ijparam = elem2param[itype][jtype][jtype];
      delr1[0] = x[j][0] - xtmp;
      delr1[1] = x[j][1] - ytmp;
      delr1[2] = x[j][2] - ztmp;
      rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];
      if (rsq1 > params[ijparam].cutsq || jtype == 0) continue;
//   sprintf(str,"j %i x %g y %g z %g",j,x[j][0],x[j][1],x[j][2]);
//   error->warning(FLERR,str);

      for (kk = jj+1; kk < jnum; kk++) {
        k = jlist[kk];
        k &= NEIGHMASK;
        ktype = map[type[k]];
        ikparam = elem2param[itype][ktype][ktype];
        ijkparam = elem2param[itype][jtype][ktype];

        delr2[0] = x[k][0] - xtmp;
        delr2[1] = x[k][1] - ytmp;
        delr2[2] = x[k][2] - ztmp;
        rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
        if (rsq2 > params[ikparam].cutsq || ktype == 0) continue;
//   sprintf(str,"k %i x %g y %g z %g",k, x[k][0], x[k][1],x[k][2]);
//   error->warning(FLERR,str);

        threebody(&params[iparam],&params[ijparam],&params[ikparam],rsq1,rsq2,delr1,delr2,fj,fk,eflag,evdwl);
//   sprintf(str,"eng %g fj1 %g fj2 %g fj3 %g fk1 %g fk2 %g fk3 %g", evdwl,fj[0],fj[1],fj[2],fk[0],fk[1],fk[2]);
//   error->warning(FLERR,str);

        f[i][0] -= fj[0] + fk[0];
        f[i][1] -= fj[1] + fk[1];
        f[i][2] -= fj[2] + fk[2];
        f[j][0] += fj[0];
        f[j][1] += fj[1];
        f[j][2] += fj[2];
        f[k][0] += fk[0];
        f[k][1] += fk[1];
        f[k][2] += fk[2];

        if (evflag) ev_tally3(i,j,k,evdwl,0.0,fj,fk,delr1,delr2);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairHK::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  map = new int[n+1];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairHK::settings(int narg, char **arg)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairHK::coeff(int narg, char **arg)
{
  int i,j,n;

  if (!allocated) allocate();

  if (narg != 3 + atom->ntypes)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // insure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // read args that map atom types to elements in potential file
  // map[i] = which element the Ith atom type is, -1 if NULL
  // nelements = # of unique elements
  // elements = list of element names

  if (elements) {
    for (i = 0; i < nelements; i++) delete [] elements[i];
    delete [] elements;
  }
  elements = new char*[atom->ntypes];
  for (i = 0; i < atom->ntypes; i++) elements[i] = NULL;

  nelements = 0;
  for (i = 3; i < narg; i++) {
    if (strcmp(arg[i],"NULL") == 0) {
      map[i-2] = -1;
      continue;
    }
    for (j = 0; j < nelements; j++)
      if (strcmp(arg[i],elements[j]) == 0) break;
    map[i-2] = j;
    if (j == nelements) {
      n = strlen(arg[i]) + 1;
      elements[j] = new char[n];
      strcpy(elements[j],arg[i]);
      nelements++;
    }
  }

  // read potential file and initialize potential parameters

  read_file(arg[2]);
  setup();

  // clear setflag since coeff() called once with I,J = * *

  n = atom->ntypes;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements

  int count = 0;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      if (map[i] >= 0 && map[j] >= 0) {
        setflag[i][j] = 1;
        count++;
      }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairHK::init_style()
{
//  if (atom->tag_enable == 0)
//    error->all(FLERR,"Pair style Stillinger-Weber requires atom IDs");
//  if (force->newton_pair == 0)
//    error->all(FLERR,"Pair style Stillinger-Weber requires newton pair on");

  // need a full neighbor list

  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairHK::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  return cutmax;
}

/* ---------------------------------------------------------------------- */

void PairHK::read_file(char *file)
{
  int params_per_line = 23;
  char **words = new char*[params_per_line+1];

  memory->sfree(params);
  params = NULL;
  nparams = maxparam = 0;

  // open file on proc 0

  FILE *fp;
  if (comm->me == 0) {
    fp = fopen(file,"r");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open Huang-Kieffer potential file %s",file);
      error->one(FLERR,str);
    }
  }

  // read each set of params from potential file
  // one set of params can span multiple lines
  // store params if all 3 element tags are in element list

  int n,nwords,ielement,jelement,kelement;
  char line[MAXLINE],*ptr;
  int eof = 0;

  while (1) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fp);
      if (ptr == NULL) {
        eof = 1;
        fclose(fp);
      } else n = strlen(line) + 1;
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // strip comment, skip line if blank

    if (ptr = strchr(line,'#')) *ptr = '\0';
    nwords = atom->count_words(line);
    if (nwords == 0) continue;

    // concatenate additional lines until have params_per_line words

    while (nwords < params_per_line) {
      n = strlen(line);
      if (comm->me == 0) {
        ptr = fgets(&line[n],MAXLINE-n,fp);
        if (ptr == NULL) {
          eof = 1;
          fclose(fp);
        } else n = strlen(line) + 1;
      }
      MPI_Bcast(&eof,1,MPI_INT,0,world);
      if (eof) break;
      MPI_Bcast(&n,1,MPI_INT,0,world);
      MPI_Bcast(line,n,MPI_CHAR,0,world);
      if (ptr = strchr(line,'#')) *ptr = '\0';
      nwords = atom->count_words(line);
    }

    if (nwords != params_per_line)
      error->all(FLERR,"Incorrect format in Huang-Kieffer potential file");

    // words = ptrs to all words in line

    nwords = 0;
    words[nwords++] = strtok(line," \t\n\r\f");
    while (words[nwords++] = strtok(NULL," \t\n\r\f")) continue;

    // ielement,jelement,kelement = 1st args
    // if all 3 args are in element list, then parse this line
    // else skip to next entry in file

    for (ielement = 0; ielement < nelements; ielement++)
      if (strcmp(words[0],elements[ielement]) == 0) break;
    if (ielement == nelements) continue;
    for (jelement = 0; jelement < nelements; jelement++)
      if (strcmp(words[1],elements[jelement]) == 0) break;
    if (jelement == nelements) continue;
    for (kelement = 0; kelement < nelements; kelement++)
      if (strcmp(words[2],elements[kelement]) == 0) break;
    if (kelement == nelements) continue;

    // load up parameter settings and error check their values

    if (nparams == maxparam) {
      maxparam += DELTA;
      params = (Param *) memory->srealloc(params,maxparam*sizeof(Param),
                                          "pair:params");
    }

    params[nparams].ielement = ielement;
    params[nparams].jelement = jelement;
    params[nparams].kelement = kelement;
    params[nparams].littleZ = atof(words[3]);
    params[nparams].littleN = atof(words[4]);
    params[nparams].bigA = atof(words[5]);
    params[nparams].littleLambda = atof(words[6]);
    params[nparams].eta = atof(words[7]);
    params[nparams].kappa = atof(words[8]);
    params[nparams].ac = atof(words[9]);
    params[nparams].bc = atof(words[10]);
    params[nparams].az = atof(words[11]);
    params[nparams].bz = atof(words[12]);
    params[nparams].C3 = atof(words[13]);
    params[nparams].C4 = atof(words[14]);
    params[nparams].A3 = atof(words[15]);
    params[nparams].A4 = atof(words[16]);
    params[nparams].gamma3 = atof(words[17]);
    params[nparams].gamma4 = atof(words[18]);
    params[nparams].thetaBar3 = atof(words[19]);
    params[nparams].thetaBar4 = atof(words[20]);
    params[nparams].cut = atof(words[21]);
    params[nparams].cutC = atof(words[22]);
    params[nparams].cdef = false;

    if (params[nparams].C3 < 0.0 || params[nparams].littleN < 0.0 || params[nparams].C4 < 0.0 || 
        params[nparams].bigA < 0.0 || params[nparams].littleLambda < 0.0 || params[nparams].A3 < 0.0 || 
        params[nparams].eta < 0.0 || params[nparams].kappa < 0.0 || params[nparams].A4 < 0.0 || 
        params[nparams].ac < 0.0 || params[nparams].az < 0.0 || params[nparams].gamma3 < 0.0 || 
        params[nparams].az < 0.0 || params[nparams].bz < 0.0 || params[nparams].gamma4 < 0.0 ||
        params[nparams].thetaBar3 < 0.0 ||  params[nparams].thetaBar4 < 0.0 || params[nparams].cut < 0.0 ||
        params[nparams].cutC < 0.0)
      error->all(FLERR,"Illegal Huang-Kieffer parameter");

    nparams++;
  }

  delete [] words;
}

/* ---------------------------------------------------------------------- */

void PairHK::setup()
{
  int i,j,k,m,n;
  double rtmp;

  // set elem2param for all triplet combinations
  // must be a single exact match to lines read from file
  // do not allow for ACB in place of ABC

  memory->destroy(elem2param);
  memory->create(elem2param,nelements,nelements,nelements,"pair:elem2param");

  for (i = 0; i < nelements; i++)
    for (j = 0; j < nelements; j++)
      for (k = 0; k < nelements; k++) {
        n = -1;
        for (m = 0; m < nparams; m++) {
          if (i == params[m].ielement && j == params[m].jelement &&
              k == params[m].kelement) {
            if (n >= 0) error->all(FLERR,"Potential file has duplicate entry");
            n = m;
          }
        }
        if (n < 0) error->all(FLERR,"Potential file is missing an entry");
        elem2param[i][j][k] = n;
      }


  // compute parameter values derived from inputs

  for (m = 0; m < nparams; m++) {
    rtmp = params[m].cut;
    params[m].cutsq = rtmp * rtmp;
    params[m].cutCsq = pow(params[m].cutC,2);
    params[m].lambda_eta = params[m].littleLambda * params[n].eta; 
    params[m].abc = params[m].ac * params[m].bc;
    params[m].eabc = exp( params[m].ac * params[m].bc );
    params[m].ebc = exp( params[m].bc );
    params[m].bebc = params[m].bc * exp( params[m].bc );
    params[m].abz = params[m].az * params[m].bz;
    params[m].eabz = exp( params[m].az * params[m].bz );
    params[m].ebz = exp( params[m].bz ); 
    params[m].bebz = params[m].bz * exp( params[m].bz );
  }

  // set cutmax to max of all params

  cutmax = 0.0;
  for (m = 0; m < nparams; m++) {
    rtmp = sqrt(params[m].cut);
    if (rtmp > cutmax) cutmax = rtmp;
  }
}

/* --------------------------------------------------------------------- */

void PairHK::calcBigC (Param *parami, Param *paramj, Param *paramij)
{
//   char str[128];
//   sprintf(str,"A %g zi %g ni %g zj %g nj %g kappa %g eta %g",paramij->bigA, parami->littleZ, parami->littleN, paramj->littleZ, paramj->littleN, paramij->kappa, paramij->eta);
//   error->warning(FLERR,str);
   
   paramij->bigC = paramij->bigA*( 1 + (parami->littleZ/parami->littleN) + (paramj->littleZ/paramj->littleN) );
   paramij->bigC_kappa_eta = paramij->bigC * paramij->kappa / paramij->eta;
}

/* ---------------------------------------------------------------------- */

void PairHK::calcBigLambda(Param * parami, double rsq1, double rsq2, double *delr1, double *delr2,
                           double *fj, double *fk, int eflag, double &eng)
{
   int i;
   double r1,r2;
   double bigLambda = 0;
   double rinv12,costheta,theta,exps;
   double delZ3,delZsq3,deltheta3,delthetasq3,expZ3,exptheta3,lambda3,consts3;
   double delZ4,delZsq4,deltheta4,delthetasq4,expZ4,exptheta4,lambda4,consts4;

   
   r1 = sqrt(rsq1);
   r2 = sqrt(rsq2);
   
   rinv12 = 1.0/(r1*r2);
   costheta = (delr1[0]*delr2[0] + delr1[1]*delr2[1] + delr1[2]*delr2[2]) * rinv12;
   theta = acos(costheta);

   delZ3 = 3 - parami->bigZ;
   delZsq3 = delZ3*delZ3;
   expZ3 = exp( -1 * parami->A3 * delZsq3 );

   deltheta3 = parami->thetaBar3 - theta;
   delthetasq3 = deltheta3*deltheta3;
   exptheta3 = exp ( -1 * parami->gamma3 * delthetasq3 );

   lambda3 = parami->C3 * expZ3 * exptheta3;
   bigLambda += lambda3;

   delZ4 = 4 - parami->bigZ;
   delZsq4 = delZ4*delZ4;
   expZ4 = exp( -1 * parami->A4 * delZsq4 );

   deltheta4 = parami->thetaBar4 - theta;
   delthetasq4 = deltheta4*deltheta4;
   exptheta4 = exp ( -1 * parami->gamma4 * delthetasq4 );

   lambda4 = parami->C4 * expZ4 * exptheta4;
   bigLambda += lambda4;
  
//   char str[128];
//   sprintf(str,"Z %g", parami->bigZ);
 //  error->warning(FLERR,str);
   
   if (eflag) eng = bigLambda;
   
   consts3 = lambda3 * 2 * parami->A3 * delZ3 * parami->bebz;
   consts4 = lambda4 * 2 * parami->A4 * delZ4 * parami->bebz;
   exps = exp(parami->az + r1) / pow( parami->eabz + parami->ebz * exp(r1) ,2) / r1;
   
   for (i=0;i<3;i++) {
      fj[i] = (consts3 * exps * delr1[i]) + (consts4 * exps * delr1[i]);
/*      if (i==2){
      }
  */ }
   //sprintf(str,"expsj %g c3 %g c4 %g delr11 %g fj1 %g delr12 %g fj2 %g recal %g", exps,consts3,consts4,delr1[1],fj[1],delr1[2],fj[2], (consts3 * exps * delr1[2]) + (consts4 * exps * delr1[2]));
   //error->warning(FLERR,str);
   
   exps = exp(parami->az + r2) / pow( parami->eabz + parami->ebz * exp(r2) ,2) / r2;

   for (i=0;i<3;i++) {
      fk[i] = (consts3 * exps * delr2[i]) + (consts4 * exps * delr2[i]);
    /*  if (i==2){
   sprintf(str,"expsk %g c3 %g c4 %g delr21 %g fk1 %g delr22 %g fk2 %g recal %g", exps,consts3,consts4,delr2[1],fk[1],delr2[i],fk[i],(consts3 * exps * delr2[i]) + (consts4 * exps * delr2[i]));
   error->warning(FLERR,str);
      }
  */ }
//   sprintf(str,"expsk %g c3 %g c4 %g delr21 %g fk1 %g delr22 %g fk2 %g recal %g", exps,consts3,consts4,delr2[1],fk[1],delr2[2],fk[2],(consts3 * exps * delr2[2]) + (consts4 * exps * delr2[2]));
//   error->warning(FLERR,str);
}

/* ---------------------------------------------------------------------- */

void PairHK::calcZeta(Param * paramij, double rsq)
{
    paramij->zeta = pow( 1 + exp( (sqrt(rsq) - paramij->ac) * paramij->bc ) , -1);
/*   char str[128];
//   if(paramij->zeta < .001){
   sprintf(str,"zeta %g rsq %g ac %g bc %g",paramij->zeta,rsq,paramij->ac,paramij->bc);
   error->warning(FLERR,str);
//   }*/
}

/* ---------------------------------------------------------------------- */

void PairHK::calcBigZ(Param * parami, int i)
{
  int j,jj,jnum;
  int itype,jtype,ijparam,jparam,iparam;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double **x = atom->x;
  int nlocal = atom->nlocal;
  int *type = atom->type;
  parami->bigZ = 0;

  ilist = list->ilist;
  numneigh = list->numneigh;
  jnum = numneigh[i];
  firstneigh = list->firstneigh;
  itype = map[type[i]];

  xtmp = x[i][0];
  ytmp = x[i][1];
  ztmp = x[i][2];

  // loop over full neighbor list of my atoms

  for (jj = 0; jj < jnum; jj++) {
    jlist = firstneigh[i];
    j = jlist[jj];
    j &= NEIGHMASK;
    delx = xtmp - x[j][0];
    dely = ytmp - x[j][1];
    delz = ztmp - x[j][2];
    rsq = delx*delx + dely*dely + delz*delz;
    
    jtype = map[type[j]];
    iparam = elem2param[itype][itype][itype];
    jparam = elem2param[jtype][jtype][jtype];
    ijparam = elem2param[itype][jtype][jtype];
    
    if (params[ijparam].cdef == false && rsq <= params[ijparam].cutsq) {
      calcBigC(&params[iparam], &params[jparam], &params[ijparam]);
      calcZeta(&params[ijparam], rsq);
      params[ijparam].cdef = true;
    }
    
    if (rsq > params[ijparam].cutCsq) continue;
    parami->bigZ += pow( 1 + exp( (sqrt(rsq) - parami->az) * parami->bz ) , -1);
  }
}

/* ---------------------------------------------------------------------- */

void PairHK::calcPhi(Param *p, double rsq, double *delr, double *f, int eflag, double &eng)
{
   double r,phi,fcon,exps;
   
   r = sqrt(rsq);
   
   phi = -1 * p->zeta * p->bigC_kappa_eta * exp( p->eta * (p->littleLambda - r) );
   if (eflag) eng = phi;
/*   char str[128];
   if(eng == 0){
   sprintf(str,"zeta %g cke %g eta %g llambda %g r %g phi %g",p->zeta,p->bigC_kappa_eta,p->eta,p->littleLambda, r, phi);
   error->warning(FLERR,str);
   }
  */ 
   exps = p->eabc + p->ebc * exp(r);
   fcon = p->bigC_kappa_eta / r / pow(exps,2) * exp( p->abc + p->lambda_eta + p->eta*r);
   fcon *= p->eta * exps + p->bebc * exp(r);
   for (int i=0;i<3;i++) f[i] = fcon * delr[i];
}

/* ---------------------------------------------------------------------- */

void PairHK::threebody(Param * parami, Param * paramij, Param * paramik,
                       double rsq1, double rsq2, double * delr1, double * delr2,
                       double * fj, double * fk, int eflag, double &eng)
{
   double bigLambda, phi_ij, phi_ik, sumPhi;
   double fphi_ij[3], fphi_ik[3], fbigLambda_ij[3], fbigLambda_ik[3];
   
   calcPhi(paramij, rsq1, delr1, fphi_ij, eflag, phi_ij);
   calcPhi(paramik, rsq2, delr2, fphi_ik, eflag, phi_ik);
   calcBigLambda(parami, rsq1, rsq2, delr1, delr2, fbigLambda_ij, fbigLambda_ik, eflag, bigLambda);
   
   sumPhi =  phi_ij + phi_ik;
   if (eflag) eng = sumPhi * bigLambda;

/*   char str[128];
   if(eng == 0){
   sprintf(str,"phij %g phik %g lambda %g",phi_ij,phi_ik,bigLambda);
   error->warning(FLERR,str);
   }*/
   
   for (int i=0;i<3;i++) {
      fj[i] = ( fphi_ij[i] + phi_ij ) * bigLambda + sumPhi * fbigLambda_ij[i];
      fk[i] = ( fphi_ik[i] + phi_ik ) * bigLambda + sumPhi * fbigLambda_ik[i];
   }
}
