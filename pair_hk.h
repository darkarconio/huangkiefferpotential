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

#ifdef PAIR_CLASS

PairStyle(hk,PairHK)

#else

#ifndef LMP_PAIR_HK_H
#define LMP_PAIR_HK_H

#include "pair.h"

namespace LAMMPS_NS {

class PairHK : public Pair {
 public:
  PairHK(class LAMMPS *);
  virtual ~PairHK();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void init_style();

 protected:
  struct Param {
    double littleZ, littleN, bigA;
    double littleLambda, eta, kappa;
    double ac, bc;
    double az, bz;
    double C3, C4, A3, A4, gamma3, gamma4, thetaBar3, thetaBar4;
    double cut,cutsq,cutC,cutCsq;
    double bigC,bigC_kappa_eta;
    bool cdef;
    double bigZ, zeta;
    double lambda_eta;
    double abc,abz;
    double eabc,ebc,bebc;
    double eabz,ebz,bebz;
    int ielement,jelement,kelement;
  };

  double cutmax;                // max cutoff for all elements
  int nelements;                // # of unique elements
  char **elements;              // names of unique elements
  int ***elem2param;            // mapping from element triplets to parameters
  int *map;                     // mapping from atom types to elements
  int nparams;                  // # of stored parameter sets
  int maxparam;                 // max # of parameter sets
  Param *params;                // parameter set for an I-J-K interaction

  void allocate();
  void read_file(char *);
  void setup();
  void threebody(Param *, Param *, Param *, double, double, double *, double *, double *, double *, int, double &);
  void calcPhi(Param *, double, double *, double *, int, double &);
  void calcBigLambda(Param *, double, double, double *, double *, double *, double *, int, double &);
  void calcZeta(Param *, double);
  void calcBigZ(Param *, int);
  void calcBigC(Param *, Param *, Param *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair style Stillinger-Weber requires atom IDs

This is a requirement to use the SW potential.

E: Pair style Stillinger-Weber requires newton pair on

See the newton command.  This is a restriction to use the SW
potential.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

E: Cannot open Stillinger-Weber potential file %s

The specified SW potential file cannot be opened.  Check that the path
and name are correct.

E: Incorrect format in Stillinger-Weber potential file

Incorrect number of words per line in the potential file.

E: Illegal Stillinger-Weber parameter

One or more of the coefficients defined in the potential file is
invalid.

E: Potential file has duplicate entry

The potential file for a SW or Tersoff potential has more than
one entry for the same 3 ordered elements.

E: Potential file is missing an entry

The potential file for a SW or Tersoff potential does not have a
needed entry.

*/
