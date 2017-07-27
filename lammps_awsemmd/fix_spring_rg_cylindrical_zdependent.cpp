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
   Contributing author: Naveen Michaud-Agrawal (Johns Hopkins U)
                        Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "fix_spring_rg_cylindrical_zdependent.h"
#include "atom.h"
#include "update.h"
#include "group.h"
#include "respa.h"
#include "domain.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixSpringRGCylindricalZdependent::FixSpringRGCylindricalZdependent(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 7) error->all(FLERR,"Illegal fix spring/rg/cylindricalZdependent command");
  scalar_flag = 1;
  global_freq = 1;

  k = force->numeric(FLERR,arg[3]);
  rg0_flag = 0;
  if (strcmp(arg[4],"NULL") == 0) rg0_flag = 1;
  else rg0 = force->numeric(FLERR,arg[4]);

  k_bin = force->numeric(FLERR,arg[5]);
  memb_b = force->numeric(FLERR,arg[6]);

  dynamic_group_allow = 1;
  respa_level_support = 1;
  ilevel_respa = 0;
}

/* ---------------------------------------------------------------------- */

int FixSpringRGCylindricalZdependent::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSpringRGCylindricalZdependent::init()
{
  masstotall = group->mass(igroup);
  force_flag = 0;
  // if rg0 was specified as NULL, compute current Rg
  // only occurs on 1st run

  if (rg0_flag) {
    double xcm[3];
    group->z_based_xcm(igroup,masstotall,xcm, memb_b, k_bin);
    rg0 = group->z_based_cylindricalgyration(igroup,masstotall,xcm, memb_b, k_bin);
    rg0_flag = 0;
  }

  if (strstr(update->integrate_style,"respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixSpringRGCylindricalZdependent::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    ((Respa *) update->integrate)->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixSpringRGCylindricalZdependent::post_force(int vflag)
{
  // compute current Rg and center-of-mass

  force_flag = 0;
  double xcm[3];
  group->z_based_xcm(igroup,masstotall,xcm, memb_b, k_bin);
  double rg;
  rg = group->z_based_cylindricalgyration(igroup,masstotall, xcm,memb_b, k_bin);
  rg_total = rg;
  double xxcm[3];
  group->xcm(igroup,masstotall,xxcm);
  double rg_old = group->cylindricalgyration(igroup,masstotall,xxcm);

  // apply restoring force to atoms in group
  // f = -k*(r-r0)*mass/masstotal

  double dx,dy,z,term1,term2;

  double **f = atom->f;
  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  double masstotal;
  double massfrac;
  double mass_theta_prime_frac;
  double dxcm_dz, dycm_dz;
  double di_dz;
  double left, right;
  double theta, theta_prime;
  double unwrap[3];
  double forceAddX,forceAddY,forceAddZ;

  masstotal = group->z_based_mass(igroup, memb_b, k_bin);
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      domain->unmap(x[i],image[i],unwrap);
      dx = unwrap[0] - xcm[0];
      dy = unwrap[1] - xcm[1];
      z = unwrap[2];
      term1 = 2.0 * k * (1.0 - rg0/rg);
      theta = 0.5*(tanh(k_bin*(z+memb_b))-tanh(k_bin*(z-memb_b)));
      if (rmass) massfrac = rmass[i]*theta/masstotall;
      else  massfrac = mass[type[i]]*theta/masstotall;
      theta_prime = 0.5*k_bin*(-pow(tanh(k_bin*(z+memb_b)),2) + pow(tanh(k_bin*(z-memb_b)),2));
      mass_theta_prime_frac = mass[type[i]]*theta_prime/masstotall;
      forceAddX = term1*dx*massfrac;
      forceAddY = term1*dy*massfrac;
      f[i][0] -= forceAddX;
      f[i][1] -= forceAddY;
      left = dx*dx+dy*dy;
      forceAddZ = term1*0.5*(left)*mass_theta_prime_frac;
      f[i][2] -= forceAddZ;
    }
    total_energy = k * (rg_total - rg0) * (rg_total - rg0);
}

/* ---------------------------------------------------------------------- */

void FixSpringRGCylindricalZdependent::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}


double FixSpringRGCylindricalZdependent::compute_scalar()
{
  // only sum across procs one time

  // if (force_flag == 0) {
  //   MPI_Allreduce(total_energy_one,total_energy,1,MPI_DOUBLE,MPI_SUM,world);
  //   force_flag = 1;
  // }
  return total_energy;
}
