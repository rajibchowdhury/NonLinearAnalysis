#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "diffeqheader.h"

double particularsol1(parameters* info, double t)
{
  double sol;
  double c1, c2, det;
  double gamma = info->gamma;
  double omega = info->omega;
  double omega_0 = info->omega_0;
  double c = info->c;

  // Solve by undetermined coefficients
  det = (-1*gamma*gamma*omega*omega) - (pow((omega_0*omega_0) - (omega*omega), 2));
  c1 = (c*gamma*omega) / det;
  c2 = (-1*c*((omega_0*omega_0) - (omega*omega))) / det;

  sol = (c1*cos(omega*t)) + (c2*sin(omega*t));
  return sol;
}

double generalsol1(parameters* info, double t)
{
  double c1, c2, r1, r2;
  double a, b;
  double sol;

  double gamma = info->gamma;
  double omega_0 = info->omega_0;
  sol = 0.0;


  // Complex Roots
  if ((gamma*gamma) < (4*omega_0*omega_0))
  {
    a = (double) (-1*gamma)/2;
    b = (double) (sqrt((4*omega_0*omega_0) - (gamma*gamma))) / 2;
    c1 = info->X_0;
    c2 = (double)((info->X_DOT_0) - ((info->X_0) * a)) / b;

    sol = (c1*exp(a*t)*cos(b*t)) + (c2*exp(a*t)*sin(b*t));
    return sol;
  }

  // Real Roots
  else if ((gamma*gamma) > (4*omega_0*omega_0))
  {
    r1 = (double)((-1*gamma) + sqrt((gamma*gamma) - (4*omega_0*omega_0))) / 2;
    r2 = (double)((-1*gamma) - sqrt((gamma*gamma) - (4*omega_0*omega_0))) / 2;

    c1 = (double) (((info->X_0)*r2) - info->X_DOT_0) / (r2 - r1);
    c2 = (double) (info->X_DOT_0 - ((info->X_0)*r1)) / (r2 - r1);

    sol = (c1*exp(r1*t)) + (c2*exp(r2*t));
    return sol;
  }

  // Real and repeated roots
  else if ((gamma*gamma) == (4*omega_0*omega_0))
  {
    r1 = (double) (-1*gamma) / 2;
    c1 = info->X_0;
    c2 = (info->X_DOT_0) - ((info->X_0) * r1);

    sol = (c1*exp(r1*t)) + (c2*t*exp(r1*t));
    return sol;
  }

  return sol;
}

double** analyticdiffeq1(parameters* info)
{
  int i, j;
  double r1, r2;
  double const1, const2;

  double gamma = info->gamma;
  double c = info->c;
  double omega = info->omega;
  double omega_0 = info->omega_0;

  double X_0 = info->X_0;
  double X_DOT_0 = info->X_DOT_0;

  double range = info->range;
  double delta_t = info->delta_t;
  int len = (int) range / delta_t;

  double **array = malloc(sizeof(double*) * 3);
  for (i = 0; i < 3; i++)
  {
    array[i] = malloc(sizeof(double) * len);
  }

  // Initial Conditions
  array[0][0] = 0.0;
  array[1][0] = X_0;
  array[2][0] = X_DOT_0;

  for (i = 1; i < len; i++)
  {
    array[0][i] = array[0][i-1] + delta_t;
    array[1][i] = generalsol1(info, array[0][i]) + particularsol1(info, array[0][i]);
  }

  return array;
}
