#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "diffeqheader.h"

// Runge-Kutta 4th Order
// Numerical solution to nonlinear oscillator with periodic driving
double f1(parameters* info, double t, double x, double xdot)
{
  return xdot;
}

double f2(parameters* info, double t, double x, double xdot)
{
  double gamma = info->gamma;
  double c = info->c;
  double omega = info->omega;
  double a = info->a;
  double omega_0 = info->omega_0;

  // Change equation here
  // c + (a*sin(omega*t)) - (gamma*xdot) - (omega_0 * omega_0 * x)
  return c + (a*sin(omega*t)) - (gamma*xdot) - sin(x);
}

double rgx2(parameters* info, double t, double x, double xdot)
{
  double k1, k2, k3, k4;
  double l1, l2, l3, l4;
  double sol;
  double delta_t = info->delta_t;

  k1 = delta_t * f1(info, t, x, xdot);
  l1 = delta_t * f2(info, t, x, xdot);

  k2 = delta_t * f1(info, t + (0.5*delta_t), x + (0.5*k1), xdot + (0.5*l2));
  l2 = delta_t * f2(info, t + (0.5*delta_t), x + (0.5*k1), xdot + (0.5*l2));

  k3 = delta_t * f1(info, t + (0.5*delta_t), x + (0.5*k2), xdot + (0.5*l2));
  l3 = delta_t * f2(info, t + (0.5*delta_t), x + (0.5*k2), xdot + (0.5*l2));

  k4 = delta_t * f1(info, t + delta_t, x + k3, xdot + l3);
  l4 = delta_t * f2(info, t + delta_t, x + k3, xdot + l3);

  sol = (double)((k1 + (2*k2) + (2*k3) + k4)/6);
  //printf("k = %f\n", sol);
  return sol;
}

double rgx3(parameters* info, double t, double x, double xdot)
{
  double k1, k2, k3, k4;
  double l1, l2, l3, l4;
  double sol;

  double delta_t = info->delta_t;

  k1 = delta_t * f1(info, t, x, xdot);
  l1 = delta_t * f2(info, t, x, xdot);

  k2 = delta_t * f1(info, t + (0.5*delta_t), x + (0.5*k1), xdot + (0.5*l2));
  l2 = delta_t * f2(info, t + (0.5*delta_t), x + (0.5*k1), xdot + (0.5*l2));

  k3 = delta_t * f1(info, t + (0.5*delta_t), x + (0.5*k2), xdot + (0.5*l2));
  l3 = delta_t * f2(info, t + (0.5*delta_t), x + (0.5*k2), xdot + (0.5*l2));

  k4 = delta_t * f1(info, t + delta_t, x + k3, xdot + l3);
  l4 = delta_t * f2(info, t + delta_t, x + k3, xdot + l3);

  sol = (double)((l1 + (2*l2) + (2*l3) + l4)/6);
  //printf("l = %f\n", sol);
  return sol;
}

double** rungekutta3(parameters* info)
{
  int i, j;

  double X_0 = info->X_0;
  double X_DOT_0 = info->X_DOT_0;

  // Parameters
  double gamma = info->gamma;
  double c = info->c;
  double omega = info->omega;
  double a = info->a;

  // Time Steps
  double delta_t = info->delta_t;
  double range = info->range;
  int len = (int) range / delta_t;

  // Allocate memory for vector
  double **array = malloc(sizeof(double*) * 3);
  for (i = 0; i < 3; i++)
  {
    array[i] = malloc(sizeof(double) * len);
  }

  // Initialization of initial conditions
  array[0][0] = 0.0; // Time
  array[1][0] = X_0; // Position
  array[2][0] = X_DOT_0; // Velocity

  // Runge Kutta Routine
  for (i = 1; i < len; i++)
  {
    array[0][i] = array[0][i-1] + delta_t;
    array[1][i] = array[1][i-1] + rgx2(info, array[0][i-1], array[1][i-1], array[2][i-1]);
    array[2][i] = array[2][i-1] + rgx3(info, array[0][i-1], array[1][i-1], array[2][i-1]);
  }

  return array;
}

// Problem here? Showing chaos?
double averageVelocity(double** array, parameters* info)
{
  int i, j = 0, start;
  double sum = 0, avg = 0;

  int len = (int) info->range / info->delta_t;

  // Start count after transient
  start = (int) len * (9.0 / 10.0);
  //printf("%d %d\n", len, start);

  for (i = 0; i < len; i++)
  {
    sum += array[2][i];
    j++;
  }

  avg = (double) sum / j;
  return avg;
}

double** ivcharacteristic(parameters* info, double c1, double c2, double precision)
{
  double **array, **rk4;
  double range, constant;
  int i;

  int len = (int) info->range / info->delta_t;

  // Josephson phase voltage relation
  constant = 1.0;

  range = (int) (c2 - c1) / precision;
  array = malloc(sizeof(double*) * 2);

  for (i = 0; i < 2; i++)
  {
    array[i] = malloc(sizeof(double) * range);
  }

  // Initialization of driving force array
  array[0][0] = c1;
  for (i = 1; i < range; i++)
  {
    array[0][i] = array[0][i-1] + precision;
  }

  // Computation
  for (i = 0; i < range; i++)
  {
    info->c = array[0][i];
    array[1][i] = averageVelocity(rk4 = rungekutta3(info), info) * constant;


    // CHECKING TRANSIENT DELETE HERE
    if ((info->c > 0.499000) && (info->c < 0.501000))
    {
      // printArray(rungekutta3(info), len);
      // printf("Avg vel = %f\n", averageVelocity(rungekutta3(info), info));
    }

    // Setting initial conditions dynamically to reduce transient
    // info->X_0 = rk4[1][len-1];
    // info->X_DOT_0 = rk4[2][len-1];

    free(rk4[0]);
    free(rk4[1]);
    free(rk4[2]);
  }

  /*
  // Hysterisis, reverse computation
  for (i = range - 1; i > -1; i--)
  {
    info->c = array[0][i]

  }

  */

  return array;
}

double** practicefunction(parameters* info)
{
  int i;
  double** array;

  double delta_t = info->delta_t;
  double range = info->range;
  int len = (int) range / delta_t;

  array = malloc(sizeof(double*) * 3);
  for (i = 0; i < 3; i++)
  {
    array[i] = malloc(sizeof(double) * len);
  }

  array[0][0] = 0.0;
  for (i = 1; i < len; i++)
  {
    array[0][i] = array[0][i-1] + delta_t;
    array[1][i] = sin(array[0][i]);
    array[2][i] = cos(array[0][i]);
  }

  return array;
}
