#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "diffeqheader.h"

// Numerical solution to set of three coupled oscillators
// Use of Runge Kutta 4th Order method
//
// x1" + gx1' + sin(x1) = C1 + A1sin((w1*t) + p1) + k1(x2 - x1)
// x2" + gx2' + sin(x2) = C2 + A2sin((w2*t) + p2) + k2(x3 - 2x2 +x1)
// x3" + gx3' + sin(x3) = C3 + A3sin((w3*t) + p3) + k3(x2 - x3)
//
// Time | X1 | X1' | X2 | X2'| X3 | X3'

// Current Functions
double current1(parameters2* info, double t)
{
  // Resolves nan error
  if (info->sigma_1 == 0.0 || info->pulse_1 == 0.0)
    return info->c_1;
    
  return info->c_1 + ((info->pulse_1 / sqrt(2*M_PI*info->sigma_1*info->sigma_1)) * (exp( (-1 * pow(t - info->tau_1 , 2)) / (2 * info->sigma_1 * info->sigma_1) )));
}

double current2(parameters2* info, double t)
{
  if (info->sigma_2 == 0.0 || info->pulse_2 == 0.0)
    return info->c_2;

  return info->c_2 + ((info->pulse_2 / sqrt(2*M_PI*info->sigma_2*info->sigma_2)) * (exp( (-1 * pow(t - info->tau_2 , 2)) / (2 * info->sigma_2 * info->sigma_2) )));
}

double current3(parameters2* info, double t)
{
  if (info->sigma_3 == 0.0 || info->pulse_3 == 0.0)
    return info->c_3;

  return info->c_3 + ((info->pulse_3 / sqrt(2*M_PI*info->sigma_3*info->sigma_3)) * (exp( (-1 * pow(t - info->tau_3 , 2)) / (2 * info->sigma_3 * info->sigma_3) )));
}

// System of ODEs
double x_1_dotf(parameters2* info, double t, double x_1 , double x_1_dot, double x_2, double x_2_dot, double x_3, double x_3_dot)
{
  return x_1_dot;
}

double x_1_dotdotf(parameters2* info, double t, double x_1 , double x_1_dot, double x_2, double x_2_dot, double x_3, double x_3_dot)
{
  return (-1*info->gamma_1*x_1_dot) - (sin(x_1)) + current1(info, t) + (info->k_1*(x_2 - x_1));
}

double x_2_dotf(parameters2* info, double t, double x_1 , double x_1_dot, double x_2, double x_2_dot, double x_3, double x_3_dot)
{
  return x_2_dot;
}

double x_2_dotdotf(parameters2* info, double t, double x_1 , double x_1_dot, double x_2, double x_2_dot, double x_3, double x_3_dot)
{
  return (-1*info->gamma_2*x_2_dot) - (sin(x_2)) + current2(info, t) + (info->k_2*(x_3 - (2*x_2) + x_1));
}

double x_3_dotf(parameters2* info, double t, double x_1 , double x_1_dot, double x_2, double x_2_dot, double x_3, double x_3_dot)
{
  return x_3_dot;
}

double x_3_dotdotf(parameters2* info, double t, double x_1 , double x_1_dot, double x_2, double x_2_dot, double x_3, double x_3_dot)
{
  return (-1*info->gamma_3*x_3_dot) - (sin(x_3)) + current3(info, t) + (info->k_3*(x_2 - x_3));
}

double* rg4(parameters2* info, double t, double x_1 , double x_1_dot, double x_2, double x_2_dot, double x_3, double x_3_dot)
{
  double* sol = malloc(sizeof(double) * 6);
  double delta_t;

  double k1_1, k1_2, k1_3, k1_4;
  double k2_1, k2_2, k2_3, k2_4;
  double k3_1, k3_2, k3_3, k3_4;

  double l1_1, l1_2, l1_3, l1_4;
  double l2_1, l2_2, l2_3, l2_4;
  double l3_1, l3_2, l3_3, l3_4;

  delta_t = info->delta_t;

  k1_1 = delta_t * x_1_dotf(info, t, x_1, x_1_dot, x_2, x_2_dot, x_3, x_3_dot);
  l1_1 = delta_t * x_1_dotdotf(info, t, x_1, x_1_dot, x_2, x_2_dot, x_3, x_3_dot);
  k2_1 = delta_t * x_2_dotf(info, t, x_1, x_1_dot, x_2, x_2_dot, x_3, x_3_dot);
  l2_1 = delta_t * x_2_dotdotf(info, t, x_1, x_1_dot, x_2, x_2_dot, x_3, x_3_dot);
  k3_1 = delta_t * x_3_dotf(info, t, x_1, x_1_dot, x_2, x_2_dot, x_3, x_3_dot);
  l3_1 = delta_t * x_3_dotdotf(info, t, x_1, x_1_dot, x_2, x_2_dot, x_3, x_3_dot);

  k1_2 = delta_t * x_1_dotf(info, t + (0.5*delta_t), x_1 + (0.5*k1_1), x_1_dot + (0.5*l1_1), x_2 + (0.5*k2_1), x_2_dot + (0.5*l2_1), x_3 + (0.5*k3_1), x_3_dot + (0.5*l3_1));;
  l1_2 = delta_t * x_1_dotdotf(info, t + (0.5*delta_t), x_1 + (0.5*k1_1), x_1_dot + (0.5*l1_1), x_2 + (0.5*k2_1), x_2_dot + (0.5*l2_1), x_3 + (0.5*k3_1), x_3_dot + (0.5*l3_1));
  k2_2 = delta_t * x_2_dotf(info, t + (0.5*delta_t), x_1 + (0.5*k1_1), x_1_dot + (0.5*l1_1), x_2 + (0.5*k2_1), x_2_dot + (0.5*l2_1), x_3 + (0.5*k3_1), x_3_dot + (0.5*l3_1));
  l2_2 = delta_t * x_2_dotdotf(info, t + (0.5*delta_t), x_1 + (0.5*k1_1), x_1_dot + (0.5*l1_1), x_2 + (0.5*k2_1), x_2_dot + (0.5*l2_1), x_3 + (0.5*k3_1), x_3_dot + (0.5*l3_1));
  k3_2 = delta_t * x_3_dotf(info, t + (0.5*delta_t), x_1 + (0.5*k1_1), x_1_dot + (0.5*l1_1), x_2 + (0.5*k2_1), x_2_dot + (0.5*l2_1), x_3 + (0.5*k3_1), x_3_dot + (0.5*l3_1));
  l3_2 = delta_t * x_3_dotdotf(info, t + (0.5*delta_t), x_1 + (0.5*k1_1), x_1_dot + (0.5*l1_1), x_2 + (0.5*k2_1), x_2_dot + (0.5*l2_1), x_3 + (0.5*k3_1), x_3_dot + (0.5*l3_1));

  k1_3 = delta_t * x_1_dotf(info, t + (0.5*delta_t), x_1 + (0.5*k1_2), x_1_dot + (0.5*l1_2), x_2 + (0.5*k2_2), x_2_dot + (0.5*l2_2), x_3 + (0.5*k3_2), x_3_dot + (0.5*l3_2));
  l1_3 = delta_t * x_1_dotdotf(info, t + (0.5*delta_t), x_1 + (0.5*k1_2), x_1_dot + (0.5*l1_2), x_2 + (0.5*k2_2), x_2_dot + (0.5*l2_2), x_3 + (0.5*k3_2), x_3_dot + (0.5*l3_2));
  k2_3 = delta_t * x_2_dotf(info, t + (0.5*delta_t), x_1 + (0.5*k1_2), x_1_dot + (0.5*l1_2), x_2 + (0.5*k2_2), x_2_dot + (0.5*l2_2), x_3 + (0.5*k3_2), x_3_dot + (0.5*l3_2));
  l2_3 = delta_t * x_2_dotdotf(info, t + (0.5*delta_t), x_1 + (0.5*k1_2), x_1_dot + (0.5*l1_2), x_2 + (0.5*k2_2), x_2_dot + (0.5*l2_2), x_3 + (0.5*k3_2), x_3_dot + (0.5*l3_2));
  k3_3 = delta_t * x_3_dotf(info, t + (0.5*delta_t), x_1 + (0.5*k1_2), x_1_dot + (0.5*l1_2), x_2 + (0.5*k2_2), x_2_dot + (0.5*l2_2), x_3 + (0.5*k3_2), x_3_dot + (0.5*l3_2));
  l3_3 = delta_t * x_3_dotdotf(info, t + (0.5*delta_t), x_1 + (0.5*k1_2), x_1_dot + (0.5*l1_2), x_2 + (0.5*k2_2), x_2_dot + (0.5*l2_2), x_3 + (0.5*k3_2), x_3_dot + (0.5*l3_2));

  k1_4 = delta_t * x_1_dotf(info, t + delta_t, x_1 + k1_3, x_1_dot + l1_3, x_2 + k2_3, x_2_dot + l2_3, x_3 + k3_3, x_3_dot + l3_3);
  l1_4 = delta_t * x_1_dotdotf(info, t + delta_t, x_1 + k1_3, x_1_dot + l1_3, x_2 + k2_3, x_2_dot + l2_3, x_3 + k3_3, x_3_dot + l3_3);
  k2_4 = delta_t * x_2_dotf(info, t + delta_t, x_1 + k1_3, x_1_dot + l1_3, x_2 + k2_3, x_2_dot + l2_3, x_3 + k3_3, x_3_dot + l3_3);
  l2_4 = delta_t * x_2_dotdotf(info, t + delta_t, x_1 + k1_3, x_1_dot + l1_3, x_2 + k2_3, x_2_dot + l2_3, x_3 + k3_3, x_3_dot + l3_3);
  k3_4 = delta_t * x_3_dotf(info, t + delta_t, x_1 + k1_3, x_1_dot + l1_3, x_2 + k2_3, x_2_dot + l2_3, x_3 + k3_3, x_3_dot + l3_3);
  l3_4 = delta_t * x_3_dotdotf(info, t + delta_t, x_1 + k1_3, x_1_dot + l1_3, x_2 + k2_3, x_2_dot + l2_3, x_3 + k3_3, x_3_dot + l3_3);

  sol[0] = (double) ((k1_1 + (2*k1_2) + (2*k1_3) + k1_4) / 6.0);
  sol[1] = (double) ((l1_1 + (2*l1_2) + (2*l1_3) + l1_4) / 6.0);
  sol[2] = (double) ((k2_1 + (2*k2_2) + (2*k2_3) + k2_4) / 6.0);
  sol[3] = (double) ((l2_1 + (2*l2_2) + (2*l2_3) + l2_4) / 6.0);
  sol[4] = (double) ((k3_1 + (2*k3_2) + (2*k3_3) + k3_4) / 6.0);
  sol[5] = (double) ((l3_1 + (2*l3_2) + (2*l3_3) + l3_4) / 6.0);

  return sol;
}
// Main function
double** coupled_RK4(parameters2* info)
{
  double** array;
  double* sol;
  double range, delta_t;
  int i, j, len;

  delta_t = info->delta_t;
  len = (int) (info->range) / (info->delta_t);

  array = malloc(sizeof(double*) * 7);
  for (i = 0; i < 7; i++)
  {
    array[i] = malloc(sizeof(double) * len);
  }

  // Initialization of intial conditions
  array[0][0] = 0.0;
  array[1][0] = info->x_1_0;
  array[2][0] = info->x_1_DOT_0;
  array[3][0] = info->x_2_0;
  array[4][0] = info->x_2_DOT_0;
  array[5][0] = info->x_3_0;
  array[6][0] = info->x_3_DOT_0;

  // Runge-Kutta
  for (i = 1; i < len; i++)
  {
    array[0][i] = array[0][i-1] + delta_t;

    sol = rg4(info, array[0][i-1], array[1][i-1], array[2][i-1], array[3][i-1],
                    array[4][i-1], array[5][i-1], array[6][i-1]);

    array[1][i] = array[1][i-1] + sol[0];
    array[2][i] = array[2][i-1] + sol[1];
    array[3][i] = array[3][i-1] + sol[2];
    array[4][i] = array[4][i-1] + sol[3];
    array[5][i] = array[5][i-1] + sol[4];
    array[6][i] = array[6][i-1] + sol[5];

    free(sol);
  }

  return array;
}

// Calculates average after specified transient period
double averagexdotcoupled(parameters2* info, double* array)
{
  int i, j, start, len;
  double sum;

  len = (int) info->range / info->delta_t;
  start = (int) len * (9.0 / 10.0);

  sum = 0, j = 0;
  for (i = start; i < len; i++)
  {
    sum += array[i];
    j++;
  }

  return (double) sum / j;
}

// Computes I-V curves for coupled oscillators
double** ivcoupled(parameters2* info, double c1, double c2, double precision)
{
  double** array, **rk4;
  double range, constant;
  int i;

  // Use for Josephson Phase-Voltage Relation
  constant = 1.0;

  // Allocate memory for time series
  range = (int) (c2 - c1) / precision;
  array = malloc(sizeof(double*) * 4);
  for (i = 0; i < 4; i++)
  {
    array[i] = malloc(sizeof(double) * range);
  }

  // Initialization of driving array
  array[0][0] = c1;
  for (i = 1; i < range; i++)
  {
    array[0][i] = array[0][i-1] + precision;
  }

  // Computation
  for (i = 0; i < range; i++)
  {
    info->c_1 = array[0][i];
    info->c_2 = array[0][i];
    info->c_3 = array[0][i];

    rk4 = coupled_RK4(info);

    array[1][i] = averagexdotcoupled(info, rk4[2]) * constant;
    array[2][i] = averagexdotcoupled(info, rk4[4]) * constant;
    array[3][i] = averagexdotcoupled(info, rk4[6]) * constant;
  }

  return array;
}
