#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Header File for all diffeq assignment files
// TODO:
// * Free memory after
// * Manage output better
// * Automatic plotting of solutions
// * Runge Kutta
//    * Automatic time stepping method

// Coupled Pendulum Code

typedef struct parameters2
{
  double gamma_1;
  double gamma_2;
  double gamma_3;

  double c_1;
  double c_2;
  double c_3;

  double a_1;
  double a_2;
  double a_3;

  double omega_1;
  double omega_2;
  double omega_3;

  double phi_1;
  double phi_2;
  double phi_3;

  double k_1;
  double k_2;
  double k_3;

  double pulse_1;
  double pulse_2;
  double pulse_3;

  double sigma_1;
  double sigma_2;
  double sigma_3;

  double tau_1;
  double tau_2;
  double tau_3;

  double x_1_0;
  double x_2_0;
  double x_3_0;

  double x_1_DOT_0;
  double x_2_DOT_0;
  double x_3_DOT_0;

  double range;
  double delta_t;


} parameters2;

double** coupled_RK4(parameters2* info);
double current1(parameters2* info, double t);
double current2(parameters2* info, double t);
double current3(parameters2* info, double t);

typedef struct parameters
{
  double gamma;
  double c;
  double omega;
  double omega_0;
  double a;

  // Initial conditions?
  double X_0;
  double X_DOT_0;

  // Pass time parameters?
  double range;
  double delta_t;

} parameters;

// DIFFEQ1
double getxdotdot(parameters* info, double time, double x, double xdot);
double** diffeq1(parameters* info);

// RUNGE KUTTA
// DIFFEQ 3
double** rungekutta3(parameters* info);
double rgx3(parameters* info, double t, double x, double xdot);
double rgx2(parameters* info, double t, double x, double xdot);
double f1(parameters* info, double t, double x, double xdot);
double f2(parameters* info, double t, double x, double xdot);
double** rungekutta2(parameters* info);
double rg2(parameters* info, double t, double x);
double** diffeq22(parameters* info);
parameters* setdiffeq3(double gamma, double c, double omega, double omega_0, double a, double x_0, double x_dot_0, double range, double delta_t);
// Holds all information about

parameters* setdiffeq3(double gamma, double c, double omega, double omega_0, double a, double x_0, double x_dot_0, double range, double delta_t);

// Differential Equation 1
double** analyticdiffeq1(parameters* info);

//FFT
double** function(parameters* info);
double** powerSpectrum(double** array1, parameters* info);
double** dft(double** array1, parameters* info);
int DFT(int dir,int m,double *x1,double *y1);


// IV
double averageVelocity(double** array, parameters* info);
double** practicefunction(parameters* info);
double** ivcharacteristic(parameters* info, double c1, double c2, double precision);
void printArray(double** array, double len);
double** ivcoupled(parameters2* info, double c1, double c2, double precision);
/*
// FIX IMPLEMENTATION OF PARAMETER SETUP
// USE #IMPORT? #ifndef?
parameters* setdiffeq3(double gamma, double c, double omega, double a, double x_0, double x_dot_0, double range, double delta_t)
{
  parameters* info = malloc(sizeof(parameters));
  info->gamma = gamma;
  info->c = c;
  info->omega = omega;
  info->a = a;

  info->X_0 = x_0;
  info->X_DOT_0 = x_dot_0;

  info->range = range;
  info->delta_t = delta_t;

  return info;
}

*/


// Nonlinear oscillator with periodic driving
double getxdoubledot3(parameters* info, double t, double x, double xdot);
double** diffeq3(parameters* info);
