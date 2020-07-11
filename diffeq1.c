#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Linear oscillator
// TODO:
// Write all parameters/initial conditions
// into a struct to be passed to the function
// Header file
// Output to a textfile

// typedef struct parameters
// {
//   double gamma;
//   double omega;
//   double cc;
//} parameters;

double getxdotdot(double time, double x, double xdot)
{
  double gamma = 1.1;
  double omega = 0.5;
  double c = 0.5;

  return (c*sin(omega*time)) - (gamma * xdot) - (omega * omega * x);
}

int main(void)
{
  int i, j;
  //parameters* constants = malloc(sizeof(parameters*));

  // Constants
  //constants->gamma = 0.3;
  //constants->omega = 0.5;
  //constants->cc = 0.3;
  //printf("C is %f\n", constants->cc);

  double gamma = 1.1;
  double omega = 0.5;
  double c = 0.5;

  // Initial Conditions
  double X_0 = 20;
  double X_DOT_0 =  0;
  double X_DOTDOT_0;

  // Time Step and Time
  double delta_t = 0.05;
  double range = 100; //s
  int len = (int) range / delta_t;

  // Allocate memory for vector
  double **array = malloc(sizeof(double*) * 3);
  for (i = 0; i < 3; i++)
  {
    array[i] = malloc(sizeof(double) * len);
  }

  // Initialization of initial conditions
  array[0][0] = X_0;
  array[1][0] = X_DOT_0;
  array[2][0] = 0.0;

  double xdot;
  double xdotdot;

  for (j = 1; j < len; j++)
  {

    array[0][j] = array[0][j-1] + (delta_t * array[1][j-1]);

    xdotdot = getxdotdot(array[2][j-1], array[0][j-1], array[1][j-1]);
    array[1][j] = array[1][j-1] + (delta_t * xdotdot);

    array[2][j] = array[2][j-1] + delta_t;
  }

  printf("Linear Oscillator (gamma = %f,omega = %f, c = %f)\n", gamma, omega, c);
  printf("Time\tx\tx'\t\n");
  for (j = 0; j < len; j++)
  {
    printf("%f\t%f\t%f\n", array[2][j], array[0][j], array[1][j]);
  }

}
