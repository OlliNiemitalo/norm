#include <stdio.h>
#include <math.h>
#include "keyboard.h"
#include "opti.hpp"

class NormProblem : public Opti::Problem {
private:
  int numParams;
  int numSamples;
  double *min;
  double *max;
  double *x;
  double *y;

public:

  double ramp;

  // numParams must be odd
  NormProblem(int numParams, int numSamples, double startX, double endX, double ramp, double* candidate = NULL) : numParams(numParams), numSamples(numSamples), ramp(ramp) {
    min = new double[numParams];
    max = new double[numParams];
    x = new double[numSamples];
    if (candidate != NULL) {
      for (int i = 0; i < numParams; i++) {
        min[i] = candidate[i]-fabs(candidate[i])*(1.0/65536);
        max[i] = candidate[i]+fabs(candidate[i])*(1.0/65536);
      }
    } else {
      for (int i = 0; i < numParams; i++) {
        min[i] = -0.5;
        max[i] = 0.5;
      }
    }
    for (int i = 0; i < numSamples; i++) {
      // x[i] = startX + (endX-startX)*i/(numSamples-1);  // Uniform sampling
      x[i] = startX + (endX-startX)*(0.5 - 0.5*cos(M_PI*i/(numSamples-1)));  // Like Chebyshev nodes but including end points, to make MSE-optimimal similar to maxabs-error-optimal
    }
  }

  double *getMin() {
    return min;
  }

  double *getMax() {
    return max;
  }

  void print(double *params) {
    printf("Printout:\n");
    for (int i = 0; i < numParams; i += 3) {
      printf("(");
      for (int j = 0; j < 3; j++) {
        printf("%.20f", params[i+j]);
        if (j < 2) {
          printf(", ");
        }
      }
      printf("),\n");
    }
    printf("\n");
    for (int i = 0; i < numParams; i += 3) {
      for (int j = 0; j < 3; j++) {
        printf("%.20f x^%d + ", params[i+j], j*2+1);
      }
      printf("\n");
    }
    printf("\n");
  }

  double costFunction(double *params, double compare) {
    double maxAbsErr = 0.0;
    double accu = 0.0;
    params[0] = fabs(params[0]);
    for (int j = 0; j < numParams; j += 3) {
      params[j] = params[0];
    }
    for (int i = 0; i < numSamples; i++) {
      double y = x[i];
      for (int j = 0; j < numParams; j += 3) {
        y = params[j]*y + params[j+1]*(y*(y*y)) + params[j+2]*(y*(y*y)*(y*y));
      }
      accu += (y - 1.0)*(y - 1.0);
      double absErr = fabs(y - 1.0);
      if (absErr > maxAbsErr) {
        maxAbsErr = absErr;
      }
    }
    double rms = sqrt(accu/numSamples);
    double retval = rms + ramp*(maxAbsErr - rms);
    return retval;
  } 

  int getNumDimensions() {
    return numParams;
  }

  ~NormProblem() {
    delete[] min;
    delete[] max;
    delete[] x;
  }
};

int main() {
  INITKEYBOARD;
  int rampStart =  500000;
  int rampEnd   = 1000000;
  NormProblem problem(3*5, 8193, 0.001, 1.010916328, 0.0);
  Opti::DERecombinator deRecombinator(0.999, 0.76);
  Opti::DE optimizer(&problem, 1000, &deRecombinator);
  for(int t = 0;; t++) {
    if (optimizer.pos == 0) {
      if (t > rampStart && t <= rampEnd) {
        problem.ramp = t < rampEnd ? double(t-rampStart)/(rampEnd-rampStart) : 1.0; // Ramp between root mean square error and max abs error
        optimizer.statistics();  // Recalculate costs for the population because we changed ramp
      }
    }
    double bestcost = optimizer.evolve();
    if (!(t % 10000)) {
      printf("gen=%d, bestcost=%.20f, average=%.20f, ramp=%.5f\n", t, bestcost, optimizer.averageCost(), problem.ramp);
      if (kbhit()) {
        printf("Parameter vector printout:\n");
        problem.print(optimizer.best());
        printf("Best cost %f\n", problem.costFunction(optimizer.best(), 0));
        if (getch() == 27) {
          break;
        }
        getch();
      }
    }
  }
  DEINITKEYBOARD;
  return 0;
}

// Compile with:
// g++ -g -O0 optimize.cpp opti.cpp
// g++ optimize.cpp opti.cpp -ffast-math -march=native -O3 -Wno-unused-result