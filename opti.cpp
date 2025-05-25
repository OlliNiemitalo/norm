#include "opti.hpp"

#include <math.h>
#include <time.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <float.h>
#include <assert.h>

namespace Opti {

  MTRand rng;
  
  // Perform Fisher-Yates shuffle 
  void shuffle(int *table, int num) {
    for (int t = 0; (t < num); t++) {
      int u = rng.randInt(num-t-1);
      int temp = table[t];
      table[t] = table[t+u];
      table[t+u] = temp;
    }
  }
  
  // Perform partial Fisher-Yates shuffle
  // Only first numshuffle entries in the table are shuffled properly with the rest of the table
  void partialShuffle(int *table, int numtotal, int numshuffle) {
    for (int t = 0; (t < numshuffle); t++) {
      int u = rng.randInt(numtotal-t-1);
      int temp = table[t];
      table[t] = table[t+u];
      table[t+u] = temp;
    }
  }
  
  // Compute square of the perpendicular (that is, shortest) distance from 
  // a point (point) to a line in a multidimensional space. The line is 
  // defined as pointonline+a*linedirection where a is a scalar and 
  // pointonline and linedirection are vectors. numdimensions is the number 
  // of dimensions.
  double squaredPerpendicularDistance(double const *pointonline, double const *linedirection, double const *point, int numdimensions)
  {
    double s2=0;
    double b=0;
    double v2=0;
    for(int i=0;i<numdimensions;i++) {
      double d=point[i]-pointonline[i];
      s2+=d*d;
      b+=linedirection[i]*d;
      v2+=linedirection[i]*linedirection[i];
    }
    return s2-b*b/v2;
  }
  
  // Randomize by system clock. Please call me!
  void randomize() {
    srand((unsigned)time(0));
  }

  /*  
  // Generate gaussian random number.
  // mean = 0, standard deviation = 1
  double normalRandom()
  {
  double R1=rand()*(1.0/(RAND_MAX+1.0)); 
  double R2=rand()*(1.0/(RAND_MAX+1.0));
  double result = sqrt((-2)*log(1-R1))*cos((2*3.1415926535897932384626433832795)*R2);
  return result;   
  }
  */

  // Print parameter vector.
  void Problem::print(double *params)
  {
    printf("%.17f",params[0]);
    for(int i=1;i<getNumDimensions();i++)
      {
	printf(",%.17f", params[i]);
      }
    printf("\n");
  }
  
  Problem::~Problem()
  {
  }
  
  Strategy::~Strategy()
  {
  }
  
  Recombinator::~Recombinator()
  {
  }
  
  
  // The PCX recombinator
    
  PCXRecombinator::PCXRecombinator(int numparents, double sd1, double sd2)
  {
    assert(numparents>0);
    this->numparents = numparents;
    meanvector = 0;
    this->sd1 = sd1;
    this->sd2 = sd2;
  }
  
  void PCXRecombinator::setNumDimensions(int numdimensions)
  {
    assert(numdimensions>0);
    this->numdimensions = numdimensions;
    meanvector = new double[numdimensions];
  }
  
  int PCXRecombinator::numParents()
  {
    return numparents;
  }
  
  PCXRecombinator::~PCXRecombinator()
  {
    delete[] meanvector;
  }
  
  void PCXRecombinator::recombine(double *dest, double const *const *parents)
  {
    int u,t;                  // variables used in for loops
        
    assert(meanvector);
      
    // parents[0] is the Chosen One.  
    // 1. Calculate vector from parents[0] to mean of parents[0..numparents-1]
    for (u = 0; (u < numdimensions); u++) {
      meanvector[u] = parents[0][u];
    }
    for (t = 1; (t < numparents); t++) {
      for (int u = 0; (u < numdimensions); u++) {
	meanvector[u] += parents[t][u];
      }
    }
    double meanvectorlengthsquared = 0;
    for (u = 0; (u < numdimensions); u++) {
      meanvector[u] *= (1.0/numparents);
      meanvector[u] -= parents[0][u];
      meanvectorlengthsquared += meanvector[u]*meanvector[u];
    }
    // 2. Calculate mean of perpendicular distances from parents to mean vector
    double meansquareddistance = 0;
    for (t = 1; (t < numparents); t++) {
      meansquareddistance += squaredPerpendicularDistance(parents[0], meanvector, parents[t], numdimensions);
    }
    meansquareddistance /= numparents-1;
      
    double rmsdistance = sqrt(meansquareddistance);
      
    // 3. Random each dimension
    if (meanvectorlengthsquared == 0) {
      for (u = 0; (u < numdimensions); u++) {
	dest[u] = rng.randNorm(parents[0][u], sd2*rmsdistance);
      }
    } else {
      double dotproduct = 0;
      double meanvectorlength = sqrt(meanvectorlengthsquared);
      for (u = 0; (u < numdimensions); u++) {		
	dest[u] = rng.randNorm(0, 1);
	dotproduct += dest[u]*meanvector[u];
      }
      for (u = 0; (u < numdimensions); u++) {
	double along = meanvector[u]*(dotproduct/meanvectorlengthsquared);
	dest[u] = parents[0][u] + along*(meanvectorlength*sd1) + 
	  (dest[u]-along)*(rmsdistance*sd2);
      }
    }
  }
  
  
  
  G3::Individual::Individual()
  {
    cost=0;
    vector=0;
  }
  
  G3::Individual::~Individual()
  {
    delete [] vector;
  }
  
  void G3::Individual::swap(Individual &other)
  {
    double *temp=vector;
    double temp2=cost;
    vector=other.vector;
    cost=other.cost;
    other.vector=temp;
    other.cost=temp2;
  }
  
  void G3::Individual::init(int dim)
  {
    cost=0;
    vector=new double[dim];
  }
  
  G3::G3(Problem *problem, int populationsize, Recombinator *recombinator, int numOffspring)
  {
    this->recombinator=recombinator;
    recombinator->setNumDimensions(problem->getNumDimensions());
    this->numParents=recombinator->numParents();
     
    this->problem=problem;
    this->populationsize=populationsize;
    this->numOffspring=numOffspring;
    this->numdimensions=problem->getNumDimensions();
     
    population=new Individual[populationsize];
     
    double *min=problem->getMin();
    double *max=problem->getMax();
     
    for(int i=0;i<populationsize;i++)
      {
	population[i].init(numdimensions);
	for(int j=0;j<numdimensions;j++)
	  {
	    population[i].vector[j]=min[j]+(max[j]-min[j])*rng.rand();
	  }
	population[i].cost=problem->costFunction(population[i].vector,DBL_MAX);
      }
    offspring.init(numdimensions);
     
    parentList=new double *[numParents+2];
  }
  
  G3::~G3()
  {
    delete recombinator;
    delete [] population;
    delete [] parentList;
  }
  
  double *G3::best()
  {
    return population[0].vector;
  }
  
  double G3::averageCost()
  {
    double sum=0;
    for(int i=0;i<populationsize;i++)
      {
	sum+=population[i].cost;
      }
    return sum/populationsize;
  }
  
  double G3::evolve()
  {
    int i;
    // population[0] contains best
        
    for(i=1;i<numParents+2;i++)
      {
	int j=i+rng.randInt(populationsize-i-1);
	population[i].swap(population[j]);
	parentList[i]=population[i].vector;
      }
    parentList[0]=population[0].vector;
    std::swap(parentList[0],parentList[rng.randInt(numParents-1)]);
      
    Individual *best=&population[numParents+0];
    Individual *nextBest=&population[numParents+1];
    if(nextBest->cost < best->cost)
      std::swap(best,nextBest);
      
    for(i=0;i<numOffspring;i++)
      {
	recombinator->recombine(offspring.vector,parentList);
	offspring.cost=problem->costFunction(offspring.vector,nextBest->cost);
	if(offspring.cost<nextBest->cost)
	  {
	    nextBest->swap(offspring);
	    if(nextBest->cost < best->cost)
	      std::swap(best,nextBest);
	  }
      }
    if(best->cost < population[0].cost)
      best->swap(population[0]);
      
    return population[0].cost;
  }
  
  // Get latest best parameter vector in population
  double *DE::best()
  {
    return _best;		
  }
  
  double DE::averageCost()
  {
    return sumcost/np;
  }
  
  // Find best parameter vector in population, and calculate sum of costs (for average cost)
  void DE::statistics()
  {
    _best = &population[0]; // Just in case
    sumcost = 0;
    bestcost = DBL_MAX;
    for (int t = 0; (t < np); t++) {
      costs[t] = problem->costFunction(&population[t*d], DBL_MAX);
      sumcost += costs[t];
      if (costs[t] < bestcost) {
	bestcost = costs[t];
	_best = &population[t*d];
      }
    }        
  }
  
  // Non-constant member functions
  // -----------------------------
      
  // Fill parameters in all population members with random values (restarts evolution)
  //
  // minx[] = parameter vector with all parameters set to minimum possible values
  // maxx[] = parameter vector with all parameters set to maximum possible values
  void DE::randomPopulation(double *minx, double *maxx)
  {
    for (int member = 0; (member < np); member++) {
      for (int param = 0; (param < d); param++) {
	population[member*d+param] = rng.rand(maxx[param]-minx[param])+minx[param];
      }
    }
    _best = NULL;
  }
  
  double DE::evolve()
  {
    // The first of the parents is the destination vector
    // (so that DERecombinator can do crossover)
    parents[0] = &population[pos*d];
    // Pick additional numparents-1 parents at random
    partialShuffle(permuter, np, numparents-1);
    for (int t = 1; t < numparents; t++) {
      parents[t] = &population[permuter[t+1]*d];
    }
    // Recombine parents into trialvector
    recombinator->recombine(trialvector, parents);
    // Get cost of trialvector
    double trialcost = problem->costFunction(trialvector, costs[pos]);
    // If better than destination vector, replace it
    if (trialcost < costs[pos]) {
      memcpy(&population[pos*d], trialvector, sizeof(double)*d);
      // Update sumcost and costs[] and possibly _best and bestcost
      sumcost -= costs[pos];
      costs[pos] = trialcost;
      sumcost += trialcost;
      if (trialcost < bestcost) {
	bestcost = trialcost;
	_best = &population[pos*d];
      }
    }
    // Update gencost (sum of costs from 0..pos)
    gencost += costs[pos];
    // If we have browsed through the whole population...
    if (++pos >= np) {
      pos = 0;
      // Reset sumcost to a stable value
      // to avoid drift in floating point accumulation
      sumcost = gencost;
      // Restart gencost
      gencost = 0;
    }
    return bestcost;
  }

  DE::DE(Problem *problem, int np, Recombinator *recombinator)
  {
    this->problem=problem;
    this->np = np;
    this->d = problem->getNumDimensions();
    recombinator->setNumDimensions(d);
    this->recombinator = recombinator;
    this->trialvector = new double[d];
    numparents = recombinator->numParents();
    parents = new double *[numparents];
    population = new double[np*d];
    costs = new double[np];
    permuter = new int[np];
    for (int member = 0; (member < np); member++) {
      permuter[member] = member;
    }
    randomPopulation(problem->getMin(), problem->getMax()); // Initialize population
    pos = 0; // Point at first parent
    gencost = 0;
    statistics();
  }
  
  // Destructor
  DE::~DE() {
    delete[] population;
    delete[] costs;
    delete[] permuter;
    delete[] parents;
    delete[] trialvector;
  }

  // Differential Evolution recombinator
  DERecombinator::DERecombinator(double cr, double c) {
    this->c = c;
    this->cr = cr;
  }

  void DERecombinator::setNumDimensions(int numDimensions) {
    d = numDimensions;
  }

  int DERecombinator::numParents() {
    return 4; // target, parent1+(parent2-parent3)
  }

  void DERecombinator::recombine(double *dest, double const *const *parents) {
    // Start at a random parameter
    int pos = rng.randInt(d-1);
    for (int count = 0; count < d;) {
      dest[pos] = parents[1][pos] + c*(parents[2][pos] - parents[3][pos]);
      if (++pos >= d) pos = 0;
      count++;
      if (rng.randExc() > cr) {
	for (;count < d; count++) {
	  dest[pos] = parents[0][pos];
	  if (++pos >= d) pos = 0;
	}
      }
    }
  };

   
}
