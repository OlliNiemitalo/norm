#ifndef OPTI_HPP
#define OPTI_HPP

// v1.1
// EVOLUTIONARY ALGORITHMS FOR THE OPTIMIZATION OF MULTIPLE REAL VARIABLES
// by minimization of an arbitrary function of those variables. Global minimum
// (perfect solution) cannot be guaranteed, but might be reached. The used 
// algorithms are outlined in [1] (DE) and [2] (G3PCX). G3PCX is a bit buggy
// and may sometimes give nans in the parameter vector.
// 
// Written in 2002-2003 by Olli Niemitalo (o@iki.fi) and Magnus Jonsson,
// and in 2019 by Olli Niemitalo.
// This work is placed in the public domain / released under CC0.
//
// Version history:
// v1.1, 2019-06-05
//      * Removed experimental optimizer GreedyMagnus and its recombinator
//      * Increased precision in parameter vector printout
// v1.0, 2003 Initial version
//
// References:
// 
// [1] Storn, R. and Price, K., "Differential Evolution - a Simple and 
// Efficient Adaptive Scheme for Global Optimization over Continuous Spaces", 
// Technical Report TR-95-012, ICSI, March 1995, ftp.icsi.berkeley.edu.
//
// [2] Deb, K , Anand, A., and Joshi, D (April, 2002). 
// A Computationally Efficient Evolutionary Algorithm for Real-Parameter 
// Optimization. KanGAL Report No. 2002003.

#include "MersenneTwister.h"

namespace Opti {

  // Mersenne twister random generator
  extern MTRand rng;
    
  // Perform Fisher-Yates shuffle 
  void shuffle(int *table, int num);
        
  // Perform partial Fisher-Yates shuffle 
  // Only first numshuffle entries in the table are shuffled properly with the rest of the table
  void partialShuffle(int *table, int numtotal, int numshuffle);

        
  // Compute square of the perpendicular (that is, shortest) distance from 
  // a point (point) to a line in a multidimensional space. The line is 
  // defined as pointonline+a*linedirection where a is a scalar and 
  // pointonline and linedirection are vectors. numdimensions is the number 
  // of dimensions.
  double squaredPerpendicularDistance(double const *pointonline, double const *linedirection, double const *point, int numdimensions);
    

  // Generate gaussian random number.
  // mean = 0, standard deviation = 1
  double normalRandom();
    

  // Optimization problem base class. This class should be inherited by a class
  // representing an actual optimization problem.
  class Problem {
  public:
    // Destructor
    virtual ~Problem();
		
    // Return number of parameters to optimize
    //
    // Must be implemented in the actual optimization problem.
    virtual int getNumDimensions()=0;
		
    // The parameters being optimized are assumed to be within a range.
    // These functions give arrays containing minimum and maximum values 
    // for the parameters. NOTE: Obtained solution may be outside these ranges,
    // but a solution will be found more easily if it is within these ranges.
    //
    // Must be implemented in the actual optimization problem.
    virtual double *getMin()=0;
    virtual double *getMax()=0;
		
    // "cost function" being minimized. Return cost for the set of parameters
    // given in the params array. The compare value is provided as an eariler 
    // cost value the current cost value is compared to. If the current cost
    // evaluation is known to give a higher value than the compare value, the
    // evaluation can be interrupted prematurely to save computing time and
    // a value higher than or equal to the compare value should be returned.
    // The cost function is allowed to modify params to ensure parameter 
    // constraints, wraparound, etc.
    // 
    // Must be implemented in the actual optimization problem.
    virtual double costFunction(double *params, double compare)=0; 
		
    // Print parameter vector to stdout.
    virtual void print(double *params);
  };
    
	
  // Optimization algorithm base class. Inherited by the actual optimization
  // algorithms (implemented later in this file).
  class Strategy {
  public:
    // Destructor
    virtual ~Strategy();
		
    // Return best parameter vector so far
    virtual double *best()=0;
		
    // Return average cost of the population
    virtual double averageCost()=0;
		
    // Evolve some...
    virtual double evolve()=0;
  };
	
    
  // Recombinator operator base class. Used in evolutionary algorithms.
  // This is almost like sex! :-)
  class Recombinator {
  public:    
    // 1: Strategy sends the number of dimensions in problem
    virtual void setNumDimensions(int numDimensions) = 0;
    // 2: Strategy asks: How many parents does this recombinator require?
    virtual int numParents() = 0;
    // 3: Strategy uses recombinator to make offspring from parents.
    // (dest being one of the parents causes undefined behaviour)
    virtual void recombine(double *dest, double const *const *parents) = 0;
    // 4: Strategy destroys its recombinator.
    virtual ~Recombinator();
  };
    

  // The PCX recombinator
  class PCXRecombinator : public Recombinator {
  private:
    int numparents;
    int numdimensions;
    double sd1, sd2;
    double *meanvector;
  public:
    void setNumDimensions(int numDimensions);
    PCXRecombinator(int numparents=3, double sd1 = 0.1, double sd2 = 0.1);
    ~PCXRecombinator();
    int numParents();
    void recombine(double *dest, double const *const *parents);
  };
	
    
  // The G3 evolution strategy using the PCX recombinator.
  class G3 : public Strategy
  {
  public:
    G3(Problem *problem, int populationsize, Recombinator *recombinator = new PCXRecombinator(), int numOffspring=2);
    ~G3();

    double *best();
    double averageCost();
    double evolve();
  private:
    class Individual
    {
    public:
      Individual();
      ~Individual();
			
      void swap(Individual &other);
      void init(int dim);

      double cost;
      double *vector;
    };

    Problem *problem;
    int numdimensions;
		
    int populationsize;
    Individual *population;
		
    int numOffspring;
    Individual offspring;
		
    int numParents;
    double **parentList;
		
    Recombinator *recombinator;
  };
    
  // Differential Evolution recombinator
  class DERecombinator : public Recombinator{
  private:
    int d;
    double cr, c;
  public:    
    void setNumDimensions(int numDimensions);
    int numParents();
    void recombine(double *dest, double const *const *parents);

    // Constructor
    // cr = Cross-over amount. 0 is unreasonable.
    // c = Weight for difference of two parents
    DERecombinator(double cr = 1.0, double c = 0.61803398875);
  };

  // Differential Evolution class
  // ----------------------------
  //
  // Tries to search for the global minimum of a cost function which gets a 
  // parameter vector as an argument
  class DE : public Strategy
  {
  public:
		
    int d;                 // Number of parameters
		
    // Get latest best parameter vector in population
    double *best();

    double averageCost();
		
    // Find best parameter vector in population
    //

    // z  = Pointer to cost function being minimized (argument: parameter vector) 
    //
    // Returns: Cost of best. Updates member variables best, bestcost and costs
    double findBest();
		
    // Non-constant member functions
    // -----------------------------
		
    // Fill parameters in all population members with random values (restarts evolution)
    //
    // minx[] = parameter vector with all parameters set to minimum possible values
    // maxx[] = parameter vector with all parameters set to maximum possible values
    void randomPopulation(double *minx, double *maxx);
		
    void statistics();

    // Evolve into next generation - try evolve(&z, 0, 0.7, 1.0, 1);
    //
    // z            = Pointer to cost function to minimize (argument: parameter vector) 
    // gainbest     = Coefficient for best population member (try 1.2)
    // gainr3       = Coefficient for 3rd random population member (try 1)
    // (Note: Coefficient for original population member = 1 - sum of the above)
    // gaindiffr1r2 = Coefficient for difference of 1st and 2nd random population member (try 0.5)
    // cr       = Crossing-over amount, 0..1
    //
    // Returns: Cost of best parameter vector in population. 
    // Member variable best becomes a pointer to the best parameter vector in population.
		
    double evolve();
		
    void init(Problem *problem, int np, Recombinator *recombinator);

    // Constructor
    DE(Problem *problem, int np, Recombinator *recombinator);
		
    // Destructor
    ~DE();

    int pos; // Where we are going in population
  private:
    int numparents; // Number of parents taken by recombinator
    double gencost; // Cost of generation
    double *costs;             // Costs of parameter vectors in population
    int np;                      // Number of population members
    double *population;        // Parameter vectors in population, one-by-one
    double *_best;        // Pointer to best parameter vector in population
    double bestcost;     // Cost of the above (invalid if best == NULL)
    double sumcost;      // sum of all costs in population, for calculation of average cost
    int *permuter;  // Shuffled parent table
    Problem *problem;
    Recombinator *recombinator;
    double **parents; // Temporary parents table for recombinator
    double *trialvector; // Temporary trial vector
  };
    	    
} // end namespace Opti

#endif
