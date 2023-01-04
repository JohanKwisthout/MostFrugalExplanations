#ifndef MFESIMHEADER 
#define MFESIMHEADER

// STL includes
#include <iostream>
#include <string>
#include <numeric>
#include <set>
#include <vector>
#include <map>
#include <cstdlib>
#include <experimental/random>
#include "dai/alldai.h"  		// Include main libDAI header file
#include "dai/factorgraph.h"
#include "dai/jtree.h"
#include "dai/varset.h"
#include "dai/index.h"

// comment for production mode, uncomment for debug messages
#define DEBUGMODE

#ifdef DEBUGMODE
	#define DEBUG(a) a;
#else
	#define DEBUG(a) ;
#endif	

double relevance(dai::FactorGraph fg, unsigned int node, std::vector<unsigned int> evidence_vars, std::vector<unsigned int> evidence_values, 
	std::vector<unsigned int> hypothesis_vars, std::vector<unsigned int> intermediate_vars, unsigned long int samples, std::mt19937 rngen);

std::vector<unsigned long int> compute_MFE(dai::FactorGraph fg, std::vector<unsigned int> evidenceVars, std::vector<unsigned int> evidenceValues,
	std::vector<unsigned int> hypothesisVars, std::vector<unsigned int> relevantVars, std::vector<unsigned int> irrelevantVars,
	bool relevanceComputation, unsigned long int samplesRel, double relThreshold, unsigned long int samples, unsigned long int cutoffTime);

std::vector<unsigned long int> get_mpe(dai::FactorGraph fg, std::vector<unsigned int> evidence_vars, std::vector<unsigned int> evidence_values);
std::vector<unsigned long int> get_map(dai::FactorGraph fg, std::vector<unsigned int> hypothesis_vars, std::vector<unsigned int> evidence_vars,
	std::vector<unsigned int> evidence_values, bool mapList);
std::vector<unsigned long int> prior_map(dai::FactorGraph fg, std::vector<unsigned int> hypothesis_vars);
std::vector<unsigned long int> local_prior_map(dai::FactorGraph fg, std::vector<unsigned int> hypothesis_vars, 
    std::vector<double>& map_scores);
std::vector<unsigned int> getIntermediateVars(dai::FactorGraph fg, std::vector<unsigned int> hypothesis_vars, 
    std::vector<unsigned int> evidence_vars);

std::vector<unsigned long int> annealed_map(dai::FactorGraph fg, std::vector<unsigned int> hypothesis_vars, std::vector<unsigned int> evidence_vars,
	std::vector<unsigned int> evidence_values, unsigned long int cutoffTime);

bool weak_map_indep(dai::FactorGraph fg, std::vector<unsigned int> evidenceVars, std::vector<unsigned int> evidenceValues, 
    std::vector<unsigned int> hypothesisVars, std::vector<unsigned int> hypothesisValues, std::vector<unsigned int> independenceTestVars, unsigned long int cutoffTime);
bool strong_map_indep(dai::FactorGraph fg, std::vector<unsigned int> evidenceVars, std::vector<unsigned int> evidenceValues, 
    std::vector<unsigned int> hypothesisVars, std::vector<unsigned int> hypothesisValues, std::vector<unsigned int> independenceTestVars, unsigned long int cutoffTime);
std::vector<unsigned long int> max_weak_map_indep(dai::FactorGraph fg, std::vector<unsigned int> evidenceVars, std::vector<unsigned int> evidenceValues, 
    std::vector<unsigned int> hypothesisVars, std::vector<unsigned int> hypothesisValues, std::vector<unsigned int> independenceTestVars, unsigned long int cutoffTime);
std::vector<unsigned long int> max_strong_map_indep(dai::FactorGraph fg, std::vector<unsigned int> evidenceVars, std::vector<unsigned int> evidenceValues, 
    std::vector<unsigned int> hypothesisVars, std::vector<unsigned int> hypothesisValues, std::vector<unsigned int> independenceTestVars, unsigned long int cutoffTime);
double weak_map_indep_measure(dai::FactorGraph fg, std::vector<unsigned int> evidenceVars, std::vector<unsigned int> evidenceValues, 
    std::vector<unsigned int> hypothesisVars, std::vector<unsigned int> hypothesisValues, std::vector<unsigned int> independenceTestVars, 
    unsigned long int cutoffTime, bool decision);
double weak_map_indep_measure(dai::FactorGraph fg, std::vector<unsigned int> evidenceVars, std::vector<unsigned int> evidenceValues, 
    std::vector<unsigned int> hypothesisVars, std::vector<unsigned int> hypothesisValues, std::vector<unsigned int> independenceTestVars, 
    unsigned long int cutoffTime, bool decision);
double strong_map_indep_measure(dai::FactorGraph fg, std::vector<unsigned int> evidenceVars, std::vector<unsigned int> evidenceValues, 
    std::vector<unsigned int> hypothesisVars, std::vector<unsigned int> hypothesisValues, std::vector<unsigned int> independenceTestVars, 
    unsigned long int cutoffTime, bool decision);


std::ostream& operator<<(std::ostream& os, const std::vector<int> &input);
std::ostream& operator<<(std::ostream& os, const std::vector<unsigned int> &input);
std::ostream& operator<<(std::ostream& os, const std::vector<long unsigned int> &input);

void iterate(unsigned int dimensions, unsigned int skip_node, std::vector<unsigned int> &ordinates, std::vector<unsigned int> maximums);
void random_sample(unsigned int dimensions, unsigned int skip_node, std::vector<unsigned int> &ordinates, std::vector<unsigned int> maximums,
	 std::mt19937 rngen);

int sample(dai::Factor fact, double rand);
double CalculateSpecHeat(const std::vector<double> &scores, const double &temperature, const double &bestScore);

#endif // defined MFESIM
