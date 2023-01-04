/************************************************************************/
/* MAP Independence algorithm implementation   					        */
/* Written by:		Johan Kwisthout                                		*/
/* Version:			1.0													*/
/* Last changed:	04-01-2023                                         	*/
/*                                                                     	*/
/* Version History:                                                    	*/
/*                                                                     	*/
/* Version Comments:                                                   	*/
/* - implementation of the algorithm described in Kwisthout (2021)     	*/
/* - cutoffTime is not yet used in the algorithms						*/
/************************************************************************/

// headers
#include "mfesim.h"
#include "combinations.hpp"		// Howard Hinnant's combinations template
#include <chrono>

bool weak_map_indep(dai::FactorGraph fg, std::vector<unsigned int> evidenceVars, std::vector<unsigned int> evidenceValues, 
    std::vector<unsigned int> hypothesisVars, std::vector<unsigned int> hypothesisValues, std::vector<unsigned int> independenceTestVars, unsigned long int cutoffTime)
{
    if (weak_map_indep_measure(fg, evidenceVars, evidenceValues, hypothesisVars, hypothesisValues, independenceTestVars, cutoffTime, true) == 1.0)
        return true;
    else
        return false;        
}

double weak_map_indep_measure(dai::FactorGraph fg, std::vector<unsigned int> evidenceVars, std::vector<unsigned int> evidenceValues, 
    std::vector<unsigned int> hypothesisVars, std::vector<unsigned int> hypothesisValues, std::vector<unsigned int> independenceTestVars, 
    unsigned long int cutoffTime, bool decision)
{
    int count = 0, different = 0;
    std::vector<unsigned long int> map;
    for (const unsigned int &e: hypothesisValues) { map.push_back((unsigned int) e); }

    // for each variable R in independenceTestVars
    for (auto varR = independenceTestVars.begin(); varR != independenceTestVars.end(); ++varR)
	{
        // for each value r of R
        for (unsigned int state = 0; state < fg.var(*varR).states(); state++)
        {
            // local copy of evidence nodes
            std::vector<unsigned int> mapTestVars(evidenceVars);
            mapTestVars.push_back(*varR);
            // set the evidence to this state
            std::vector<unsigned int> mapTestValues(evidenceValues);
            mapTestValues.push_back(state);
            // get the map
            std::vector<unsigned long int> best = get_map(fg, hypothesisVars, mapTestVars, mapTestValues, false);
            if (best == map)
            {
                DEBUG(std::cout << "Same for R = " << *varR << " and r = " << state << std::endl;)
            }
            else
            {
                DEBUG(std::cout << "Different for R = " << *varR << " and r = " << state << std::endl;)
                different++;
                if (decision) return 1.0;
            }
            count++;                
        }
    }
    DEBUG(std::cout << "Quantified weak MAP independence:  " << ((double) different / (double) count) << std::endl;)
    return ((double) different / (double) count);
}

std::vector<unsigned long int> max_weak_map_indep(dai::FactorGraph fg, std::vector<unsigned int> evidenceVars, std::vector<unsigned int> evidenceValues, 
    std::vector<unsigned int> hypothesisVars, std::vector<unsigned int> hypothesisValues, std::vector<unsigned int> independenceTestVars, unsigned long int cutoffTime)
{
	// we simply test for each of the variables in independenceTestVars whether they are
	// weakly map independent and if so, we add them to the set 'weak'

	std::vector<unsigned long int> weak;

	for (auto varR = independenceTestVars.begin(); varR != independenceTestVars.end(); ++varR)
	{
		std::vector<unsigned int> varVec(1, *varR);
		if (weak_map_indep(fg, evidenceVars, evidenceValues, hypothesisVars, hypothesisValues, varVec, cutoffTime))
			weak.push_back(*varR);
	}
    return weak;
}

bool strong_map_indep(dai::FactorGraph fg, std::vector<unsigned int> evidenceVars, std::vector<unsigned int> evidenceValues, 
    std::vector<unsigned int> hypothesisVars, std::vector<unsigned int> hypothesisValues, std::vector<unsigned int> independenceTestVars, unsigned long int cutoffTime)
{
    if (strong_map_indep_measure(fg, evidenceVars, evidenceValues, hypothesisVars, hypothesisValues, independenceTestVars, cutoffTime, true) == 1.0)
        return true;
    else
        return false;        
}

double strong_map_indep_measure(dai::FactorGraph fg, std::vector<unsigned int> evidenceVars, std::vector<unsigned int> evidenceValues, 
    std::vector<unsigned int> hypothesisVars, std::vector<unsigned int> hypothesisValues, std::vector<unsigned int> independenceTestVars, 
    unsigned long int cutoffTime, bool decision)
{
    int count = 0, different = 0, nr_vars = 0;
    unsigned long int iteration = 0, max_iterations = 1;

    // actual map (h*) and running MAP values
    std::vector<unsigned long int> map, best;
    for (const unsigned int &e: hypothesisValues) { map.push_back((unsigned int) e); }

    // copy of all evidence and indep test variables
    std::vector<unsigned int> mapTestVars;

    // iterator over the joint value assignments to R
    std::vector<unsigned int> independenceValues;

    // maximum values of these joint value assignments
    std::vector<unsigned int> independenceMaxValues;

    // initialize these lists
    for (auto inter: independenceTestVars)
    {
        mapTestVars.push_back(inter);

        unsigned int st = fg.var(inter).states();
         
        independenceValues.push_back(0);
        independenceMaxValues.push_back(st - 1);
        max_iterations *= st;
        nr_vars++;
    }
    DEBUG(std::cout << "Testing R = " << mapTestVars << std::endl;)

    // vector of all test+evidence variables (add evidence here here)
    std::copy(evidenceVars.begin(), evidenceVars.end(), back_inserter(mapTestVars));

    // for each joint value assignment over independenceTestVars
    for (iteration = 1; iteration <= max_iterations; iteration++)
    {
        // collate evidence (add evidence + jva)
        std::vector<unsigned int> mapTestValues;
        std::copy(independenceValues.begin(), independenceValues.end(), back_inserter(mapTestValues));
        std::copy(evidenceValues.begin(), evidenceValues.end(), back_inserter(mapTestValues));

        DEBUG(std::cout << "Testing " << mapTestVars << " with value " << mapTestValues << std::endl;)

        // find MAP for this value
        best = get_map(fg, hypothesisVars, mapTestVars, mapTestValues, false);

        if (best == map)
        {
            // (best == h*)
            DEBUG(std::cout << "Same for r = " << independenceValues  << std::endl;)
        }
        else
        {
            // (best != h*)
            DEBUG(std::cout << "Different for r = " << independenceValues << std::endl;)
            different++;
            if (decision) return 1.0;
        }
      	// next value in iteration
        iterate(nr_vars, -1, independenceValues, independenceMaxValues);
        count++;                
    }	
    DEBUG(std::cout << "Quantified strong MAP independence:  " << ((double) different / (double) count) << std::endl;)
    return ((double) different / (double) count);
}

std::vector<unsigned long int> max_strong_map_indep(dai::FactorGraph fg, std::vector<unsigned int> evidenceVars, std::vector<unsigned int> evidenceValues, 
    std::vector<unsigned int> hypothesisVars, std::vector<unsigned int> hypothesisValues, std::vector<unsigned int> independenceTestVars, unsigned long int cutoffTime)
{
	// This is a very time-consuming algorithm: we iterate over all subsets of independenceTestVars,
	// run strong_map_indep over this subset, and if it answers 'yes' we keep track of the largest size
	// we make use of Howard Hinnant's combinations and permutations; in particular, we use the code at
	// https://stackoverflow.com/questions/25984609/iteratively-calculate-the-power-set-of-a-set-or-vector

    std::vector<unsigned long int> strong;

	unsigned int max = 0;

    for (std::size_t k = 0; k <= independenceTestVars.size(); ++k)
	{
        for_each_combination(independenceTestVars.begin(), independenceTestVars.begin()+k, independenceTestVars.end(),
			[&](std::vector<unsigned int>::const_iterator first, std::vector<unsigned int>::const_iterator last)
        {
			std::vector<unsigned int> testVars (first, last);
			DEBUG(std::cout << "Testing set " << testVars << std::endl;)

			if (strong_map_indep(fg, evidenceVars, evidenceValues, hypothesisVars, hypothesisValues,
				testVars, cutoffTime))
			{
				DEBUG(std::cout << "This is now the largest set of size " << k << std::endl;)
    
				max = k;
				strong.clear();
				std::copy(testVars.begin(), testVars.end(), back_inserter(strong));
				return true;		// we have a set of size k, so try k+1
			}
			else
			{
				return false;		// keep iterating combinatons
			}
        });	// end of for_each_combination template
	}

    return strong;
}


