/************************************************************************/
/* MFE algorithm implementation              					        */
/* Written by:		Johan Kwisthout                                		*/
/* Version:			1.0                 								*/
/* Last changed:	26-05-2022                                         	*/
/*                                                                     	*/
/* Version History:                                                    	*/
/*                                                                     	*/
/* Version Comments:                                                   	*/
/* - implementation of the algorithm described in Kwisthout (2015)     	*/
/************************************************************************/

// headers
#include "mfesim.h"
#include <chrono>

std::vector<unsigned long int> compute_MFE(dai::FactorGraph fg, std::vector<unsigned int> evidenceVars, std::vector<unsigned int> evidenceValues,
	std::vector<unsigned int> hypothesisVars, std::vector<unsigned int> relevantVars, std::vector<unsigned int> irrelevantVars,
	bool relevanceComputation, unsigned long int samplesRel, double relThreshold, unsigned long int samples, unsigned long int cutoffTime)
{
    unsigned long int timeBound = cutoffTime * 1000000000UL;

	std::vector<unsigned long int> MFE;

   // maximum values of these variables
    std::vector<unsigned int> irrelevant_max_values;

	// sample over irrelevant variables
    std::vector<unsigned int> irrelevant_sample;

	// random number generator
	std::random_device rd;                      // note C++11 code
    std::mt19937 gen(rd());						// random numbers by Mersenne twister algorithm

	// if relevanceComputation is true, we need to populate the relevant intermediate variables (all is currently in irrelevant, we rebuild them)
	if (relevanceComputation)
	{
		std::vector<unsigned int> intermediateVars(irrelevantVars);
		irrelevantVars.clear();
		
	    for (auto inter = intermediateVars.begin(); inter != intermediateVars.end(); ++inter)
	    {
			double rel = relevance(fg, *inter, evidenceVars, evidenceValues, hypothesisVars, intermediateVars, samplesRel, gen);
            DEBUG(std::cout << "relevance of " << *inter << " is " << rel << std::endl;)
			
			if (rel >= relThreshold)
			{
				// move *inter to relevantVars
				relevantVars.push_back(*inter);
			}
			else
			{
				// move *inter to irrelevantVars
				irrelevantVars.push_back(*inter);
			}
		}
	}

    // set maximum values per variable
    for (auto inter: irrelevantVars)
    {
		unsigned int st = fg.var(inter).states();
		irrelevant_sample.push_back(0);
        irrelevant_max_values.push_back(st - 1);
    }

	// get the current MAP
	std::vector<unsigned long int> map;
	std::map<std::vector<unsigned long int>, int> map_counts;
	std::map<std::vector<unsigned long int>, int>::iterator map_it;
	
	// MAIN loop (comment lines match the algorithm description):

    // internal time keeping to cut off computation after time bound
    auto start = std::chrono::steady_clock::now();
		
	// for n = 1 to N do
	for (unsigned long int n = 0; n < samples; n++)
	{
		// Choose i \in I- at random
        random_sample(irrelevantVars.size(), -1, irrelevant_sample, irrelevant_max_values, gen);
		
		std::vector<unsigned int> combined_evidence;
		std::vector<unsigned int> combined_evidence_values;
		std::copy(evidenceVars.begin(), evidenceVars.end(), std::back_inserter(combined_evidence));
		std::copy(irrelevantVars.begin(), irrelevantVars.end(), std::back_inserter(combined_evidence));
		std::copy(evidenceValues.begin(), evidenceValues.end(), std::back_inserter(combined_evidence_values));
		std::copy(irrelevant_sample.begin(), irrelevant_sample.end(), std::back_inserter(combined_evidence_values));


		// Determine h = argmax_h Pr(H = h, i, e)
		map = get_map(fg, hypothesisVars, combined_evidence, combined_evidence_values, false);
		
		// Collate the joint value assignments h (std::map<<vector>,int>) -- if <vector> does not exist, add it (int = 1) otherwise int++
		map_it = map_counts.find(map);
		
		if (map_it == map_counts.end()) 
		{
			map_counts.insert(std::pair<std::vector<unsigned long int>, int>(map, 1));
		}
		else
		{
			map_it->second++;
		}

        if (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count() > timeBound)
        {
            std::cout << "stopping computation - time bound" << std::endl;
            n = samples;
        }
	}
	
	// Decide upon the joint value assignment hmaj that was picked most often
	int mfe_max = 0;
	for (map_it = map_counts.begin(); map_it != map_counts.end(); ++map_it)
	{
		DEBUG(std::cout << "assignment " << map_it->first << " count was " << map_it->second << std::endl;)
		if (map_it->second > mfe_max)
		{
			mfe_max = map_it->second;
			MFE = map_it->first;
		}
	}

	return MFE;
}

