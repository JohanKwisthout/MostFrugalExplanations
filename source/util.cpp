/************************************************************************/
/* MFE generic tools and functions             					        */
/* Written by:		Johan Kwisthout                                		*/
/* Version:			0.1 (not finished yet)								*/
/* Last changed:	08-01-2020                                         	*/
/*                                                                     	*/
/* Version History:                                                    	*/
/*                                                                     	*/
/* Version Comments:                                                   	*/
/*                                                                     	*/
/************************************************************************/

// headers
#include <chrono>
#include <ctime>
#include "mfesim.h"
#include "dai/alldai.h"
#include "dai/jtree.h"

// from www.techiedelight.com/print-vector-cpp/
std::ostream& operator<<(std::ostream& os, const std::vector<int> &input)
{
	for (auto const& i: input) { os << i << " "; }
	return os;
}

std::ostream& operator<<(std::ostream& os, const std::vector<unsigned int> &input)
{
	for (auto const& i: input) { os << i << " "; }
	return os;
}

std::ostream& operator<<(std::ostream& os, const std::vector<long unsigned int> &input)
{
	for (auto const& i: input) { os << i << " "; }
	return os;
}

// adapted after https://stackoverflow.com/questions/26844032/fastest-way-to-iterate-over-n-dimensional-array-of-arbitrary-extents
void iterate(unsigned int dimensions, unsigned int skip_node, std::vector<unsigned int> &ordinates, std::vector<unsigned int> maximums)
{
    // iterate over dimensions in reverse...
    for (int dimension = dimensions - 1; dimension >= 0; dimension--)
    {
        // skip node itself (set skip_node to -1 to ignore this step)
        if ((skip_node >= 0) && (dimension == skip_node))
		{
			if (dimension > 0)
	            dimension--;
			else			// last node, we're done!
				return;
		}

        if (ordinates[dimension] < maximums[dimension])
        {
            // if this dimension can handle another increment... then done.
            ordinates[dimension]++;
            break;
        }

        // otherwise, reset this dimension and bubble up to the next dimension to take a look
        ordinates[dimension] = 0;
    }
}

void random_sample(unsigned int dimensions, unsigned int skip_node, std::vector<unsigned int> &ordinates, std::vector<unsigned int> maximums,
	std::mt19937 rngen)
{
	auto start = std::chrono::steady_clock::now();
   // iterate over dimensions in reverse...
    for (int dimension = dimensions - 1; dimension >= 0; dimension--)
    {
        // skip node itself (set skip_node to -1 to ignore this step)
        if ((skip_node >= 0) && (dimension == skip_node))
		{
			if (dimension > 0)
	            dimension--;
			else			// last node, we're done!
				return;
		}

		// get a random value within this dimension
		std::uniform_int_distribution<> dist(0, maximums[dimension]);
		ordinates[dimension] = dist(rngen);
	}
	auto end = std::chrono::steady_clock::now();
	std::cout << "Taking a sample " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << " ns" << std::endl;
	DEBUG(std::cout << "Taking a sample " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << " ns" << std::endl;)
}

std::vector<unsigned long int> get_mpe(dai::FactorGraph fg, std::vector<unsigned int> evidence_vars, std::vector<unsigned int> evidence_values)
{
	// returns the mpe, the joint value assignment to the hypothesis vars that has maximum posterior probability given the evidence

    std::vector<unsigned long int> mpe;

    for (int i = 0; i < evidence_vars.size(); i++)
    {
        fg.clamp(evidence_vars[i], evidence_values[i], false);
    }
    
    dai::PropertySet opts;
    dai::JTree jt = dai::JTree(fg, opts("updates",std::string("HUGIN"))("inference",std::string("MAXPROD")));
    jt.init();
    jt.run();
    mpe = jt.findMaximum();

    return mpe;
}

std::vector<unsigned long int> get_map(dai::FactorGraph fg, std::vector<unsigned int> hypothesis_vars, std::vector<unsigned int> evidence_vars,
	std::vector<unsigned int> evidence_values, bool mapList)
{
	// returns the map, the joint value assignment to the hypothesis vars that has maximum posterior probability given the evidence
	// while marginalizing over the (relevant) intermediate variables. As libDAI has no MAP function we just compute the distribution
	// over the MAP variables and select the state with maximum value from the posterior, which is the MAP assignment

	// when used in MFE function, the evidence is the actual 'real' evidence plus the sampled irrelevant intermediate nodes

    std::vector<unsigned long int> map;
    std::vector<unsigned long int> h_vars(begin(hypothesis_vars), end(hypothesis_vars));    // needs cast to long

	dai::VarSet hypSet = fg.inds2vars(h_vars);

	auto start = std::chrono::steady_clock::now();
    for (int i = 0; i < evidence_vars.size(); i++)
    {
        fg.clamp(evidence_vars[i], evidence_values[i], false);
    }
	auto end = std::chrono::steady_clock::now();
	DEBUG(std::cout << "Clamping evidence " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << " ns" << std::endl;)

	start = std::chrono::steady_clock::now();
    dai::PropertySet opts;
    dai::JTree jt = dai::JTree(fg, opts("updates",std::string("HUGIN"))("inference",std::string("SUMPROD")));
    jt.init();
    jt.run();
	end = std::chrono::steady_clock::now();
	DEBUG(std::cout << "JT run " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << " ns" << std::endl;)

	start = std::chrono::steady_clock::now();
	dai::Factor hypFact = jt.calcMarginal(hypSet);
	end = std::chrono::steady_clock::now();
	DEBUG(std::cout << "Marginal time " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << " ns" << std::endl;)
    
	// find element with maximum value ( = MAP explanation)
	start = std::chrono::steady_clock::now();
	double max = 0.0;
	int entry = 0; 
    for (int i = 0; i < hypFact.nrStates(); i++)
    {
        if (mapList)
        {
            std::cout << "entry ";
            for (auto const& j: dai::calcState(hypFact.vars(), i))
                std::cout << j.second;
            std::cout << " has probability " << hypFact.p()[i] << std::endl;
        }
		if (hypFact.p()[i] > max)
		{
    	    max = hypFact.p()[i];
			entry = i;
		}
    }

	// transform index to map of <Var, value> pairs
	std::map<dai::Var, size_t> mapValues = dai::calcState(hypFact.vars(), entry);
	end = std::chrono::steady_clock::now();
	DEBUG(std::cout << "MAP time " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << " ns" << std::endl;)

	// now set map accordingly to the values in hypothesis_vars
	for (auto const& i: mapValues)
	{
		map.push_back(i.second);
	}
    DEBUG(std::cout << "map " << map << " has probability " << max << std::endl;)

	return map;
}

std::vector<unsigned long int> prior_map(dai::FactorGraph fg, std::vector<unsigned int> hypothesis_vars)
{
	// returns the joint value assignment to the hypothesis vars with maximum *prior* probability
	std::vector<unsigned int> evidence;			// empty evidence
	std::vector<unsigned int> evidence_values;	// empty evidence
	return get_map(fg, hypothesis_vars, evidence, evidence_values, false);
}	

std::vector<unsigned long int> local_prior_map(dai::FactorGraph fg, std::vector<unsigned int> hypothesis_vars, 
    std::vector<double>& map_scores)
{
	// returns the assignments to the hypothesis vars which each individually have maximum *prior* probability
    std::vector<unsigned long int> local_map;

    dai::PropertySet opts;
    dai::JTree jt = dai::JTree(fg, opts("updates",std::string("HUGIN"))("inference",std::string("SUMPROD")));
    jt.init();
    jt.run();
	
	for (auto const& var: hypothesis_vars)
	{
		dai::Factor hypFact = jt.belief(fg.var(var));
		double max = 0.0;
		int entry = 0; 
        for (int i = 0; i < hypFact.nrStates(); i++)
		{
			if (hypFact.p()[i] > max)
			{
				max = hypFact.p()[i];
				entry = i;
			}
		}
        map_scores.push_back(max);
		local_map.push_back(entry);
	}
	return local_map;
}	

int sample(dai::Factor fact, double rand)
{
    // return a sample of the factor according to its potentials
    double sum = 0.0;
    int entry = 0;

    rand = rand * fact.p().sum();               // marginalize

    while (sum < rand)                          // sample
    {
        sum = sum + fact.p()[entry];
        entry++;
    }
    entry--;

    DEBUG(std::cout << "sampling from " << fact << " gives entry " << entry << " with prob " << fact.p()[entry] << std::endl;)
		
    return entry;                            // sample'd entry from factor
}

std::vector<unsigned int> getIntermediateVars(dai::FactorGraph fg, std::vector<unsigned int> hypothesis_vars, 
    std::vector<unsigned int> evidence_vars)
{
	// populate intermediate vars by matching all variables in fg with hypothesis and evidence variables

    std::vector<unsigned int> intermediateVars;
	std::vector<unsigned int>::iterator it;

	for (size_t i = 0; i < fg.nrVars(); i++ )
	{
		unsigned int var = fg.var(i).label();
		it = find(hypothesis_vars.begin(), hypothesis_vars.end(), var);
		if (it == hypothesis_vars.end())
		{
			it = find(evidence_vars.begin(), evidence_vars.end(), var);
			if (it == evidence_vars.end())
			{
				intermediateVars.push_back(var);
			}
		}
	}
    return intermediateVars;
}
