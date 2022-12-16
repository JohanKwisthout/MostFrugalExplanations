/************************************************************************/
/* MFE relevance computation implementation    					        */
/* Written by:		Johan Kwisthout                                		*/
/* Version:			1.0													*/
/* Last changed:	18-01-2020                                         	*/
/*                                                                     	*/
/* Version History:                                                    	*/
/*                                                                     	*/
/* Version Comments:                                                   	*/
/* - Computes the relevance of an intermediate variable relative to  	*/
/*   evidence and hypothesis nodes.									  	*/
/************************************************************************/

// headers
#include "mfesim.h"
#include "dai/alldai.h"
#include "dai/jtree.h"

// compute the relevance of a variable
double relevance(dai::FactorGraph fg, unsigned int node, std::vector<unsigned int> evidence_vars, std::vector<unsigned int> evidence_values, 
	std::vector<unsigned int> hypothesis_vars, std::vector<unsigned int> intermediate_vars, unsigned long int samples, std::mt19937 rngen)
{
	// if samples = 0, relevance is computed exactly, otherwise by that amount of samples over the intermediate variables
	// algorithm: compute (approximate) the fraction of joint value assignments to the intermediate variables (other than node)
	// for which the value of node changes the MPE.

    unsigned long int non_equals = 0, iteration = 0, max_iterations = 1;
    unsigned int nr_int_vars = intermediate_vars.size();
    bool mpe_equal;

    // iterator (or sampler) over the intermediate values
    std::vector<unsigned int> intermediate_values;

    // maximum values of these variables
    std::vector<unsigned int> intermediate_max_values;

    // vector containing the MPEs (must be long because of libDAI)
    std::vector<unsigned long int> mpe, mpe_cmp;

    // vector of all evidence+intermediate variables (initialize this here)
    std::vector<unsigned int> ev_vars;
    for (auto inter: intermediate_vars)
    {
        ev_vars.push_back(inter);
    }
    std::copy(evidence_vars.begin(), evidence_vars.end(), back_inserter(ev_vars));

    // find out the node index of the variable for which we want to compute the evidence
    unsigned int node_index = std::distance(intermediate_vars.begin(), find(intermediate_vars.begin(), intermediate_vars.end(), node));
    if (node_index == nr_int_vars)
    {
        std::cerr << "Node " << node << " not found in intermediate_vars!" << std::endl;
        return 0.0;
    }

    // initialize values to zeros, set maximum values per variable, and set max_interations to total joint value assignments;
    for (auto inter: intermediate_vars)
    {
        unsigned int st = fg.var(inter).states();
         
        intermediate_values.push_back(0);
        intermediate_max_values.push_back(st - 1);

	    if (samples == 0)   // exact computation: iterate over all instantiations of intermediate vars, compute total iterations
		{
	        if (inter != node)
	        {
	            max_iterations *= st;
	        }
		}
    }

    if (samples != 0)
	{
		// set maximum iterations to the number of samples we want
		max_iterations = samples;
    }

    for (iteration = 1; iteration <= max_iterations; iteration++)
    {
        // set the first MPE
        // set value of node to test and copy the values of this sample to the total evidence
        std::vector<unsigned int> ev_values;
        std::copy(intermediate_values.begin(), intermediate_values.end(), back_inserter(ev_values));
        std::copy(evidence_values.begin(), evidence_values.end(), back_inserter(ev_values));

        ev_values[node_index] = 0;
        mpe_cmp = get_mpe(fg, ev_vars, ev_values);
                        
        // test whether MPE are all equal or not
        for (unsigned int i = 1; i <= intermediate_max_values[node_index]; i++)
        {
            mpe_equal = true;

            // set value of node to test and copy the values of this sample to the total evidence
            ev_values[node_index] = i;
            mpe = get_mpe(fg, ev_vars, ev_values);

            for (auto it: hypothesis_vars)
            {
				if (mpe[it] != mpe_cmp[it])
                {
       	            mpe_equal = false;
                    break;
                }
            }
        }

        if (mpe_equal == false)
        {
            non_equals++;
        }

		if (samples == 0)
		{
        	// next value in iteration
	        iterate(nr_int_vars - 1, node_index, intermediate_values, intermediate_max_values);
		}
		else
		{
        	// random sample
	        random_sample(nr_int_vars - 1, node_index, intermediate_values, intermediate_max_values, rngen);
		}
    }

	return ((double) non_equals) / ((double) (iteration - 1));
}


