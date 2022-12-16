/************************************************************************/
/* MFE algorithm implementation              					        */
/* Written by:		Johan Kwisthout                                		*/
/* Version:			1.0													*/
/* Last changed:	22-05-2022                                         	*/
/*                                                                     	*/
/* Version History:                                                    	*/
/*                                                                     	*/
/* Version Comments:                                                   	*/
/* - implementation of the algorithm described in Yuan et al. (2004)   	*/
/************************************************************************/

// headers
#include "mfesim.h"
#include <cmath>
#include <chrono>

// constant parameters used in the paper
const double Tinit = 0.99;   	// initial temperature
const double alpha = 0.8;     	// cooling rate
const double K = 0.1;         	// used in formula (7) in the paper
const int iStopSteps = 20;	  	// stop if 20 rounds no change
const int iReheatSteps = 10;	// reheat if 10 rounds no improvement
const double kProb = 1.0;

const int iterations = 1000;    // if fixed i is used

//#define ITER
#define RFC

// implementation of algorithm AnnealedMAP in Yuan et al. (2004)
// code partially based on SimulatedMAP.cxx provided by Changhe Yuan,
// adjusted to work with libDAI rather than Genie (+ some simplifications)

std::vector<unsigned long int> annealed_map(dai::FactorGraph fg, std::vector<unsigned int> hypothesis_vars, std::vector<unsigned int> evidence_vars,
	std::vector<unsigned int> evidence_values, unsigned long int cutoffTime)
{
    unsigned long int timeBound = cutoffTime * 1000000000UL;

    std::vector<unsigned long int> map;			// current best map
    std::vector<double> map_scores;	    		// current best map
	
	double score = 1.0, cscore = 0.0;			// score of current best / working map

    std::random_device rd;                      // note C++11 code
    std::mt19937 gen(rd());						// random numbers by Mersenne twister algorithm
    std::uniform_real_distribution<> dis(0, 1); // uniform distribution between 0 and 1
    double u;                                   // current random number
    double T;                                   // current temperature
	double d;									// comparison value
	double specHeat = 0.0;						// from RFC technique
    double specT = Tinit;
	bool stopping = false;						// stopping criterion
    bool noChange = true;                       // used in RFC decision
    bool noIncrease = true;
	int	noChangeIterations = 0;
    int noIncreaseIterations = 0;
	int noIncreaseStop = 0;

	int i;                                      // iterations
	std::vector<double> currentScores;			// keep track of scores

    // internal time keeping to cut off computation after time bound
    auto start = std::chrono::steady_clock::now();

	// numbers relate to steps in the algorithm

    // 1. initialize X0, T0, and set i = 0
    T = Tinit;
    i = 0;
	map = local_prior_map(fg, hypothesis_vars, map_scores);
    for (auto const& s: map_scores)
		score *= s;

    DEBUG(std::cout << "prior MAP of Hyp " << map << " score " << score << std::endl;)
		
	cscore = score;
	currentScores.push_back(cscore);				// prior value

    // clamp evidence
    for (int j = 0; j < evidence_vars.size(); j++)
    {
        fg.clamp(evidence_vars[j], evidence_values[j], false);
    }
	
    // 2. while stopping rule is not satisfied
	while (!stopping)
	{	
        int xj_index = 0;
        dai::FactorGraph *fgc = fg.clone();      // make a local copy of fg to work with

		// 3. for each variable xj in hypothesis_vars do
		for (auto const& xj: hypothesis_vars)
		{   
			// 4. sample u ~ U[0,1]
			u = (double) dis(gen);

			// 5. sample xj proportional to its parents and evidence (using inference)
            dai::PropertySet opts;
            dai::JTree jt = dai::JTree(*fgc, opts("updates",std::string("HUGIN"))("inference",std::string("SUMPROD")));
            jt.init();
            jt.run();
            dai::Factor xjFact = jt.belief(fgc->var(xj));
            int xj_val = sample(xjFact, ((double)dis(gen)));

			// 6. accept sample according to temperature and probability
			d = (xjFact.p()[xj_val] / xjFact.p()[map[xj_index]]);

            noChange = true;
			noIncrease = true;

			if (d < 1)
			{
                DEBUG(std::cout << "u " << u << " d " << d << " d^(1/T-1) " << pow(d, 1/T -1) << std::endl;)

				if (u < pow(d, 1/T -1))
				{	

                    DEBUG(std::cout << "stochastic updating " << map[xj_index] << " to " << xj_val << std::endl;)
                    // adjust score and map value
        			cscore = cscore / map_scores[xj_index] * xjFact.p()[xj_val];
					map[xj_index] = xj_val;
                    map_scores[xj_index] = xjFact.p()[xj_val];
					noChange = false;
				}
			}
			else if (d > 1)
			{
                DEBUG(std::cout << "improvement: updating " << map[xj_index] << " to " << xj_val << std::endl;)
                 // adjust score and map value
     			cscore = cscore / map_scores[xj_index] * xjFact.p()[xj_val];
				map[xj_index] = xj_val;
                map_scores[xj_index] = xjFact.p()[xj_val];
				noIncrease = false;
				noChange = false;
			}
            else
            {
                DEBUG(std::cout << "nothing to update: d = 1 " << std::endl;)
            }
            fgc->clamp(xj, map[xj_index], false);
            
			// 7. keep track of best configuration so far (implicit in map[])

            xj_index++;
		}
		// end for
		if (noChange)
			noChangeIterations++;
		else
			noChangeIterations = 0;
		
		if (noIncrease)
		{
			noIncreaseIterations++;
			noIncreaseStop++;
		}
		else
		{
			noIncreaseIterations = 0;
			noIncreaseStop=0;
		}

        DEBUG(std::cout << i << ": current best " << map << " current max " << score << " actual score " << cscore << std::endl;)

		currentScores.push_back(cscore);
		if (cscore > score) score = cscore;

		// calculate heat
		double tmpSpecHeat = CalculateSpecHeat(currentScores, T, score);

        DEBUG(std::cout << "tmpSpecHeat " << tmpSpecHeat << " specHeat " << specHeat << " specT " << specT << std::endl;)

		if (tmpSpecHeat > specHeat) 
		{
			specHeat = tmpSpecHeat;
			specT = T;
		}

		// 8. adjust T using chosen annealing scheme
		T = T * alpha;			// linear decrease

		// reheating - RFC technique
		if (noIncreaseIterations >= iReheatSteps)
		{
			T = K*(1 - score) + specT;
			if (T >= 1.0) 
				T = Tinit;
            DEBUG(std::cout << "reheating " << std::endl;)

			noIncreaseIterations = 0;
		}
        #ifdef RFC
        DEBUG(std::cout << "noIncreaseStop " << noIncreaseStop << " iStopSteps " << iStopSteps << std::endl;)
		stopping = (noIncreaseStop > iStopSteps);
		#endif

        #ifdef ITER 
		stopping = (i > iterations);
		#endif

        if (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count() > timeBound)
        {
            std::cout << "stopping computation - time bound" << std::endl;
            stopping = true;
        }

        DEBUG(std::cout << "temperature is now " << T << std::endl;)
		
		// 9. increase i
		i++;
	}	
	// 10 end while

    return map;
}

double CalculateSpecHeat(const std::vector<double> &scores, const double &temperature, const double &bestScore)
{
	std::vector<double>::const_iterator itr;
	std::vector<double> probabilities;
	double aveCost = 0.0;
	double aveSqrCost = 0.0;
	double totalProb = 0.0;

	for (itr = scores.begin(); itr != scores.end(); itr++)
	{
		double tmpProb = exp(-1*(bestScore-(*itr))/(kProb*temperature));
		probabilities.push_back(tmpProb);
		totalProb += tmpProb;
	}

	for (int i = 0; i < probabilities.size(); i++)
	{
		probabilities[i] = probabilities[i]/totalProb;
		aveCost += (bestScore-scores[i])*probabilities[i];
		aveSqrCost += (bestScore-scores[i])*(bestScore-scores[i])*probabilities[i]; 
	}

	double var = aveSqrCost - aveCost*aveCost;
	if (var < 0) 
		var = 0;
	return var/(temperature*temperature);
}
