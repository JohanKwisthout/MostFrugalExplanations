/************************************************************************/
/* MFE, MAP independence, and Annealed MAP experimentation code	        */
/* Written by:		Johan Kwisthout                                		*/
/* Version:			1.2                             					*/
/* Last changed:	01-07-2022                                         	*/
/*                                                                     	*/
/* Version History:                                                    	*/
/* 1.2 Max Independence (weak and strong)                     		*/
/* 1.1 This version also implements MAP independence (Kwisthout, 2021)  */
/* 1.0 This version contains the MFE simulation code as well as an      */
/*     implementation of Annealed MAP.                                  */
/*                                                                    	*/
/* Version Comments:                                                   	*/
/*     Due to the limits of the current version of the libDAI library,  */
/*     the MAP is computed by computing the posterior over the          */
/*     hypothesis nodes and then finding the MAP brute-force.		    */
/*     For fair comparison, the number of MAP nodes should be limited.  */
/*                                                                    	*/
/************************************************************************/

// STL includes
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <chrono>
#include <ctime>

// headers
#include "mfesim.h"
#include "cxxopts.hpp"

// global values (with default values)
std::string inputfile = "./alarm.fg";
std::string outputfile = "./results";
std::vector<unsigned int> independenceTestVars;
std::vector<unsigned int> hypothesisVars;
std::vector<unsigned int> evidenceVars;
std::vector<unsigned int> evidenceValues;
std::vector<unsigned int> relevantVars;
std::vector<unsigned int> irrelevantVars;
std::vector<unsigned int> intermediateVars;

bool weakMapIndep = false;
bool strongMapIndep = false;
bool quantifiedMapIndep = false;
bool maxMapIndep = false;
bool exampleComputation = false;
bool mapComputation = false;
bool mapList = false;
bool annealedComputation = false;
bool relevanceComputationStandalone = false;
bool relevanceComputation = false;
bool mfeComputation = false;
unsigned long int cutoffTime = 3600;
unsigned long int samples = 100;
unsigned long int samplesRel = 10;
double relThreshold = 0.1;

int versionMajor = 1;
int versionMinor = 2;

// function prototypes
int main(int argc, char *argv[]);

// helper function for the command-line-operated program
cxxopts::ParseResult parse(int argc, char* argv[])
{
    const char *shortdes = "MAP, MFE, and Annealed MAP experimental simulation";
    try
    {
        cxxopts::Options options(argv[0], shortdes);
        options.add_options()
            ("i,input", "factor graph to run simulations on", cxxopts::value<std::string>())
            ("o,output", "output file for simulation results", cxxopts::value<std::string>())
            ("H,hypothesis-variables", "hypothesis variables", cxxopts::value<std::vector<unsigned int>>())
            ("E,evidence-variables", "evidence variables", cxxopts::value<std::vector<unsigned int>>())
            ("e,evidence-values", "values of the evidence variables", cxxopts::value<std::vector<unsigned int>>())
            ("D,independence-test", "variables 'R' to run independence test on", cxxopts::value<std::vector<unsigned int>>())
            ("R,relevant-variables", "relevant variables", cxxopts::value<std::vector<unsigned int>>())
            ("I,irrelevant-variables", "irrelevant variables", cxxopts::value<std::vector<unsigned int>>())
            ("r,relevance-computation", "do not use explicit relevant variables in MFE but assess them")
            ("S,relevance-samples", "number of samples to assess intermediate variables for relevance (0 = exact)", 
				cxxopts::value<unsigned long int>())
            ("t,relevance-threshold", "relevance threshold for inclusion", cxxopts::value<double>())
            ("s,samples", "number of samples to take from irrelevant variables", cxxopts::value<unsigned long int>())
            ("T,time", "cutoff time in seconds (0 = will run until big freeze", cxxopts::value<unsigned long int>())
            ("O,relevance-test", "run relevance test independent of MFE heuristic")
            ("A,annealed", "run Annealed MAP using reported parameters")
            ("M,map", "run exact MAP computation")
            ("m,map-list", "output all explanations with their probability")
            ("F,mfe", "run MFE heuristic")
            ("d,strong", "run Strong MAP-independence test")
            ("W,weak", "run Weak MAP-independence test")
            ("Q,quantified-indep", "run quantified independence tests")
            ("q,max-indep", "find maximum independent sets")
            ("x,example", "run MFE example with alarm.fg network")
            ("h,help", "display this help and exit")
            ("v,version", "output version information and exit")
        ;

        if (argc == 1)
        {
          std::cout << shortdes << std::endl;
          exit(0);
        }
    
        auto result = options.parse(argc, argv);

        if (result.count("help"))
        {
          std::cout << options.help({"", "Group"}) << std::endl;
          exit(0);
        }

        if (result.count("version"))
        {
            std::cout << "MAP-Indep, MFE and Annealed MAP experimental simulation version " << versionMajor << "." << versionMinor << std::endl << std::endl;
            std::cout << "The simulations implement the Most Frugal Explanation" << std::endl;
            std::cout << "heuristic (Kwisthout, 2015) and the Annealed MAP algorithm" << std::endl;
            std::cout << "(Yuan, Lu, and Druzdzel, 2004), as well as the MAP independence test" << std::endl;
            std::cout << "(Kwisthout, 2021) using the LibDAI library (Mooij, 2010)" << std::endl << std::endl;
            std::cout << "The programme and its source code are governed by a BSD-style license" << std::endl;
            std::cout << "that can be found in the LICENSE file." << std::endl;
            exit(0);
        }

        if (result.count("example"))
        {
            exampleComputation = true;  
            DEBUG(std::cout << "Example using the alarm network" << std::endl)
        }

        if (result.count("map"))
        {
            mapComputation = true;  
            DEBUG(std::cout << "Exact computation using MAP" << std::endl)
        }

        if (result.count("map-list"))
        {
            mapList = true;  
            DEBUG(std::cout << "Outputting all explanations" << std::endl)
        }

        if (result.count("relevance-test"))
        {
            relevanceComputationStandalone = true;  
            DEBUG(std::cout << "Assessing relevance of intermediate variables" << std::endl)
        }

        if (result.count("mfe"))
        {
            mfeComputation = true;  
            DEBUG(std::cout << "Heuristic using Most Frugal Explanation" << std::endl)
        }

        if (result.count("strong"))
        {
            strongMapIndep = true;  
            DEBUG(std::cout << "Running Strong MAP independence tests" << std::endl)
        }

        if (result.count("weak"))
        {
            weakMapIndep = true;  
            DEBUG(std::cout << "Running Weak MAP independence tests" << std::endl)
        }

        if (result.count("quantified-indep"))
        {
            quantifiedMapIndep = true;  
            DEBUG(std::cout << "Running Quantified MAP independence tests" << std::endl)
        }

        if (result.count("max-indep"))
        {
            maxMapIndep = true;  
            DEBUG(std::cout << "Finding maximum MAP independent sets" << std::endl)
        }

        if (result.count("annealed"))
        {
            annealedComputation = true;  
            DEBUG(std::cout << "Approximating using Annealed MAP" << std::endl)
        }

        if (result.count("input"))
        {
            inputfile = result["input"].as<std::string>();
            DEBUG(std::cout << "Input file: " << inputfile << std::endl)
        }

        if (result.count("output"))
        {
            outputfile = result["output"].as<std::string>();
            DEBUG(std::cout << "Output file: " << outputfile << std::endl)
        }

        if (result.count("relevance-computation"))
        {
            relevanceComputation = true;  
            DEBUG(std::cout << "Computing relevance explicitly" << std::endl)
        }

        if (result.count("relevance-samples"))
        {
            samplesRel = result["relevance-samples"].as<unsigned long int>();  
            DEBUG(std::cout << "Computing relevance using " << samplesRel << " samples" << std::endl)
        }

        if (result.count("time"))
        {
            cutoffTime = result["time"].as<unsigned long int>();  
            DEBUG(std::cout << "Cutoff time " << time << " seconds" << std::endl)
        }

        if (result.count("relevance-threshold"))
        {
            relThreshold = result["relevance-threshold"].as<double>();  
            DEBUG(std::cout << "Deciding relevance using threshold " << relThreshold << " for inclusion" << std::endl)
        }

        if (result.count("samples"))
        {
            samples = result["samples"].as<unsigned long int>();  
            DEBUG(std::cout << "Computing MFE using " << samples << " samples" << std::endl)
        }

        if (result.count("hypothesis-variables"))
        {  
            hypothesisVars = result["hypothesis-variables"].as<std::vector<unsigned int>>();
            DEBUG(
                std::cout << "Hypothesis variables: ";
                for (auto i = hypothesisVars.begin(); i != hypothesisVars.end(); ++i) std::cout << *i << ' ';
                std::cout << std::endl;
                 )
        }

        if (result.count("evidence-variables"))
        {  
            evidenceVars = result["evidence-variables"].as<std::vector<unsigned int>>();
            DEBUG(
                std::cout << "Evidence variables: ";
                for (auto i = evidenceVars.begin(); i != evidenceVars.end(); ++i) std::cout << *i << ' ';
                std::cout << std::endl;
                 )
        }

        if (result.count("evidence-values"))
        {  
            evidenceValues = result["evidence-values"].as<std::vector<unsigned int>>();
            DEBUG(
                std::cout << "Evidence values: ";
                for (auto i = evidenceValues.begin(); i != evidenceValues.end(); ++i) std::cout << *i << ' ';
                std::cout << std::endl;
                 )
        }

        if (result.count("relevant-variables"))
        {  
            relevantVars = result["relevant-variables"].as<std::vector<unsigned int>>();
            DEBUG(
                std::cout << "Relevant variables: ";
                for (auto i = relevantVars.begin(); i != relevantVars.end(); ++i) std::cout << *i << ' ';
                std::cout << std::endl;
                 )
        }

        if (result.count("irrelevant-variables"))
        {  
            irrelevantVars = result["irrelevant-variables"].as<std::vector<unsigned int>>();
            DEBUG(
                std::cout << "Irrelevant variables: ";
                for (auto i = irrelevantVars.begin(); i != irrelevantVars.end(); ++i) std::cout << *i << ' ';
                std::cout << std::endl;
                 )
        }

        if (result.count("independence-test"))
        {  
            independenceTestVars = result["independence-test"].as<std::vector<unsigned int>>();
            DEBUG(
                std::cout << "MAP independence test variables: ";
                for (auto i = independenceTestVars.begin(); i != independenceTestVars.end(); ++i) std::cout << *i << ' ';
                std::cout << std::endl;
                 )
        }

        return result;
    } 
    catch (const cxxopts::OptionException& e)
    {
        std::cout << "error parsing options: " << e.what() << std::endl;
        exit(1);
    }
}

int main(int argc, char *argv[])
{
    auto result = parse(argc, argv);
    auto arguments = result.arguments();

   	// random number generator
   	std::random_device rd;
    std::mt19937 gen(rd());						// random numbers by Mersenne twister algorithm

    // run an example of the computaions
	if (exampleComputation)
	{
    	dai::FactorGraph fg;
    	fg.ReadFromFile("./alarm.fg");

	    std::vector<unsigned int> ex_evidenceVars =        { 0, 1, 2, 8, 9,11,14,15,17,18,20,21,25,27,35,36};
	    std::vector<unsigned int> ex_evidenceValues =      { 1, 1, 2, 1, 1, 1, 1, 2, 2, 1, 1, 1, 2, 1, 1, 2};
	    std::vector<unsigned int> ex_hypothesisVars =      { 3, 5,12,13,16,22,24,26};
	    std::vector<unsigned int> ex_intermediateVars =    { 4, 6, 7,10,19,23,28,29,30,31,32,33,34};
	    std::vector<unsigned int> ex_indepTestVars =       {10,19,23,28,29};
		
		std::vector<unsigned int> ex_relevantVars;
		std::vector<unsigned int> ex_irrelevantVars;

	    for (auto inter = ex_intermediateVars.begin(); inter != ex_intermediateVars.end(); ++inter)
		{
    	    std::cout << "Relevance of " << *inter << " using 1000 samples equals ";
			auto start = std::chrono::steady_clock::now();
			double rel = relevance(fg, *inter, ex_evidenceVars, ex_evidenceValues, ex_hypothesisVars, ex_intermediateVars, 1000UL, gen);
			if (rel > 0.01) ex_relevantVars.push_back(*inter);
			else ex_irrelevantVars.push_back(*inter);
			auto end = std::chrono::steady_clock::now();
    	    std::cout << rel << std::endl;
			std::cout << "Computation took " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << " ns" << std::endl;
		}
		std::cout << "Relevant variables (threshold 0.01): " << ex_relevantVars << std::endl;
		std::cout << "Irrelevant variables (threshold 0.01): " << ex_irrelevantVars << std::endl;
		auto start = std::chrono::steady_clock::now();
		std::vector<unsigned long int> MAP = get_map(fg, ex_evidenceVars, ex_evidenceValues, ex_hypothesisVars, false); 
		auto end = std::chrono::steady_clock::now();
		std::cout << "MAP: " << MAP << std::endl;
		std::cout << "Computation took " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << " ns" << std::endl;
		start = std::chrono::steady_clock::now();
		std::vector<unsigned long int> MFE = compute_MFE(fg, ex_evidenceVars, ex_evidenceValues, ex_hypothesisVars, ex_relevantVars,
			ex_irrelevantVars, false, 0, 0, 2000, 3600); 
		end = std::chrono::steady_clock::now();
		std::cout << "MFE heuristic gives: " << MFE << std::endl;
		std::cout << "Computation took " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << " ns" << std::endl;
		start = std::chrono::steady_clock::now();
		std::vector<unsigned long int> ANN = annealed_map(fg, ex_evidenceVars, ex_evidenceValues, ex_hypothesisVars, 3600); 
		end = std::chrono::steady_clock::now();
		std::cout << "MAP approximated by Annealed MAP algorithm gives: " << ANN << std::endl;
		std::cout << "Computation took " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << " ns" << std::endl;
		start = std::chrono::steady_clock::now();

        std::vector<unsigned int> map;
        for (const unsigned long int &e: MAP) { map.push_back((unsigned int) e); }

		std::vector<unsigned long int> SIndep = max_strong_map_indep(fg, ex_evidenceVars, ex_evidenceValues, ex_hypothesisVars, map, ex_indepTestVars, 3600); 
		end = std::chrono::steady_clock::now();
		std::cout << "Strong MAP independent variables: " << SIndep << std::endl;
		std::cout << "Computation took " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << " ns" << std::endl;
		start = std::chrono::steady_clock::now();
		std::vector<unsigned long int> WIndep = max_weak_map_indep(fg, ex_evidenceVars, ex_evidenceValues, ex_hypothesisVars, map, ex_indepTestVars, 3600); 
		end = std::chrono::steady_clock::now();
		std::cout << "Weak MAP independent variables: " << WIndep << std::endl;
		std::cout << "Computation took " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << " ns" << std::endl;
        return 0;
	}

	time_t now = time(0);
   	dai::FactorGraph fg;
   	fg.ReadFromFile(inputfile.c_str());

	std::ofstream ofs;
	ofs.open (outputfile.c_str(), std::ofstream::out | std::ofstream::app);

	ofs << std::endl << "command: ";
    for (int i = 0; i < argc; i++)
        ofs << argv[i] << " ";
    ofs << std::endl;

	ofs << inputfile << " simulation results " << ctime(&now) << std::endl;
	ofs << "hypothesis vars " << hypothesisVars << std::endl;
	ofs << "evidence vars " << evidenceVars << " values " << evidenceValues << std::endl;
    intermediateVars = getIntermediateVars(fg, hypothesisVars, evidenceVars);
	ofs << "intermediate vars " << intermediateVars << std::endl;
    if ((strongMapIndep) || (weakMapIndep))
    	ofs << "independence test vars " << independenceTestVars << std::endl;

    // compute strong independence
    if (strongMapIndep)
    {
    	ofs << std::endl << "[STRONG] Strong MAP independence of subset of intermediate vars" << std::endl;
        std::vector<unsigned long int> map = get_map(fg, hypothesisVars, evidenceVars, evidenceValues, false);
        std::vector<unsigned int> hypValues;
        for (const unsigned long int &e: map) { hypValues.push_back((unsigned int) e); }
        std::vector<unsigned long int> strong;
        double q = 0.0;

   	    ofs << "[STRONG] ";

   		auto start = std::chrono::steady_clock::now();
        if ((quantifiedMapIndep == true) && (maxMapIndep == false))
        {
            q = strong_map_indep_measure(fg, evidenceVars, evidenceValues, hypothesisVars, hypValues, independenceTestVars, cutoffTime, false);
            ofs << "quantified: " << q;
        }
        else if ((quantifiedMapIndep == false) && (maxMapIndep == true))
        {
            strong = max_strong_map_indep(fg, evidenceVars, evidenceValues, hypothesisVars, hypValues, independenceTestVars, cutoffTime);
            ofs << "maximum independent set " << strong;
        }        
        else if ((quantifiedMapIndep == false) && (maxMapIndep == false))
        {
            ofs << strong_map_indep(fg, evidenceVars, evidenceValues, hypothesisVars, hypValues, independenceTestVars, cutoffTime);
        }
        else
        {
      		ofs << "illegal combination of switches!" << std::endl;
        }   		
   		auto end = std::chrono::steady_clock::now();

        ofs << std::endl;
  		ofs << "[STRONG] Computation took " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << " ns" << std::endl;
    }

    // compute weak independence
    if (weakMapIndep)
    {
    	ofs << std::endl << "[WEAK] Weak MAP independence of subset of intermediate vars" << std::endl;
        std::vector<unsigned long int> map = get_map(fg, hypothesisVars, evidenceVars, evidenceValues, false);
        std::vector<unsigned int> hypValues;
        for (const unsigned long int &e: map) { hypValues.push_back((unsigned int) e); }
        std::vector<unsigned long int> weak;
        double q = 0.0;

   	    ofs << "[WEAK] ";

   		auto start = std::chrono::steady_clock::now();
        if ((quantifiedMapIndep == true) && (maxMapIndep == false))
        {
            q = weak_map_indep_measure(fg, evidenceVars, evidenceValues, hypothesisVars, hypValues, independenceTestVars, cutoffTime, false);
            ofs << "quantified: " << q;
        }
        else if ((quantifiedMapIndep == false) && (maxMapIndep == true))
        {
            weak = max_weak_map_indep(fg, evidenceVars, evidenceValues, hypothesisVars, hypValues, independenceTestVars, cutoffTime);
            ofs << "maximum independent set " << weak;
        }        
        else if ((quantifiedMapIndep == false) && (maxMapIndep == false))
        {
            ofs << weak_map_indep(fg, evidenceVars, evidenceValues, hypothesisVars, hypValues, independenceTestVars, cutoffTime);
        }
        else
        {
      		ofs << "illegal combination of switches!" << std::endl;
        }   		
   		auto end = std::chrono::steady_clock::now();

        ofs << std::endl;
  		ofs << "[WEAK] Computation took " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << " ns" << std::endl;
    }

    // compute relevance (independent of MFE heuristic)
    if (relevanceComputationStandalone)
    {
    	ofs << std::endl << "[REL] Relevance assessment of intermediate vars" << std::endl;

        for (auto inter = intermediateVars.begin(); inter != intermediateVars.end(); ++inter)
		{
    	    ofs << "[REL] Relevance of " << *inter << " using " << samplesRel << " samples equals ";
			auto start = std::chrono::steady_clock::now();
			double rel = relevance(fg, *inter, evidenceVars, evidenceValues, hypothesisVars, intermediateVars, 
                samplesRel, gen);
			auto end = std::chrono::steady_clock::now();
    	    ofs << rel << std::endl;
			ofs << "[REL] Computation took " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << " ns" << std::endl;
		}
    }

    // compute exact MAP
    if (mapComputation)
    {
 		ofs << std::endl << "[MAP] MAP explanation of the hypotheses given the evidence is: ";

   		auto start = std::chrono::steady_clock::now();
   		std::vector<unsigned long int> map = get_map(fg, hypothesisVars, evidenceVars, evidenceValues, mapList);
   		auto end = std::chrono::steady_clock::now();

   	    ofs << map << std::endl;
  		ofs << "[MAP] Computation took " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << " ns" << std::endl;
    }

    // compute annealed MAP (using parameters reported in Yuan et al., 2004)
    if (annealedComputation)
    {
  		ofs << std::endl << "[ANN] Annealed MAP approximation gives: ";

   		auto start = std::chrono::steady_clock::now();
   		std::vector<unsigned long int> a_map = annealed_map(fg, hypothesisVars, evidenceVars, evidenceValues, cutoffTime);
   		auto end = std::chrono::steady_clock::now();

   	    ofs << a_map << std::endl;
   		ofs << "[ANN] Computation took " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << " ns" << std::endl;
    }

    // compute Most Frugal Explanation (Kwisthout, 2015)
    if (mfeComputation)
    {
		if (relevanceComputation)
		{
			ofs << std::endl << "[MFE] relevance of intermediate vars assessed using " << samplesRel << " samples" << std::endl;
			
			// all intermediate variables are irrelevant for now
            irrelevantVars = intermediateVars;
		}
		else
		{
			ofs << std::endl << "[MFE] relevant vars " << relevantVars << std::endl;
			ofs << "[MFE] irrelevant vars " << irrelevantVars << std::endl;
		}

	    ofs << std::endl << "[MFE] MFE of the hypotheses given the evidence based on " << samples << " samples is: ";

   		auto start = std::chrono::steady_clock::now();
    	std::vector<unsigned long int> mfe = compute_MFE(fg, evidenceVars, evidenceValues, hypothesisVars, relevantVars, irrelevantVars,
    		relevanceComputation, samplesRel, relThreshold, samples, cutoffTime);
   		auto end = std::chrono::steady_clock::now();

        ofs << mfe << std::endl;
    	ofs << "[MFE] Computation took " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << " ns" << std::endl;
	}
 
    ofs << std::endl;
	ofs.close();
    return 0;
}
