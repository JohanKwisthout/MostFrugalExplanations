/************************************************************************/
/* bif2fg .bif to .fg format translater      					        */
/* Written by:		Johan Kwisthout                                		*/
/* Version:			1.0													*/
/* Last changed:	17-01-2020                                         	*/
/*                                                                     	*/
/* Version History:                                                    	*/
/*                                                                     	*/
/* Version Comments:                                                   	*/
/* - use bif2fg network.bif network.fg                                 	*/
/* - quick and dirty, assuming correct .bif file, but does catch a few 	*/
/*   inconsistencies in the .bif file									*/
/* - the probs should change in the order P(a1|b2,c3,d4) as no check is */
/*   done on the actual values of the parents per line!                 */ 
/************************************************************************/

// comment for production mode, uncomment for debug messages
//#define DEBUGMODE

#ifdef DEBUGMODE
	#define DEBUG(a) a;
#else
	#define DEBUG(a) ;
#endif	

// STL includes
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <boost/tokenizer.hpp>

// data prototypes
std::vector<int> values;
std::vector<std::string> variables;
std::vector<std::vector<int>> factors;
std::vector<std::vector<double>> probabilities;

int variableCount = 0;
int index = 0;

// function prototypes
int main(int argc, char *argv[]);

int main(int argc, char *argv[])
{

    if (argc != 3)
    {
        std::cout << "bif2fg: translates .bif network into libDAI factor graph format." << std::endl;
        std::cout << "Use of this programme is governed by a BSD-style license" << std::endl;
        std::cout << "that can be found in the LICENSE file." << std::endl;
        std::cout << "Use: " << argv[0] << " network.bif network.fg" << std::endl;
        return 0;
    }

    try
    {
        std::ifstream inFile(argv[1]);
        std::ofstream outFile;
        outFile.open(argv[2], std::ios::binary | std::ios::trunc);
        outFile << "# file created with bif2fg utility; source file " << argv[1] << std::endl;

        while (inFile)
        {
			std::string currentLine;
            getline(inFile, currentLine);
			index = currentLine.find("variable ");

			// variable keyword found: add the name to the 'names' list
			if (index != std::string::npos)
			{
				index = index + sizeof("variable");
				std::string name = currentLine.substr(index, currentLine.find(' ', index) - index);
				DEBUG(std::cout << "Name found: " << name << std::endl;)
				variables.push_back(name);
				
				// find number of values for this variable
	            getline(inFile, currentLine);
				index = currentLine.find("[");
				int card = std::stoi(currentLine.substr(index + 1, currentLine.find("]", index) - index), nullptr, 10);
				DEBUG(std::cout << "Number of values found: " << card << std::endl;)
				values.push_back(card);
			}
			else
			{
				// probability keyword found: find the indices of the variables involved
				index = currentLine.find("probability");
				if (index != std::string::npos)
				{
					index = index + sizeof("probability");
					std::string vars = currentLine.substr(index, currentLine.find('{', index) - index);
					DEBUG(std::cout << "String to parse: " << vars << std::endl;)

					// scan the line for separators) which denote end of variable name
					boost::char_separator<char> sep("(,|);{} ");
					typedef boost::tokenizer< boost::char_separator<char> > t_tokenizer;
					t_tokenizer tok_var(vars, sep);

					std::vector<int> fact;

					for (t_tokenizer::iterator beg = tok_var.begin(); beg != tok_var.end(); ++beg)
					{
						auto it = std::find(variables.begin(), variables.end(), *beg);
						if (it == variables.end())
						{
							std::cerr << "Error: in CPT " << vars << ", variable name " << *beg << " not found" << std::endl;
							exit(1);
						}

						// find the index of the variable name and add it to the factor
						int var_no = std::distance(variables.begin(), it);
						fact.push_back(var_no);
						DEBUG(std::cout << *beg << " (" << var_no << ")" << std::endl;)
					}
					factors.push_back(fact);
					DEBUG(std::cout << "Factor: "; for (auto i: fact) std::cout << i << ' '; std::cout << std::endl;)

					// now read in the probabilities. Read from the input file until we get a ';' and just find all
					// doubles in the resulting string, ignoring other stuff like value names.

					std::string probs;
					do
					{
			            getline(inFile, currentLine);
						probs.append(currentLine);
					}
					while (currentLine.find('}') ==  std::string::npos);
					DEBUG(std::cout << "Probs string: " << probs << std::endl;)

					t_tokenizer tok_prob(probs, sep);
					std::vector<double> prob;

					for (t_tokenizer::iterator beg = tok_prob.begin(); beg != tok_prob.end(); ++beg)
					{
						try 
						{
							double p = std::stod(*beg);
							prob.push_back(p);
							DEBUG(std::cout << *beg << " (" << p << ")" << std::endl;)
						}
						catch (const std::invalid_argument& ia) 
						{
							DEBUG(std::cout << "Ignored: " << *beg << std::endl;)
						}
					}
					probabilities.push_back(prob);
					DEBUG(std::cout << "Potentials: "; for (auto i: prob) std::cout << i << ' '; std::cout << std::endl;)

				}
			}	
        }
		
		// we have now parsed the .bif input file, let's write the .fg output!
        outFile << factors.size() << std::endl;		// number of factors

		for (int i = 0; i < factors.size(); i++)
		{
			// start with empty line
	        outFile << std::endl;
			// number of variables in the factor
			outFile << factors[i].size() << std::endl;
			// variables involved in the factor	
			for (auto j: factors[i]) { outFile << j << ' '; } outFile << std::endl; 
			// cardinality of these variables
			for (auto j: factors[i]) { outFile << values[j] << ' '; } outFile << std::endl;
			// count the number of zero probabilities in the factor
			int zeros = std::count(probabilities[i].begin(), probabilities[i].end(), 0.0);
			outFile << probabilities[i].size() - zeros << std::endl;
			// put the non-zero probabilities in the factor
			int k = 0;
			for (std::vector<double>::iterator p = probabilities[i].begin(); p != probabilities[i].end(); ++p)
			{			
				if (*p != 0.0) outFile << k << ' ' << *p << std::endl;
				k++;
			}
		}

    }
    catch (std::ios_base::failure &e)
    {
        std::cerr << "File IO failure: " << e.what() << std::endl;
        exit(1);
    }
}
