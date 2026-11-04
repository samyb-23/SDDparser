#ifndef SDD_PARSER_H
#define SDD_PARSER_H

#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <ostream>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "SDDtypes.h"

namespace karyoke
{

class SDDparser
{

public:

    bool load(const std::string& filename); 				// Function to check if initial SDD input file is loaded correctly.

    const Header& getHeader() const; 					// Function to get header objects

    const std::vector<Exposure>& getExposures() const; 			// Function to get Exposure objects 

    void printSummary(std::ostream& out) const; 			// Function that will run printHeaderSummary and printDataEntrySummary

    std::map<int, ChromosomeDamageSummary> summarizeChromosomeDamage(const Exposure& exposure) const; 		// Function that maps numBaseDamages, numSingleStrandBreaks, and numDoubleStranBreaks to each chromosome in a given exposure. 
													// Public so the web viewer's chromosome popups show the same counts as the summary file.


private:

    Header header;							// Instatiating the object header of struct Header

    bool parseHeader(std::ifstream& file); 				// Function to check if header entries parsed correctly

    bool parseDamageEntries(std::ifstream& file); 			// Function to check if exposure data entries parsed correctly

    std::vector<Exposure> exposures;					// Stores the chromosome ID and its corresponding damages per exposure.

    void printHeaderSummary(std::ostream& out) const; 			// Function to print the SDD header summary to an output file 

    void printExposureSummary(std::ostream& out) const; 		// Function to print exposure data summary to an output file.


};

} // namespace karyoke

#endif
