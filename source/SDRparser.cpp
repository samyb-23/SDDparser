#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>

#include "SDRparser.h"
#include "generalUtilities.h"
#include "SDDparserVersion.h"
#include "SDRutilities.h"
#include "SDDutilities.h"			// cellCyclePhaseMeaning, shared with the SDD summary


namespace sddparser
{

bool SDRparser::parseFile(							// Main parser function that parses the SDR file
    const std::string& filename)						// Takes as input SDR file path
{
     std::ifstream file(filename);

    if (!file.is_open())							// If SDR file cannot open, exit with error
    {
        std::cerr << "Error: Could not open SDR file: "
                  << filename << std::endl;

        return false;
    }

    // Clear previous parsed data in case the same parser object is reused.
    masterHeader = SDRmasterHeader{};
    subHeaders.clear();
    hasPendingLine = false;			// Reset in case parser object gets reused across files.
    currentLineNumber = 0;
    currentLineText.clear();

    try						// Try/catch for extra layer of safety
    {
    	if (!parseMasterHeader(file))		// If master header cannot be parsed, exit with error
    	{
            return false;
    	}

    	std::string line;			// Define each line in SDR file to be looped through

	while (nextLine(file, line))
    	{
	    line = trim(line);			// Remove all whitespaces

	    if (line.empty())			// Skip empty lines
	    {
	    	continue;
	    }

	    if (line.rfind("***subheader - cell", 0) == 0)				// When subheader marker is found, run parseSubHeader function.
	    {
	    	SDRsubHeader subHeader{};						// Instantiate SDRsubHeader object subHeader.

	        if (!parseSubHeader(file, subHeader))					// If cell subheader cannot be parsed, exit with error, otherwise parse subheader.
	    	{
		    return false;
	    	}

	    	if (!parseCellData(file, subHeader))					// If cell data cannot be parsed exit with error, otherwise parseCellData
	    	{
		    return false;
	    	}

	    	subHeaders.push_back(subHeader);					// store subHeader object in subHeaders vector

	    }

    	}
    }
    catch (const std::exception& error)							// If problem arises while parsing, return general error
    {
	std::cerr << "ERROR: Unexpected error while parsing SDR file '" << filename << "': " << error.what() << "\n";
	printErrorLocation();
	return false;
    }

    return true;
}













// --------------------------------------------------------- //
// Getter functions for the master header and subheader to be
// summarized.
// --------------------------------------------------------- //
const SDRmasterHeader& SDRparser::getMasterHeader() const
{
    return masterHeader;
}

const std::vector<SDRsubHeader>& SDRparser::getSubHeaders() const
{
    return subHeaders;
}













// -------------------------------------------------------- //
// Parsing functions for individual SDR header, cell subheader and cell data fields
// -------------------------------------------------------- //
bool SDRparser::parseMasterHeader(
    std::ifstream& file)							// Takes as input the SDR file path
{

 std::string line;

    while (nextLine(file, line))						// Loop through every SDR file line
    {
        line = trim(line);

        if (line == "***end of master header***")				// If end of master header encountered, stop parsing master header.
        {
            return true;
        }

        if (line.empty())							// Skip all blank lines
        {
            continue;
        }

        const std::size_t delimiterPos = line.find(',');			// Split the header line into (field name, field value)

        if (delimiterPos == std::string::npos)					// If delimiter not found exit with error
        {
            std::cerr << "Error: Invalid SDR master header line (expected 'field, value;').\n";
	    printErrorLocation();
            return false;
        }

	const std::string field = normalizeHeaderKey(line.substr(0, delimiterPos)); // header key is before the ':' delimiter

	const std::vector<std::string> valueVec = split(line.substr(delimiterPos + 1), ';');
        const std::string value = valueVec.empty() ? "" : valueVec[0];		// header value is after the ':' delimiter

	// Associate SDR header keys with their corresponding values
        if (field == "sdr version")
        {
            masterHeader.sdrVersion = value;
        }
        else if (field == "author")
        {
            masterHeader.author = value;
        }
        else if (field == "associated sdd file")
        {
            masterHeader.associatedSDDFile = value;
        }
        else if (field == "cell cycle phase")					// Optional field (added in SDRv2.1); older SDR files without it still parse.
        {
            if (!value.empty())
            {
                try
                {
                    masterHeader.cellCyclePhase = parseDoubleList(split(value, ','));
                }
                catch (const std::exception& error)
                {
                    std::cerr << "ERROR: Invalid 'Cell Cycle Phase' value in SDR master header: " << value << " (" << error.what() << ")\n";
                    printErrorLocation();
                    return false;
                }
            }
        }
        else if (field == "intact chromosome sizes")
        {
            try
	    {
	    	masterHeader.intactChromosomeSizes = parseDoubleList(split(value, ','));	// parseDoubleList has no error catching, so introduce try/catch incase non-integer elements are encountered
	    }
	    catch (const std::exception& error)							// If invalid value encountered, exit with error.
	    {
		std::cerr << "ERROR: Invalid 'Intact Chromosome Sizes' value in SDR master header: " << value << " (" << error.what() << ")\n";
		printErrorLocation();
		return false;
	    }

        }
        else
        {
            std::cerr << "Warning: Unknown SDR master header field: " << field << "\n";
	    printErrorLocation();
        }
    }

    std::cerr << "Error: SDR master header terminator was not found." << "\n"; 		    // Reached EOF without finding the master-header terminator.
    printEndOfFileLocation();
    return false;
}













// Function to parse each SDR subheader associated with each cell in a sample
bool SDRparser::parseSubHeader(
    std::ifstream& file, 									// Takes as input the SDR input file path
    SDRsubHeader& subHeader)									// Parses a given cell's subheader section
{

    std::string line;

    while (nextLine(file, line))							// Loop through SDR file lines
    {

	line = trim(line);								// Remove all whitespaces

	if (line.empty())								// Skip blank lines
	{
	    continue;
	}

        if (line.rfind("***data - cell", 0) == 0)					// The subheader ends when the data section begins.
        {
            return true;
        }

	const std::size_t delimiterPos = line.find(',');				// Key , value pairs are separated by ',' delimiter

	if (delimiterPos == std::string::npos)						// If delimiter not found, exit with error
	{
	    std::cerr << "ERROR: Invalid SDR subheader line (expected 'field, value;').\n";
	    printErrorLocation();
	    return false;
	}


	const std::string field = normalizeHeaderKey(line.substr(0, delimiterPos));	// field is the string before the delimiter

	const std::vector<std::string> valueVec = split(line.substr(delimiterPos + 1), ';');
	const std::string value = valueVec.empty() ? "" : valueVec[0];			// value is the string after the delimiter

	if (field == "cell id")
	{
	    const std::vector<std::string> values = split(value, ',');			// Split value string by removing ',' delimiter and storing to a vector.

	    if (values.size() != 1)							// Each subheader should present data for only one cell ID.
	    {
		std::cerr << "ERROR: Invalid Cell ID in SDR subheader.\n";
		printErrorLocation();
		return false;
	    }

	    try
	    {
		subHeader.cellID = parseIntList(values)[0];				// try/catch for parseIntList which does not accept non-integer elements in a vector.
	    }
	    catch (const std::exception& error)
	    {
		std::cerr << "ERROR: Invalid 'Cell ID' value in SDR subheader: " << value << " (" << error.what() << ")\n";
		printErrorLocation();
	    	return false;
	    }

	}
	else if (field == "mutated chromosome sizes")
	{
	    if (value.empty())
	    {
		// Not measured, leave mutatedChromosomeSizes as its default-constructed empty map.
	    }

	    else if (!parseMutatedChromosomeSizes(value, subHeader.mutatedChromosomeSizes))		// Use helper function parseMutatedChromosomeSizes to see if the entries are formatted correctly.
	    {
		return false;
	    }

	}
	else if (field == "intact strands id")
	{
	    if (value.empty())
	    {
		// Not measured - leave intactStrandIDs as its default-constructed empty vector.
	    }
	    else
	    {
	    	try
	    	{
	    	    subHeader.intactStrandIDs = parseIntList(split(value, ','));		// try/catch for parseIntList in case non-integer is encountered in vector.
	    	}
	    	catch (const std::exception& error)
 	    	{
		    std::cerr << "ERROR: Invalid 'Intact Strands ID' value in SDR subheader: " << value << " (" << error.what() << ")\n";
		    printErrorLocation();
		    return false;
	    	}
	    }

	}
	else if (field == "total dsb count")
	{
	    if (value.empty())
	    {
		subHeader.totalDSBcount = SDR_FIELD_NOT_MEASURED;			// Blank/absent field means not measured, not an error.
	    }
	    else
	    {
	        try
		{
		    subHeader.totalDSBcount = parseIntList(split(value, ',')).at(0);	// try/catch for parseIntList in case non-integer is encountered in vector
		}
		catch (const std::exception& error)
		{
		    std::cerr << "ERROR: Invalid 'Total DSB count' value in SDR subheader: " << value << " (" << error.what() << ")\n";
		    printErrorLocation();
		    return false;
		}
	    }
	}
	else if (field == "total misrepair count")
	{
	    if (value.empty())
	    {
		subHeader.totalMisrepairCount = SDR_FIELD_NOT_MEASURED;			// Blank/absent field means not measured, not an error.
	    }
	    else
	    {
		try
		{
	    	    subHeader.totalMisrepairCount = parseIntList(split(value, ',')).at(0); // try/catch for parseIntList in case non-integer is encountered in vector.
	    	}
		catch (const std::exception& error)
		{
		    std::cerr << "ERROR: Invalid 'Total Misrepair Count' value in SDR subheader: " << value << " (" << error.what () << ")\n";
		    printErrorLocation();
		    return false;
		}
	    }
	}
	else if (field == "medras-mc log")
	{

/*	    if (value.empty())
    	    {
        	// Not measured, leave medrasMClog as its default-constructed empty vector.
    	    }
	    else
	    {
		try									// Blank value already handled through split, no special error
	    	{
	            subHeader.medrasMClog = parseIntList(split(value, ','));		// try/catch for parseIntList in case non-integer is encountered in vector
	    	}
	    	catch (const std::exception& error)
	    	{
		    std::cerr << "ERROR: Invalid 'MEDRAS-MC Log' value in SDR subheader: " << value << " (" << error.what() << ")\n";
		    printErrorLocation();
		    return false;
		}
	    }
*/
		subHeader.medrasMClog = value;

	}
	else										// If unknown subheader key encountered, exit with error.
	{
	    std::cerr << "Warning: Unknown SDR subheader field: " << field << "\n";
	    printErrorLocation();
	}

    }

    std::cerr << "ERROR: SDR data section was not found for cell " << subHeader.cellID << "\n"; // If ***data - cell not found, return error.
    printEndOfFileLocation();
    return false;
}











// Parse subHeader works in tandem with parseCellData, because each cell has both a subheader and a data block
bool SDRparser::parseCellData(
    std::ifstream& file, 								// Takes as input the SDR file path
    SDRsubHeader& subHeader)								// Parses the subheader associated with each cell in the sample, if subheader exists, data block should too
{
    std::string line;

    while (nextLine(file, line))							// Loop through SDR file lines
    {
	line = trim(line);								// Remove line whitespaces

	if(line.empty())								// Skip empty lines
	{
	    continue;
	}

	if (line.rfind("***subheader - cell", 0) == 0)					// If we read into the next cell's subheader, do not consume that line, give it back and stop reading.
	{
	    pushBackLine(line);
	    return true;
	}


	SDRdataRecord record{};								// Initialize SDRdataRecord object called 'record'

	try
	{
	    if (!parseDataRecord(line, record))						// If cannot parse SDR data, exit with error
	    {
	    	return false;
	    }
	}
	catch (const std::exception& error)
	{
	    std::cerr << "ERROR: Invalid SDR data record (" << error.what() << ")\n";
	    printErrorLocation();
	    return false;
	}

	subHeader.dataRecords.push_back(record);					// If parsed correctly, store the data record

    }

    return true;
}











// Function to parse the specific fields of the SDR data records/entries
bool SDRparser::parseDataRecord(
    const std::string& line, 								// Takes as input each line within a cell's data block
    SDRdataRecord& record)								// Data record containing new and old strand IDs, fragment ranges, centromere presence, and fragment shape
{

    // Fields separated by ';'
    // Field 1: cell ID
    // Field 2: New strand ID
    // Field 3: Fragment Strand - Old Strand ID / Old Strand Start / Old Strand End / has centromere
    // Field 4: new strand is linear (1) or is circular (0)

    // ---------------------------------------------
    // Split Data entries into the four fields delimited by ';'
    // ---------------------------------------------
    const std::vector<std::string> fields = split(line, ';');

    if (fields.size() != 4)								// SDR data fields must have exactly four fields, otherwise error.
    {
	std::cerr << "ERROR: Invalid SDR data record, 4 fields expected.\n";
	printErrorLocation();
	return false;
    }


    // ---------------------------------------------
    // Check if Field 1/2 in data block are present and are valid
    // ---------------------------------------------
    const std::vector<int> cellIDvals = parseIntList({trim(fields[0])});		// Field 1 in each record - remove whitespaces
    const std::vector<int> strandIDvals = parseIntList({trim(fields[1])});		// Field 2 in each record - remove whitespaces

    if (cellIDvals.size() != 1 || strandIDvals.size() != 1)				// Field 1 and 2 must each have exactly one element
    {
	std::cerr << "ERROR: Invalid cell ID / new strand ID in SDR data record. \n";
	printErrorLocation();
	return false;
    }

    record.cellID = cellIDvals[0];							// If cellIDvals are formatted correctly, append cell IDs to the record object.
    record.newStrandID = strandIDvals[0];						// If strandIDvals are formatted correctly, append strand IDs to the record object


    // -------------------------------------------------
    // Split subfields of field 3 fragment strand
    // -------------------------------------------------
    const std::vector<std::string> fragmentTextVec = split(fields[2], ',');		// Split Field 3 fragments subfields by delimiter ',' which separate multiple fragments.
    record.fragments.clear();								// Clear previous fragments records

    for (const std::string& fragmentText: fragmentTextVec)				// Loop through all fragments in a given Field 3 data record/entry.
    {
	if (trim(fragmentText).empty())							// Skip empty lines
	{
	    continue;
	}

	SDRfragment fragment{};								// Initialize SDRfragment object called 'fragment'

	if (!parseFragment(fragmentText, fragment))					// Check if field 3 is formatted correctly/can be parsed, if not, exit with error.
	{
	    return false;
	}

	record.fragments.push_back(fragment);						// If fragment records parsed correctly store them in fragments vector.

    }

    if (record.fragments.empty())
    {
	std::cerr << "ERROR: SDR data record has no fragments.\n";			// If no fragments found, exit with error.
	printErrorLocation();
	return false;
    }


    // ---------------------------------------------------
    // Check if correct Field 4 data entry
    // ---------------------------------------------------
    const std::vector<int> isLinearVals = parseIntList({trim(fields[3])});		// Field 4 - remove whitespaces

    if (isLinearVals.size() != 1 || (isLinearVals[0] != 0 && isLinearVals[0] != 1))	// Field 4 must be of size 1 and can only be a value of either 0 or 1
    {
	std::cerr << "ERROR: Invalid linear/circular flag in SDR data record Field 4.\n";
	printErrorLocation();
	return false;
    }

    record.linear = (isLinearVals[0] == 1);						// If field 4 formatted correctly, store if fragments are linear (1) or circular (0).

    return true;
}








// Helper function to parse SDR subheader 'Mutated Chromosome Sizes' key
bool SDRparser::parseMutatedChromosomeSizes(
    const std::string& value, 								// Reads each mutated chromosome value entry in the cell subheader
    std::map<int, double>& chromosomeSizes)						// Maps the chromosome/strand ID to the mutated chromosome size
{
    chromosomeSizes.clear();								// Clear any previously stored chromosome sizes

    std::string contents = trim(value);							// Remove white spaces

    if (contents.size() < 2 || contents.front() != '{' || contents.back() != '}')	// If not formatted correctly, error
    {
	std::cerr << "ERROR: Invalid Mutated Chromosome Sizes format: " << value << "\n";
	printErrorLocation();
	return false;
    }

    contents = contents.substr(1, contents.size() - 2);					// Remove first and last curly brace '{', '}'
    contents = trim(contents);								// Remove white spaces

    if (contents.empty())
    {
	return true;
    }

    const std::vector<std::string> entryVec = split(contents, ',');			// Split the dictionary grouping chromosome ID and chromosome size

    for (const std::string& entry: entryVec)						// Loop through every chromosome ID and chromosome size pair
    {
	const size_t delimiterPos = entry.find(':');

	if (delimiterPos == std::string::npos)						// If no delimiter ':' found or incorrect delimiter, exit with error.
	{
	    std::cerr << "ERROR: Invalid mutated chromosome size entry: " << entry << "\n";
	    printErrorLocation();
	    return false;
	}

	const std::string chromID = trim(entry.substr(0, delimiterPos));		// Chromosome IDs are before the ':' delimiter
	const std::string chromSize = trim(entry.substr(delimiterPos + 1));		// Chromosome sizes are after the ':' delimiter
	if (chromID.empty() || chromSize.empty())					// If either are empty, exit with error
	{
	    std::cerr << "ERROR: Invalid muated chromosome size entry: " << entry << "\n";
	    printErrorLocation();
	    return false;
	}

	std::vector<int> IDs;
	std::vector<double> sizes;

	try
	{
	    IDs = parseIntList({chromID});						// If parseIntList encounters non-integer element, exit with error
	    sizes = parseDoubleList({chromSize});					// If parseDoubleList encouters non-double element, exit with error
	}
	catch (const std::exception& error)
	{
	    std::cerr << "ERROR: Invalid mutated  chromosome size entry: " << entry << " (" << error.what() << ")\n";
	    printErrorLocation();
	    return false;
	}


	if (IDs.size() != 1 || sizes.size() != 1)					// Each pair should have 1 chrom ID and 1 chrom size, otherwise, error.
	{
	    std::cerr << "ERROR: Could not parse mutated chromosome size entry: " << entry << "\n";
	    printErrorLocation();
	    return false;
	}

	chromosomeSizes[IDs[0]] = sizes[0];						// Pair the corresponding chromosome ID to its chromosome size
    }

    return true;
}










// Function specifically to parse SDR data field 3 fragments
bool SDRparser::parseFragment(
    const std::string& text, 								// Takes as input the raw field 3 string
    SDRfragment& fragment)								// Takes as input the fragment object containing old strand ID, fragment start and end positions (Mbp), and if the fragment has centromere
{
    // Subfields separated by '/'
    // oldStrandID / oldStartPosition / oldEndPosition / hasCentromere
    const std::vector<std::string> subfields = split(text, '/');
    if (subfields.size() != 4)								// Require exactly four subfields delimited by '/' in field 3, otherwise error.
    {
        std::cerr << "ERROR: Invalid SDR fragment (expected 4 subfields): " << text << "\n";
	printErrorLocation();
        return false;
    }

    const std::vector<int> oldStrandIDvals = parseIntList({trim(subfields[0])});	// First subfield entry is old strand ID
    const std::vector<double> startVals = parseDoubleList({trim(subfields[1])});	// Second subfield entry is old strand start position in Mbp
    const std::vector<double> endVals = parseDoubleList({trim(subfields[2])});		// Third subfield entry is old strand end position in Mbp
    const std::vector<int> centromereVals = parseIntList({trim(subfields[3])});		// Fourth subfield entry is if the fragment has a centromere or not

    // Require exactly one of each subfield entry, no more, no less, otherwise error.
    if (oldStrandIDvals.size() != 1 || startVals.size() != 1 || endVals.size() != 1 || centromereVals.size() != 1)
    {
        std::cerr << "ERROR: Could not parse SDR fragment: " << text << "\n";
	printErrorLocation();
        return false;
    }

    // hasCentromere must be either a 0 or 1, otherwise exit with error.
    if (centromereVals[0] != 0 && centromereVals[0] != 1)
    {
        std::cerr << "ERROR: Invalid has-centromere flag in SDR fragment: " << text << "\n";
	printErrorLocation();
        return false;
    }

    fragment.oldStrandID = oldStrandIDvals[0];
    fragment.oldStartPosition = startVals[0];
    fragment.oldEndPosition = endVals[0];
    fragment.hasCentromere = (centromereVals[0] == 1);

    return true;
}














// ------------------------------------------------- //
// Write Summary functions for the master header
// the cell subheaders, and cell data fields, and
// the final summary file.
// ------------------------------------------------- //

bool SDRparser::writeSummary(
    const std::string& outputFilename							// Takes as input the output file summary path
    ) const
{
    std::ofstream output(outputFilename);

    if (!output.is_open())
    {
        std::cerr << "Error: Could not open SDR summary file: "
                  << outputFilename << std::endl;

        return false;
    }

    writeMasterHeaderSummary(output);

    if (!subHeaders.empty())							// Explains every cell's "Mutations by chromosome" list, so it is written once here rather than per cell.
    {
        output << "Note on \"Mutations by chromosome\": a mutation involving several chromosomes is listed under each of them,\n"
               << "except dicentrics, which are listed only under the chromosome drawn with two centromeres.\n"
               << "\"Complex Genomic Rearrangements\" counts clusters with 4 or more pieces of a chromosome out of order or with broken-up\n"
               << "material exchanged between chromosomes, plus one for every chromosome with a change that none of the other mutation\n"
               << "types accounts for (e.g. an inversion reaching a chromosome end), even when it also has other mutations.\n\n";
    }

    for (const SDRsubHeader& subHeader : subHeaders)
    {
        writeSubHeaderSummary(output, subHeader);
    	writeCellDataSummary(output, subHeader);
    }

    return true;
}











void SDRparser::writeMasterHeaderSummary(
    std::ofstream& output) const						// Takes as input the output file summary path
{
    output << "===================================\n";
    output << "SDR MASTER HEADER SUMMARY\n";
    output << "===================================\n";
    output << "Generated by " << SDDPARSER_VERSION_TEXT << "\n\n";

    output << "SDR Version: " << masterHeader.sdrVersion << "\n";
    output << "Author: " << masterHeader.author << "\n";
    output << "Associated SDD File: " << masterHeader.associatedSDDFile << "\n";
    output << "Cell Cycle Phase: " << cellCyclePhaseMeaning(masterHeader.cellCyclePhase) << "\n\n";	// Same interpretation as the SDD summary

    if (masterHeader.intactChromosomeSizes.empty())
    {
        output << "Intact Chromosome Sizes: none specified.\n\n";
        return;
    }

    // Index 0 is the declared chromosome count, not a size.
    const double declaredCount = masterHeader.intactChromosomeSizes[0];
    const std::size_t numSizes = masterHeader.intactChromosomeSizes.size() - 1;

    output << "Declared Chromosome Count: " << declaredCount << "\n";
    output << "Number Of Chromosome Sizes Listed: " << numSizes << "\n";
    output << "Intact Chromosome Sizes (Mbp): ";

    for (std::size_t i = 1; i < masterHeader.intactChromosomeSizes.size(); ++i)
    {
        output << masterHeader.intactChromosomeSizes[i];

        if (i + 1 < masterHeader.intactChromosomeSizes.size())
        {
            output << ", ";
        }
    }

    output << "\n\n";
}










void SDRparser::writeSubHeaderSummary(
    std::ofstream& output, 							// Takes as input the output file summary path
    const SDRsubHeader& subHeader						// Accesses SDR cell subheader to summarize subheader entries
    ) const
{
    output << "-----------------------------------\n";
    output << "CELL " << subHeader.cellID << " SUMMARY\n";
    output << "-----------------------------------\n\n";

    output << "Total DSB Count: ";
    if (subHeader.totalDSBcount == SDR_FIELD_NOT_MEASURED)
    {
        output << "N/A";
    }
    else
    {
        output << subHeader.totalDSBcount;
    }
    output << "\n";

    output << "Total Misrepair Count: ";
    if (subHeader.totalMisrepairCount == SDR_FIELD_NOT_MEASURED)
    {
        output << "N/A";
    }
    else
    {
        output << subHeader.totalMisrepairCount;
    }
    output << "\n\n";

}











void SDRparser::writeCellDataSummary(
    std::ofstream& output, 							// Takes as input the input file summary path
    const SDRsubHeader& subHeader						// Accepts the SDR subheader to detect mutations and summarize them
    ) const
{

    // New-strand records numbered below this threshold are baseline restatements of an original chromosome, not rearrangement/mutated outcomes.
    const int numOriginalStrands = masterHeader.intactChromosomeSizes.empty() ? 0 : static_cast<int>(masterHeader.intactChromosomeSizes[0]);

    // Declare mutation types to summarize the types of each mutation fed into SDDparser
    const std::vector<SDRdeletionEvent> deletions = detectDeletions(subHeader, numOriginalStrands, masterHeader);
    const std::vector<SDRinversionEvent> inversions = detectInversions(subHeader, numOriginalStrands, masterHeader);
    const std::vector<SDRtranslocationEvent> translocations = detectTranslocations(subHeader, numOriginalStrands, masterHeader);
    const std::vector<SDRecDNAevent> ecDNAs = detectECDNA(subHeader, numOriginalStrands, masterHeader);
    const std::vector<SDRdeletionInversionEvent> delInvs = detectDeletionInversions(subHeader, numOriginalStrands, masterHeader);
    const std::vector<SDRdeletionTranslocationEvent> delTras = detectDeletionTranslocations(subHeader, numOriginalStrands, masterHeader);
    const std::vector<SDRdeletionInsertionEvent> delInsertions = detectDeletionInsertions(subHeader, numOriginalStrands);
    const std::vector<SDRchromoplexyEvent> chromoplexy = detectChromoplexy(subHeader, numOriginalStrands, masterHeader);
    const std::vector<SDRchromothripsisEvent> chromothripsis = detectChromothripsis(subHeader, numOriginalStrands);
    const std::size_t complexRearrangements = countComplexGenomicRearrangements(subHeader, masterHeader);	// Includes changes no other type accounts for
    const std::vector<SDRchromosomeRingEvent> chromosomeRings = detectChromosomeRings(subHeader, numOriginalStrands);
    const std::vector<SDRdicentricEvent> dicentrics = detectDicentrics(subHeader, numOriginalStrands);

    output << "Mutation Summary:\n";
    output << "  Long Deletions: " << deletions.size() << "\n";
    output << "  Balanced Inversions: " << inversions.size() << "\n";
    output << "  Balanced Translocations: " << translocations.size() << "\n";
    output << "  Extrachromosomal DNA (ecDNA): " << ecDNAs.size() << "\n";
    output << "  Deletion-Inversions: " << delInvs.size() << "\n";
    output << "  Deletion-Translocations: " << delTras.size() << "\n";
    output << "  Deletion-Insertions: " << delInsertions.size() << "\n";
    output << "  Chromoplexy: " << chromoplexy.size() << "\n";
    output << "  Chromothripsis: " << chromothripsis.size() << "\n";
    output << "  Complex Genomic Rearrangements: " << complexRearrangements << "\n";
    output << "  Ring Chromosomes: " << chromosomeRings.size() << "\n";
    output << "  Dicentrics: " << dicentrics.size() << "\n";

    // Per-chromosome breakdown, in the same layout as the SDD summary's "Damage by chromosome". Only chromosomes with at
    // least one mutation are listed, and only the mutation types they have. Uses the same counting as the web viewer's
    // chromosome popups (see countMutationsPerOriginalStrand), so the file and the viewer always agree.
    const std::map<int, std::vector<std::pair<std::string, int>>> mutationsPerChromosome =
        countMutationsPerOriginalStrand(subHeader, masterHeader);

    if (mutationsPerChromosome.empty())
    {
        output << "\nNo chromosome-associated mutations found.\n";
    }
    else
    {
        output << "\nMutations by chromosome:\n\n";				// Counting rule explained once at the top of the file, see writeSummary.

        for (const auto& [chromosomeNumber, counts] : mutationsPerChromosome)
        {
            output << "Chromosome " << chromosomeNumber << "\n";

            for (const auto& [type, count] : counts)
            {
                output << "  " << type << ": " << count << "\n";
            }
        }
    }

    output << "\n";
}













// ----------------------------------------------------------- //
// Pushback helper functions for stream positioning
// ----------------------------------------------------------- //

bool SDRparser::nextLine(
    std::ifstream& file, 						// Takes as input input file path
    std::string& line)							// Takes as input raw line in the input file
{
    if (hasPendingLine)
    {
        line = pendingLine;
        hasPendingLine = false;
	++currentLineNumber;        					// Re-counting the line that pushBackLine() un-counted
        currentLineText = line;
        return true;
    }

    if (std::getline(file, line))
    {
        ++currentLineNumber;
        currentLineText = line;
        return true;
    }

    return false;
}



void SDRparser::pushBackLine(
    const std::string& line)						// Takes as input raw line in the input file
{
    pendingLine = line;
    hasPendingLine = true;
    --currentLineNumber;						// The caller hasn't consumed this line yet
}


void SDRparser::printErrorLocation() const
{
    std::cerr << "  --> SDR file line " << currentLineNumber << ": " << trim(currentLineText) << "\n";
}


void SDRparser::printEndOfFileLocation() const
{
    std::cerr << "  --> reached end of SDR file after line " << currentLineNumber << "\n";
}





} // namespace sddparser

