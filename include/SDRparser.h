#ifndef SDR_PARSER_H
#define SDR_PARSER_H

#include <fstream>
#include <string>
#include <vector>

#include "SDRtypes.h"

namespace karyoke
{

class SDRparser
{
public:

    SDRparser() = default;

    bool parseFile(const std::string& filename);				// General parent function to parse entire SDR file.

    const SDRmasterHeader& getMasterHeader() const;				// Getter function for parseMasterHeader function
    const std::vector<SDRsubHeader>& getSubHeaders() const;			// Getter function for parseSubHeader function

    bool writeSummary(const std::string& outputFilename) const;			// General parent function to write out the summary of the SDR file.


private:

    SDRmasterHeader masterHeader;
    std::vector<SDRsubHeader> subHeaders;

    // One-line pushback buffer, so if the function reads too far (next cell Subheader tag) can track the read position for other functions.
    std:: string pendingLine;
    bool hasPendingLine = false;
    bool nextLine(std::ifstream& file, std::string& line);
    void pushBackLine(const std::string& line);

    // Tracks which SDR file line is being processed, so error messages can point to it.
    int currentLineNumber = 0;          // 1-indexed line most recently returned by nextLine()
    std::string currentLineText;        // Raw text of that line
    void printErrorLocation() const;    // Prints "  --> SDR file line N: <text>"
    void printEndOfFileLocation() const;// For errors that happen when the file ends unexpectedly

    bool parseMasterHeader(std::ifstream& file);				// Master header parsing.

    bool parseSubHeader(std::ifstream& file, SDRsubHeader& subHeader);		// Cell subheader parsing.

    bool parseCellData(std::ifstream& file, SDRsubHeader& subHeader);		// Cell data parsing.

    bool parseMutatedChromosomeSizes(const std::string& value, std::map<int, double>& chromosomeSizes);		// Function to parse cell subheader mutated chromosome sizes entry.

    bool parseDataRecord(const std::string& line, SDRdataRecord& record);	// Parse one SDR data record/entry.

    bool parseFragment(const std::string& text, SDRfragment& fragment);		// Function to help parse SDR data record field 3 fragment information.

    // Write summary helper functions.
    void writeMasterHeaderSummary(std::ofstream& output) const;
    void writeSubHeaderSummary(std::ofstream& output, const SDRsubHeader& subHeader) const;
    void writeCellDataSummary(std::ofstream& output, const SDRsubHeader& subHeader) const;

};

} //namespace karyoke
#endif
