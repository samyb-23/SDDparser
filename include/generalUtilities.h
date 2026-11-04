#ifndef GENERAL_UTILITIES_H
#define GENERAL_UTILITIES_H

#include <cstddef>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <cctype>





namespace karyoke
{

// -------------------------------------------------- //
// ----------- GENERIC HELPER FUNCTIONS ------------- //
// -------------------------------------------------- //

// Trimming all white spaces in a line of an SDD/SDR file
inline std::string trim(const std::string& s)
{
    size_t first = s.find_first_not_of(" \t\r\n");

    if(first == std::string::npos)
    {
        return "";
    }

    size_t last = s.find_last_not_of(" \t\r\n");

    return s.substr(first,last-first+1);
}

// Splitting all SDD/SDR file lines by a delimiter of choice ',', ';', '/'.
inline std::vector<std::string> split(const std::string& line,char delimiter)
{
    std::vector<std::string> tokens;

    std::stringstream ss(line);

    std::string item;

    while(std::getline(ss,item,delimiter))
    {
	tokens.push_back(trim(item));
    }

    return tokens;
}

// Useful for header field comparison converting field names to lower case
inline std::string toLower(const std::string& input)
{
    std::string result = input;

    std::transform(result.begin(), 
	    result.end(), 
	    result.begin(), 
	    [](unsigned char c)
        {
            return static_cast<char>(std::tolower(c));
        }
    );
    return result;
}

// Removes all whitespaces in header field names and changes names to lower case for comparison.
inline std::string normalizeHeaderKey(const std::string& input)
{
    return toLower(trim(input));
}

// Helper function to parse vectors consisting of doubles
inline std::vector<double> parseDoubleList(const std::vector<std::string>& values)
{
    std::vector<double> output;

    for(auto& v : values)
    {
        if (v.empty())
	{    
	     continue;
	}
	output.push_back(std::stod(v));
    }

    return output;
}

// Helper function to parse vectors consiting of ints
inline std::vector<int> parseIntList(const std::vector<std::string>& values)
{
    std::vector<int> output;

    for(auto& v : values)
    {
	if (v.empty())
	{
	    continue;
	}
        output.push_back(std::stoi(v));
    }

    return output;
}

} // namespace karyoke

#endif
