#ifndef SDDPARSER_VERSION_H
#define SDDPARSER_VERSION_H

// SDDparser version, shown by "SDDparser -v" / "--version", at the start of every run, in the summary files and in the web viewer.
// Change it here only; everything else reads it from this header.
#define SDDPARSER_VERSION "1.0"

namespace sddparser
{

inline constexpr const char* SDDPARSER_VERSION_TEXT = "SDDparser v" SDDPARSER_VERSION;	// "SDDparser v1.0"

} // namespace sddparser

#endif
