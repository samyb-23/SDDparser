#include <iostream>
#include <fstream>
#include <filesystem>
#include <functional>
#include <vector>
#include <utility>
#include <map>
#include <string>
#include <cstdio>

#if defined(_WIN32)
#include <io.h>                                                 // _isatty, _fileno
#else
#include <unistd.h>                                             // isatty, fileno
#endif

#include "KaryoKEversion.h"
#include "SDDparser.h"
#include "SDRparser.h"
#include "Karyogram.h"
#include "KaryogramViewer.h"
#include "SDRutilities.h"

using namespace karyoke;

namespace
{
    void printUsage()
    {
        std::cout
            << KARYOKE_VERSION_TEXT << "\n\n"
            << "Usage:\n"
            << "  ./karyoke --sdd <SDD file> [-kh|-ko [--save-png] [--save-svg] [--save-zip]] [<output directory>]\n"
            << "  ./karyoke --sdr <SDR file> [-kh|-ko [--save-png] [--save-svg] [--save-zip]] [<output directory>]\n"
            << "  ./karyoke --help\n"
            << "  ./karyoke --version\n\n"
            << "The arguments must be given in this order:\n"
            << "  1. --sdd <file>   Input is an SDD file (DNA damage), or\n"
            << "     --sdr <file>   Input is an SDR file (mutations)\n"
            << "  2. -kh            Draw one karyogram per cell, human genome (realistic centromere positions), or\n"
            << "     -ko            Draw one karyogram per cell, other genome (generic centromere positions)\n"
            << "                    Optional. Without it, only the summary file is written. With it, karyograms open\n"
            << "                    in a web viewer in your browser unless a --save option follows.\n"
            << "  3. --save-png     Write PNG files to the output directory instead of opening the web viewer\n"
            << "     --save-svg     Write SVG files to the output directory instead of opening the web viewer\n"
            << "     --save-zip     Put all karyograms in one zip file in the output directory instead\n"
            << "                    (PNG by default; add --save-svg for SVG, or both flags for both)\n"
            << "                    Optional, only after -kh/-ko. The three can be combined, in any order among themselves.\n"
            << "  4. <directory>    Output folder for the summary file and karyograms, created if it does not exist.\n"
            << "                    Optional, always the last argument. Default: the folder of the input file.\n\n"
            << "Help:\n"
            << "  -h, --help        Show this help\n"
            << "  -v, --version     Show the KaryoKE version\n\n"
            << "Examples:\n"
            << "  ./karyoke --sdd exampleSDD.txt\n"
            << "  ./karyoke --sdd exampleSDD.txt -kh\n"
            << "  ./karyoke --sdr exampleSDR.txt\n"
            << "  ./karyoke --sdr exampleSDR.txt -ko\n"
            << "  ./karyoke --sdr exampleSDR.txt -kh --save-png\n"
            << "  ./karyoke --sdr exampleSDR.txt -kh --save-png --save-svg\n"
            << "  ./karyoke --sdr exampleSDR.txt -kh --save-zip\n"
            << "  ./karyoke --sdr exampleSDR.txt -kh --save-png results/\n"
            << "  ./karyoke --sdd exampleSDD.txt results/\n";
    }

    // Everything the user asked for on the command line.
    struct CommandLine
    {
        bool showHelp = false;					// -h or --help
        bool showVersion = false;				// -v or --version
        bool sdrInput = false;					// --sdr (true) or --sdd (false)
        std::string filename;					// Input file given after --sdd/--sdr
        std::string outputDirectory;				// Where the summary and karyograms go; empty = the input file's folder
        bool drawKaryogram = false;				// -kh or -ko
        std::string genomeType;					// "human" (-kh) or "other" (-ko)
        bool savePNG = false;					// --save-png
        bool saveSVG = false;					// --save-svg
        bool saveZip = false;					// --save-zip
    };

    bool isSaveOption(const std::string& arg)
    {
        return arg == "--save-png" || arg == "--save-svg" || arg == "--save-zip";
    }

    // Prints why an option that KaryoKE knows was found in the wrong place (or is from an earlier version), e.g.
    // "--save-png must come right after -kh or -ko". Anything else is reported as an unknown option.
    void reportMisplacedOption(const std::string& arg, const CommandLine& commandLine, bool afterOutputDirectory)
    {
        const std::string example = (commandLine.sdrInput ? "--sdr " : "--sdd ") + commandLine.filename;	// The user's own input, for the examples

        if (arg == "-h" || arg == "--help" || arg == "-v" || arg == "--version")
        {
            std::cerr << "ERROR: " << arg << " is used on its own: ./karyoke " << arg << "\n";
        }
        else if (arg == "--sdd" || arg == "--sdr")
        {
            std::cerr << "ERROR: Give one input file, with either --sdd or --sdr, as the first argument.\n";
        }
        else if (arg == "-kh" || arg == "-ko")
        {
            std::cerr << (commandLine.drawKaryogram && !afterOutputDirectory
                ? std::string("ERROR: Give -kh or -ko only once.\n")
                : "ERROR: " + arg + " must come right after the input file, e.g. " + example + " " + arg + "\n");
        }
        else if (isSaveOption(arg) && afterOutputDirectory)
        {
            std::cerr << "ERROR: " << arg << " must come before the output directory, which is always the last argument, e.g.\n"
                      << example << (commandLine.genomeType == "other" ? " -ko " : " -kh ") << arg << " results/\n";
        }
        else if (isSaveOption(arg))
        {
            std::cerr << "ERROR: " << arg << " was given more than once.\n";
        }
        else if (arg == "-sdr" || arg == "-sdd")
        {
            std::cerr << "ERROR: " << arg << " is now -" << arg << " (two dashes).\n";
        }
        else if (arg == "--karyogram")
        {
            std::cerr << "ERROR: --karyogram human|other is now -kh (human) or -ko (other).\n";
        }
        else if (arg == "-k" || arg == "-o" || arg == "-hk" || arg == "-ok")		// Near misses of -kh/-ko
        {
            std::cerr << "ERROR: Unknown command-line option: " << arg << "\n"
                      << "Karyograms are drawn with -kh (human genome) or -ko (other genome), written as one option.\n";
        }
        else
        {
            std::cerr << "ERROR: Unknown command-line option: " << arg << "\n";
        }
    }

    // Reads the command line into commandLine. The order is fixed:
    //   ./karyoke --sdd|--sdr <input file> [-kh|-ko [--save-png] [--save-svg] [--save-zip]] [<output directory>]
    //   ./karyoke -h|--help
    //   ./karyoke -v|--version
    // -kh/-ko write karyograms (human or other genome; always one option, so -h stays free for help). The --save options may only follow
    // -kh/-ko, in any order among themselves. The output directory, if given, is always the last argument.
    // The word after --sdd/--sdr is always taken as the file name, so "-" stays available for reading the data from the terminal later.
    // Prints an error and returns false on any malformed input. Expects argc >= 2.
    bool parseCommandLine(int argc, char* argv[], CommandLine& commandLine)
    {
        commandLine = CommandLine();
        const std::string first = argv[1];

        // Help and version: on their own
        if (first == "-h" || first == "--help" || first == "-v" || first == "--version")
        {
            if (argc > 2)
            {
                std::cerr << "ERROR: " << first << " is used on its own: ./karyoke " << first << "\n";
                return false;
            }

            commandLine.showHelp = (first == "-h" || first == "--help");
            commandLine.showVersion = !commandLine.showHelp;
            return true;
        }

        // 1. --sdd <file> or --sdr <file>
        if (first != "--sdd" && first != "--sdr")
        {
            std::error_code error;

            if (first[0] != '-' && std::filesystem::is_regular_file(first, error))	// An input file without --sdd/--sdr
            {
                std::cerr << "ERROR: Say what kind of file " << first << " is, as the first argument: --sdd " << first << " or --sdr " << first << "\n";
            }
            else if (first == "-sdr" || first == "-sdd")
            {
                std::cerr << "ERROR: " << first << " is now -" << first << " (two dashes).\n";
            }
            else
            {
                std::cerr << "ERROR: The first argument must be --sdd <SDD file> or --sdr <SDR file> (found " << first << ").\n";
            }

            return false;
        }

        if (argc < 3 || (argv[2][0] == '-' && std::string(argv[2]) != "-"))
        {
            std::cerr << "ERROR: " << first << " requires a file path, e.g. " << first << " example" << (first == "--sdd" ? "SDD" : "SDR") << ".txt\n";
            return false;
        }

        commandLine.sdrInput = (first == "--sdr");
        commandLine.filename = argv[2];
        int i = 3;

        // 2. Optional -kh or -ko
        if (i < argc && (std::string(argv[i]) == "-kh" || std::string(argv[i]) == "-ko"))
        {
            commandLine.drawKaryogram = true;
            commandLine.genomeType = (std::string(argv[i]) == "-kh") ? "human" : "other";
            i++;
        }

        // 3. Optional --save-png, --save-svg, --save-zip, only after -kh/-ko, in any order among themselves
        for (; i < argc && isSaveOption(argv[i]); i++)
        {
            const std::string arg = argv[i];

            if (!commandLine.drawKaryogram)
            {
                std::cerr << "ERROR: " << arg << " saves karyograms, so it needs -kh or -ko right before it, e.g. "
                          << first << " " << commandLine.filename << " -kh " << arg << "\n";
                return false;
            }

            bool& saveFlag = (arg == "--save-png") ? commandLine.savePNG : (arg == "--save-svg") ? commandLine.saveSVG : commandLine.saveZip;

            if (saveFlag)
            {
                std::cerr << "ERROR: " << arg << " was given more than once.\n";
                return false;
            }

            saveFlag = true;
        }

        // 4. Optional output directory, always the last argument
        if (i < argc)
        {
            const std::string arg = argv[i];

            if (arg[0] == '-')
            {
                reportMisplacedOption(arg, commandLine, false);
                return false;
            }

            std::error_code error;

            if (std::filesystem::is_regular_file(arg, error))		// Most likely a second input file, not a folder to write into
            {
                std::cerr << "ERROR: " << arg << " is a file, not an output directory.\n";
                return false;
            }

            commandLine.outputDirectory = arg;
            i++;
        }

        if (i < argc)							// Anything after the output directory
        {
            const std::string arg = argv[i];

            if (arg[0] == '-')
            {
                reportMisplacedOption(arg, commandLine, true);
            }
            else
            {
                std::cerr << "ERROR: Unexpected argument after the output directory " << commandLine.outputDirectory << ": " << arg << "\n"
                          << "Give only one output directory, as the last argument.\n";
            }

            return false;
        }

        return true;
    }

    // Turns the chromosome labels drawn on a karyogram into clickable areas for the web viewer, each with a summary
    // per homolog (SDR mutations or SDD damages). Only the types that are present are listed.
    std::vector<KaryogramHotspot> buildLabelHotspots(
	const std::vector<KaryogramLabelRegion>& labelRegions,	// From Karyogram::getLabelRegions(), after drawing the cell
	const std::map<int, std::vector<std::pair<std::string, int>>>& countsPerStrand)	// Strand/chromosome ID -> (type, count), nonzero only
    {
        std::vector<KaryogramHotspot> hotspots;

        for (const KaryogramLabelRegion& region : labelRegions)
        {
            KaryogramHotspot hotspot;
            const std::string name = region.label.rfind("Chr ", 0) == 0 ? region.label.substr(4) : region.label;	// "Chr 1" -> "1", "X" stays "X"
            hotspot.title = "Chromosome " + name;
            hotspot.x = region.x;
            hotspot.y = region.y;
            hotspot.width = region.width;
            hotspot.height = region.height;

            // Homologous layouts list both homologs of the pair; non-homologous layouts have one strand per label,
            // which needs no heading.
            for (size_t h = 0; h < region.oldStrandIDs.size(); h++)
            {
                KaryogramPopupSection section;

                if (region.oldStrandIDs.size() > 1)
                {
                    section.heading = "Homolog " + std::to_string(h + 1);
                }

                const auto found = countsPerStrand.find(region.oldStrandIDs[h]);

                if (found != countsPerStrand.end())
                {
                    section.counts = found->second;
                }

                hotspot.sections.push_back(std::move(section));
            }

            hotspots.push_back(std::move(hotspot));
        }

        return hotspots;
    }

    // Base damages, single-strand breaks and double-strand breaks per SDD chromosome number for one exposure, using the
    // same per-chromosome totals as the summary file. Sister chromatids of a homolog share one chromosome number, so
    // their damages are added together. Types with no damage are left out.
    std::map<int, std::vector<std::pair<std::string, int>>> countDamagesPerChromosome(
	const SDDparser& parser,
	const Exposure& exposure)
    {
        std::map<int, std::vector<std::pair<std::string, int>>> result;

        for (const auto& [chromosomeNumber, summary] : parser.summarizeChromosomeDamage(exposure))
        {
            std::vector<std::pair<std::string, int>>& counts = result[chromosomeNumber];

            if (summary.numBaseDamages > 0)        counts.emplace_back("Base damages", summary.numBaseDamages);
            if (summary.numSingleStrandBreaks > 0) counts.emplace_back("Single-strand breaks", summary.numSingleStrandBreaks);
            if (summary.numDoubleStrandBreaks > 0) counts.emplace_back("Double-strand breaks", summary.numDoubleStrandBreaks);
        }

        return result;
    }

    // Writes all karyograms of one input file into a single HTML viewer and opens it in the browser.
    bool showInWebViewer(
	const std::vector<KaryogramImage>& images,		// Karyograms generated in memory, one per cell
	const std::filesystem::path& inputPath,			// Input SDD/SDR file, used to name the viewer
	const std::filesystem::path& outputFolder,		// Output directory (default: the input file's folder)
	const std::string& genomeType)				// human or other, kept in the file name like the PNG outputs
    {
        const bool browserAvailable = canOpenBrowser();

        // Where to write the viewer: a "temp" folder in the output directory, created if needed. It keeps that folder
        // uncluttered while staying easy to find, reopen, or copy back from a remote machine (scp), unlike the system temp folder.
        std::filesystem::path viewerFolder = outputFolder / "temp";

        std::error_code error;
        std::filesystem::create_directories(viewerFolder, error);		// No error if the folder already exists

        if (error || !std::filesystem::is_directory(viewerFolder, error))	// Safeguard: e.g. a read-only folder, or a file named "temp"
        {
            std::cout << "Could not create the folder " << viewerFolder << ", writing the viewer to " << outputFolder << " instead.\n";
            viewerFolder = outputFolder;
        }

        // Absolute, so the printed path can be reopened from any folder and the browser (or Windows, from WSL) gets an unambiguous path.
        const std::filesystem::path viewerPath = std::filesystem::absolute(viewerFolder /
            (inputPath.stem().string() + "_karyograms_" + genomeType + ".html"), error);

        if (!writeKaryogramViewer(images, inputPath.filename().string(), viewerPath))
        {
            return false;
        }

        // Always print the path first, so the user can open the viewer by hand whatever happens next.
        std::cout << "Karyogram viewer written to: " << viewerPath << "\n";

        if (!browserAvailable)
        {
            std::cout << "No browser can be opened from this session (for example over SSH).\n"
                      << "Copy the file above to your own computer and open it there, or rerun with --save-png, --save-svg or --save-zip to get image files.\n";
        }
        else if (!openInBrowser(viewerPath))
        {
            std::cout << "Could not open a browser automatically, open the file above in your browser.\n";
        }

        return true;
    }


    // Shows "Processing cell 18/100 (Cell ID 17)..." while the cells are drawn. In a terminal the line is rewritten in place, so a file with
    // 100 cells does not fill the screen; when the output goes to a file or a pipe (e.g. a cluster job log), each cell gets its own line instead.
    class ProgressLine
    {
    public:
        ProgressLine()
#if defined(_WIN32)
            : inTerminal(_isatty(_fileno(stdout)) != 0)
#else
            : inTerminal(isatty(fileno(stdout)) != 0)
#endif
        {
        }

        void show(std::size_t position, std::size_t total, int cellID)
        {
            const std::string text = "Processing cell " + std::to_string(position) + "/" + std::to_string(total) + " (Cell ID " + std::to_string(cellID) + ")...";

            if (inTerminal)
            {
                std::cout << "\r" << text << std::string(lastLength > text.size() ? lastLength - text.size() : 0, ' ') << std::flush;
                lastLength = text.size();
            }
            else
            {
                std::cout << text << "\n";
            }
        }

        // Ends the line being rewritten, so the next message (or an error) starts on its own line.
        void finish()
        {
            if (inTerminal && lastLength > 0)
            {
                std::cout << "\n";
                lastLength = 0;
            }
        }

    private:
        bool inTerminal;
        std::size_t lastLength = 0;
    };


    // Karyograms of one input file, collected cell by cell, and how they are delivered: in the web viewer (default), as one PNG/SVG file per cell
    // (--save-png/--save-svg), or together in one zip file (--save-zip).
    struct KaryogramOutput
    {
        bool savePNG = false;
        bool saveSVG = false;
        bool saveZip = false;

        std::vector<KaryogramImage> viewerImages;						// Web viewer
        std::vector<std::pair<std::string, std::vector<unsigned char>>> zipEntries;	// --save-zip: {file name, contents}
        int pngFilesWritten = 0;
        int svgFilesWritten = 0;

        bool saveFiles() const { return savePNG || saveSVG || saveZip; }			// Any save option replaces the web viewer
        bool zipPNG() const { return saveZip && (savePNG || !saveSVG); }		// --save-zip holds PNGs unless only --save-svg was asked for
        bool zipSVG() const { return saveZip && saveSVG; }
    };


    bool writeBytesToFile(const std::filesystem::path& path, const std::vector<unsigned char>& bytes)
    {
        std::ofstream file(path, std::ios::binary);
        file.write(reinterpret_cast<const char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
        return static_cast<bool>(file);
    }


    // Saves one cell's karyogram as asked: PNG/SVG files next to the input file, or entries for the zip file. generate draws the cell into the
    // given bytes in the given format. Returns false (after printing what failed) if a karyogram could not be drawn or written.
    bool saveCellKaryogram(
	KaryogramOutput& output,
	const std::filesystem::path& pngPath,						// e.g. <input>_cell17_karyogram_other.png; the SVG uses the same name
	const std::function<bool(std::vector<unsigned char>&, KaryogramFormat)>& generate)
    {
        const bool wantPNG = output.saveZip ? output.zipPNG() : output.savePNG;
        const bool wantSVG = output.saveZip ? output.zipSVG() : output.saveSVG;

        for (const KaryogramFormat format : {KaryogramFormat::PNG, KaryogramFormat::SVG})
        {
            if ((format == KaryogramFormat::PNG && !wantPNG) || (format == KaryogramFormat::SVG && !wantSVG))
            {
                continue;
            }

            const std::filesystem::path path = (format == KaryogramFormat::PNG) ? pngPath : std::filesystem::path(pngPath).replace_extension(".svg");
            std::vector<unsigned char> bytes;

            if (!generate(bytes, format))
            {
                return false;
            }

            if (output.saveZip)
            {
                output.zipEntries.emplace_back(path.filename().string(), std::move(bytes));
            }
            else if (!writeBytesToFile(path, bytes))
            {
                std::cerr << "ERROR: Failed to write karyogram file " << path << "\n";
                return false;
            }
            else
            {
                ++(format == KaryogramFormat::PNG ? output.pngFilesWritten : output.svgFilesWritten);
            }
        }

        return true;
    }


    // Delivers the collected karyograms: opens the web viewer, writes the zip file, or reports where the PNG/SVG files went.
    bool deliverKaryograms(KaryogramOutput& output, const std::filesystem::path& inputPath, const std::filesystem::path& outputFolder, const std::string& genomeType)
    {
        if (!output.saveFiles())
        {
            return output.viewerImages.empty() || showInWebViewer(output.viewerImages, inputPath, outputFolder, genomeType);
        }

        if (output.saveZip)
        {
            if (output.zipEntries.empty())
            {
                return true;
            }

            const std::filesystem::path zipPath = outputFolder / (inputPath.stem().string() + "_karyograms_" + genomeType + ".zip");

            if (!writeZipArchive(output.zipEntries, zipPath))
            {
                return false;
            }

            std::cout << "Karyograms written to: " << zipPath << " (" << output.zipEntries.size() << " files)\n";
            return true;
        }

        const std::filesystem::path folder = std::filesystem::absolute(outputFolder.empty() ? std::filesystem::path(".") : outputFolder).lexically_normal();

        if (output.pngFilesWritten > 0)
        {
            std::cout << "PNG karyograms written to: " << folder << " (" << output.pngFilesWritten << " files)\n";
        }

        if (output.svgFilesWritten > 0)
        {
            std::cout << "SVG karyograms written to: " << folder << " (" << output.svgFilesWritten << " files)\n";
        }

        return true;
    }


    // "17, 24, 28"
    std::string joinCellIDs(const std::vector<int>& cellIDs)
    {
        std::string text;

        for (std::size_t i = 0; i < cellIDs.size(); ++i)
        {
            text += (i > 0 ? ", " : "") + std::to_string(cellIDs[i]);
        }

        return text;
    }


    // Last lines of every run: the cells without karyograms and why, then either the success message or what went wrong.
    // Returns the program's exit code (0 = success).
    int reportFinalStatus(std::size_t cellCount, const std::vector<int>& absentDataCells, const std::vector<int>& failedCells, bool otherError)
    {
        if (!absentDataCells.empty())
        {
            std::cout << "Karyograms not generated for the following cells due to absent SDR data: " << joinCellIDs(absentDataCells) << "\n";
        }

        if (!failedCells.empty())
        {
            std::cerr << "ERROR: Karyograms could not be generated for the following cells: " << joinCellIDs(failedCells) << "\n";
        }

        if (!failedCells.empty() || otherError)
        {
            std::cerr << "Finished with errors, see the messages above.\n";
            return 1;
        }

        std::cout << "All " << cellCount << " cells processed successfully.\n";
        return 0;
    }
}










int main(int argc, char* argv[]) 					// Command line, in this order: --sdd <file> or --sdr <file>, optionally -kh (human) or -ko (other) to draw karyograms, optionally --save-png/--save-svg/--save-zip, optionally the output directory. See printUsage().
{

    // ---------------------------------------
    // Check Command-line arguments
    // ---------------------------------------

    if(argc < 2)							// Must have at least the input type and file.
    {
	std::cerr << "Usage: ./karyoke --sdd <SDD file> [-kh|-ko [--save-png] [--save-svg] [--save-zip]] [<output directory>]\n"
		   << "       ./karyoke --sdr <SDR file> [-kh|-ko [--save-png] [--save-svg] [--save-zip]] [<output directory>]\n"
		   << "       ./karyoke --help\n"
		   << "       ./karyoke --version\n";
	return 1;
    }

    CommandLine commandLine;

    if (!parseCommandLine(argc, argv, commandLine))
    {
        std::cerr << "Run ./karyoke --help for usage.\n";
        return 1;
    }

    if (commandLine.showHelp)
    {
        printUsage();
        return 0;
    }

    if (commandLine.showVersion)
    {
        std::cout << KARYOKE_VERSION_TEXT << "\n";
        return 0;
    }

    std::cout << KARYOKE_VERSION_TEXT << "\n";				// First line of every run, so logs show which version produced the outputs

    const std::string filename = commandLine.filename;
    const bool drawKaryogram = commandLine.drawKaryogram;
    const std::string genomeType = commandLine.genomeType;
    const bool humanGenome = (genomeType == "human");

    KaryogramOutput output;
    output.savePNG = commandLine.savePNG;
    output.saveSVG = commandLine.saveSVG;
    output.saveZip = commandLine.saveZip;

    // Summary file and karyograms go to the output directory if one was given (created if needed), else next to the input file.
    std::filesystem::path outputFolder = std::filesystem::path(filename).parent_path();

    if (!commandLine.outputDirectory.empty())
    {
        outputFolder = commandLine.outputDirectory;
        std::error_code error;
        std::filesystem::create_directories(outputFolder, error);		// No error if the folder already exists

        if (error || !std::filesystem::is_directory(outputFolder, error))
        {
            std::cerr << "ERROR: Could not create the output directory " << outputFolder << (error ? " (" + error.message() + ")" : std::string()) << "\n";
            return 1;
        }
    }

    // --------------------------------------
    // SDR MODE
    // --------------------------------------

    if (commandLine.sdrInput)
    {
        std::cout << "SDR file processing: " << filename << "...\n";

        SDRparser parser;

        if (!parser.parseFile(filename))
        {
            std::cerr << "ERROR: Failed to load SDR file " << filename << "\n";
            return 1;
        }

	std::filesystem::path inputPath(filename);
        const std::vector<SDRsubHeader>& subHeaders = parser.getSubHeaders();

        std::vector<int> absentDataCells;				// Cells with no SDR data records: nothing to draw, not an error
        std::vector<int> failedCells;					// Cells whose karyogram could not be drawn or written: an error

	// --------------------------------------------- //
	// SDR karyogram drawing, one cell at a time
	// --------------------------------------------- //
        if (drawKaryogram)						// If user specifies -kh/-ko then draw
        {
	    Karyogram karyogram;					// Instantiate the Karyogram object karyogram
	    ProgressLine progress;

	    for (std::size_t i = 0; i < subHeaders.size(); ++i)		// Draw one karyogram per cell subheader present in SDR file
	    {
		const SDRsubHeader& subHeader = subHeaders[i];
		progress.show(i + 1, subHeaders.size(), subHeader.cellID);

		// A cell without data records has nothing to draw. It is listed at the end rather than reported as an error,
		// since the file itself is fine (generateSDRkaryogram() would otherwise refuse it with an error).
		if (subHeader.dataRecords.empty())
		{
		    absentDataCells.push_back(subHeader.cellID);
		    continue;
		}

		std::filesystem::path cellKaryogramPath =			// Modify name of output karyogram
            		outputFolder /
            		(inputPath.stem().string() + "_cell" + std::to_string(subHeader.cellID) +
             		"_karyogram_" + genomeType + ".png");

		const auto generate = [&](std::vector<unsigned char>& bytes, KaryogramFormat format)
		{
		    return karyogram.generateSDRkaryogram(parser.getMasterHeader(), subHeader, humanGenome, bytes, format);
		};

		if (!output.saveFiles())					// Default: generate in memory and show in the web viewer
		{
		    KaryogramImage image;
		    image.label = "Cell " + std::to_string(subHeader.cellID);
		    image.filename = cellKaryogramPath.filename().string();

		    if (!generate(image.svgBytes, KaryogramFormat::SVG)	// SVG stays sharp when zooming in the viewer
			|| !generate(image.pngBytes, KaryogramFormat::PNG))	// PNG for the viewer's PNG download, identical to --save-png
		    {
			failedCells.push_back(subHeader.cellID);
			continue;
		    }

		    image.hotspots = buildLabelHotspots(karyogram.getLabelRegions(),					// Clickable chromosome labels
			countMutationsPerOriginalStrand(subHeader, parser.getMasterHeader()));
		    image.popupNote = "Events involving several chromosomes are counted on each of them.";
		    output.viewerImages.push_back(std::move(image));
		    continue;
		}

		if (!saveCellKaryogram(output, cellKaryogramPath, generate))	// --save-png / --save-svg / --save-zip
		{
		    failedCells.push_back(subHeader.cellID);
		}
	    }

	    progress.finish();
        }
        else
        {
            std::cout << "Summarizing " << subHeaders.size() << " cells...\n";
        }

	// --------------------------------------
	// SDR summary file
	// --------------------------------------
        std::filesystem::path summaryPath = outputFolder /
            (inputPath.stem().string() + "_summary.txt");

        bool otherError = false;

        if (!parser.writeSummary(summaryPath.string()))
        {
            std::cerr << "ERROR: Failed to create SDR summary file " << summaryPath << "\n";
            otherError = true;
        }
        else
        {
            std::cout << "Summary written to: " << summaryPath << "\n";
        }

        if (drawKaryogram && !deliverKaryograms(output, inputPath, outputFolder, genomeType))
        {
            otherError = true;
        }

        return reportFinalStatus(subHeaders.size(), absentDataCells, failedCells, otherError);
    }


    // --------------------------------------
    // SDD MODE
    // --------------------------------------

    // -------------------------------------
    // Load SDD file
    // -------------------------------------

    std::cout << "SDD file processing: " << filename << "...\n";

    SDDparser parser; 							// Instantiate SDDparser object called parser

    if(!parser.load(filename))						// Check if SDD input file loads correctly
    {
        std::cerr << "ERROR: Failed to load SDD file " << filename << "\n";
        return 1;
    }

    std::filesystem::path inputPath(filename);				// Names the outputs; they go to outputFolder
    const std::vector<Exposure>& exposures = parser.getExposures();

    std::vector<int> failedCells;					// Cells whose karyogram could not be drawn or written: an error

    // ---------------------------------------------
    // SDD karyogram drawing, one cell at a time
    // ---------------------------------------------

    if (drawKaryogram)							// Check if user passed '-kh'/'-ko', meaning drawKaryogram == true.
    {
        const Header& header = parser.getHeader();			// Get header to pass the chromosome sizes and cell cycle phase to the generate Karyogram function.

        Karyogram karyogram;						// Instantiate the Karyogram object called karyogram.
        ProgressLine progress;

    	for (std::size_t i = 0; i < exposures.size(); ++i)		// Draw a karyogram for each individual cell exposure in a given SDD file
    	{
    	    const Exposure& exposure = exposures[i];
    	    const int cellID = exposure.exposureID - 1;		// Cells are numbered from 0, exposure IDs from 1
    	    progress.show(i + 1, exposures.size(), cellID);

    	    // Create a vector containing only this exposure.
    	    std::vector<Exposure> singleExposure;
    	    singleExposure.push_back(exposure);

    	    // Create a separate output filename for this cell.
    	    std::filesystem::path karyogramPath = outputFolder /
            (inputPath.stem().string() + "_cell" +
	    std::to_string(cellID) + "_karyogram_" +
	    genomeType + ".png");

    	    const auto generate = [&](std::vector<unsigned char>& bytes, KaryogramFormat format)
    	    {
    	        return karyogram.generateSDDkaryogram(
            	    header.chromosome_sizes,
            	    header.cell_cycle_phase,
            	    header.dose_or_fluence,
            	    header.incident_particles,
            	    singleExposure,
            	    humanGenome,
            	    bytes,
            	    format);
    	    };

    	    if (!output.saveFiles())					// Default: generate in memory and show in the web viewer
    	    {
    	        KaryogramImage image;
    	        image.label = "Cell " + std::to_string(cellID);
    	        image.filename = karyogramPath.filename().string();

    	        if (!generate(image.svgBytes, KaryogramFormat::SVG)		// SVG stays sharp when zooming in the viewer
    	            || !generate(image.pngBytes, KaryogramFormat::PNG))	// PNG for the viewer's PNG download, identical to --save-png
    	        {
            	    failedCells.push_back(cellID);
            	    continue;
    	        }

    	        image.hotspots = buildLabelHotspots(karyogram.getLabelRegions(),			// Clickable chromosome labels
    	            countDamagesPerChromosome(parser, exposure));
    	        image.popupEmptyText = "No damage detected";					// SDD files record damage, not mutations
    	        output.viewerImages.push_back(std::move(image));
    	        continue;
    	    }

    	    if (!saveCellKaryogram(output, karyogramPath, generate))	// --save-png / --save-svg / --save-zip
    	    {
    	        failedCells.push_back(cellID);
    	    }
    	}

    	progress.finish();
    }
    else
    {
        std::cout << "Summarizing " << exposures.size() << " cells...\n";
    }

    // -------------------------------------
    // Create summary file
    // -------------------------------------

    std::filesystem::path summaryPath = outputFolder / 		// Summary file in the output directory (default: next to the input SDD file)
	(inputPath.stem().string() + "_summary.txt");

    bool otherError = false;
    std::ofstream summaryFile(summaryPath);				// Output to desired summary file path

    if(!summaryFile)							// Make sure summary file created successfully.
    {
        std::cerr << "ERROR: Failed to create summary file " << summaryPath << "\n";
        otherError = true;
    }
    else
    {
        parser.printSummary(summaryFile); 				// Returns Summary of SDD header fields in a summary file at the desired path
        summaryFile.close();
        std::cout << "Summary written to: " << summaryPath << "\n";	// Let user know summary file was created successfully
    }

    if (drawKaryogram && !deliverKaryograms(output, inputPath, outputFolder, genomeType))
    {
        otherError = true;
    }

    return reportFinalStatus(exposures.size(), {}, failedCells, otherError);
}
