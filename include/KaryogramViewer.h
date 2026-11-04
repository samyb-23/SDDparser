#ifndef KARYOGRAM_VIEWER_H
#define KARYOGRAM_VIEWER_H

#include <string>
#include <vector>
#include <filesystem>
#include <utility>

namespace karyoke
{

struct KaryogramPopupSection                    // One part of a chromosome's popup, e.g. "Homolog 1".
{
    std::string heading;                        // Shown above the counts; may be empty when the chromosome has only one homolog.
    std::vector<std::pair<std::string, int>> counts;   // Mutation type and count, only the types that are present.
};

struct KaryogramHotspot                         // A clickable area of the karyogram (a chromosome label) that opens a popup.
{
    std::string title;                          // Popup title, e.g. "Chromosome 1".
    double x, y, width, height;                 // Area in image pixels, from Karyogram::getLabelRegions().
    std::vector<KaryogramPopupSection> sections;
};

struct KaryogramImage                           // One karyogram to show in the web viewer.
{
    std::string label;                          // Shown in the cell list, e.g. "Cell 0".
    std::string filename;                       // Name used for downloads from the viewer, e.g. "cell0_karyogram.png" (the extension is swapped to .svg for SVG downloads).
    std::vector<unsigned char> svgBytes;        // Vector image (KaryogramFormat::SVG): shown in the viewer, stays sharp at any zoom level, offered as SVG download.
    std::vector<unsigned char> pngBytes;        // Raster image (KaryogramFormat::PNG): offered as PNG download, and shown only when svgBytes is empty.
                                                // Either one may be left empty; the viewer then only offers the format it has.
    std::vector<KaryogramHotspot> hotspots;     // Clickable chromosome labels with mutation or damage summaries (may be empty).
    std::string popupNote;                      // Small print at the bottom of each popup (may be empty).
    std::string popupEmptyText = "No mutations detected";   // Shown for a homolog with nothing to list, e.g. "No damage detected" for SDD.
};


// Writes a single self-contained HTML file with every karyogram embedded in it (no server or extra files needed).
bool writeKaryogramViewer(
    const std::vector<KaryogramImage>& images,  // Karyograms to show, in order.
    const std::string& title,                   // Page title, e.g. the input file name.
    const std::filesystem::path& htmlPath       // Where to write the HTML file.
);

// Writes several files (e.g. every karyogram of an input file) into one .zip archive, stored without compression since PNG and
// SVG output is already compact enough. Entries are {file name inside the zip, file contents}. Returns false if it could not write.
bool writeZipArchive(
    const std::vector<std::pair<std::string, std::vector<unsigned char>>>& entries,
    const std::filesystem::path& zipPath
);

// Returns false when no browser can be shown on this machine, e.g. over SSH, on an HPC cluster or in a container,
// or in WSL with Windows interop switched off. Used to decide where to write the viewer and whether to open it.
bool canOpenBrowser();

// Opens a file in the system's default web browser (xdg-open on Linux, open on macOS, start on Windows, the Windows
// browser from WSL). Returns false if it could not, in which case the caller should tell the user where the file is.
bool openInBrowser(const std::filesystem::path& htmlPath);

} // namespace karyoke

#endif
