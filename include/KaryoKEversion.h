#ifndef KARYOKE_VERSION_H
#define KARYOKE_VERSION_H

// KaryoKE version, shown by "karyoke -v" / "--version", at the start of every run, in the summary files and in the web viewer.
// Change it here only; everything else reads it from this header.
#define KARYOKE_VERSION "1.1"

namespace karyoke
{

inline constexpr const char* KARYOKE_VERSION_TEXT = "KaryoKE v" KARYOKE_VERSION;

} // namespace karyoke

#endif
