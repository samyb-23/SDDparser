#ifndef SDD_TYPES_H
#define SDD_TYPES_H

#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>

namespace karyoke
{

// Use for checking if any optional header fields are absent, to not raise an error, but report as not measured.
constexpr int SDD_INT_FIELD_NOT_MEASURED = -1;
constexpr double SDD_DOUBLE_FIELD_NOT_MEASURED = -1.0;

// The SDD header may contain any of the following header fields.
// If header field provided is not in this list, exit with error.
struct Header
{
    std::string sdd_version;
    std::string software;
    std::string author;
    std::string simulation_details; 
    std::string source;	
    int source_type = SDD_INT_FIELD_NOT_MEASURED; 
    std::vector<int> incident_particles; 
    std::vector<double> mean_particle_energy; 
    std::string energy_distribution; 
    std::vector<double> particle_fraction; 
    std::vector<double> dose_or_fluence; 
    double dose_rate = SDD_DOUBLE_FIELD_NOT_MEASURED;
    std::string irradiation_target;
    std::vector<double> volumes; 
    std::vector<double> chromosome_sizes;
    double DNA_density = SDD_DOUBLE_FIELD_NOT_MEASURED; 
    std::vector<double> cell_cycle_phase; 
    std::vector<int> DNA_structure; 
    int in_vitro_or_in_vivo = SDD_INT_FIELD_NOT_MEASURED; 
    std::vector<std::string> proliferation_status;
    std::vector<double> microenvironment;
    std::vector<double> damage_definition;
    double time = SDD_DOUBLE_FIELD_NOT_MEASURED;
    std::vector<int> damage_and_primary_count;
    std::vector<int> data_entries;
    std::string additional_information;

};


// The following structures follow the SDD data field structure specifications

//--------------------------------------------------
// Data Field Entry 1 - Classification
//--------------------------------------------------
struct Classification
{
    int exposureMarker = 0; 		// 2 for new exposure, 1 for new primary particle from same exposure, 0 for same primary particle as previous row
    int eventID = 0;	    		// The ID number of the primary particle that was simulated
};

//--------------------------------------------------
// Data Entry 2 - Spatial Position
//--------------------------------------------------
struct Position 			// Records the center and extent of each recorded damage within the bounding box specified in Header Field 14 (in micrometers)
{
    double x = 0.0;
    double x_min = 0.0; 
    double x_max = 0.0;
    double y = 0.0;
    double y_min = 0.0;
    double y_max = 0.0;
    double z = 0.0;
    double z_min = 0.0;
    double z_max = 0.0;
};

//--------------------------------------------------
// Data Entry 3 - Chromosome IDs
//--------------------------------------------------
struct ChromosomeID 			// Stores the identity of the chromatid where the damage occurred 
{
    int dnaStructure = 0; 		// 0 = unspecified, 1 = heterochromatin, 2 = euchromatin, 3 = free DNA, 4 = mtDNA/bacterial DNA/viral DNA.
    int chromosomeNumber = 0; 		// 1-46 for human chromosome
    int chromatidNumber = 0; 		// 1 for unduplicated chromosomes, 1 or 2 for duplicated chromosomes in the S/G2/M phases
    int chromosomeArm = 0; 		// 0 for short p-arm, 1 for long q-arm.
};

//--------------------------------------------------
// Data Entry 4 - Chromosome Position
//--------------------------------------------------
struct ChromosomePosition 		// Indicates damage position along the chromosome's genetic length defined from the start of the short p arm to the end of the long q arm
{
    double position = 0.0; 		// Stored as either a fraction of the total length between 0 and 1, or the length in base pairs if greater than or equal to 1. 
    bool isFractional = false; 		// Flag to check if position entered is in base pairs or fractional.
};

//--------------------------------------------------
// Data Entry 5 - Cause
//--------------------------------------------------
struct DamageCause 
{
    int cause = 0; 			// 0 = direct physical damage, 1 = indirect chemical damage, 2 = combination of direct + indirect, 3 = charge migration damage
    int numDirectDamages = 0; 		// Counts the number of direct damages
    int numIndirectDamages = 0; 	// Counts the number of indirect damages
};

//--------------------------------------------------
// Data Entry 6 - Damage Types
//--------------------------------------------------
struct DamageType			// Specifies the type of damage present at a given site. Damages separated by less than the minimum distance of base pairs specified in header field 22 are scored in a single data block.
{
    int numBaseDamages = 0;		
    int numSingleBackboneBreaks = 0;
    int presenceOfDoubleStrandBreaks = 0; // Binary flag 0 for no, 1 for yes.
};

//--------------------------------------------------
// Data Entry 7 - Full break spec
//--------------------------------------------------
struct FullBreakSpec
{
    std::vector<int> strand;		// 1 = 5' to 3' backbone, 2 = 5' to 3' bases, 3 = 3' to 5' bases, 4 = 3' to 5' backbone

    std::vector<int> base;		// Identifies the DNA base where damage occurs

    std::vector<int> baseDamageType; 	// 0 = no damage, 1 = direct damage, 2 = indirect damage, 3 = direct + indirect damage
};


// Exposure entry structure, need to include optional fields 8-14 later.
struct DamageEntry
{
    Classification classification;
    Position position;
    ChromosomeID chromosomeID;
    ChromosomePosition chromosomePosition;
    DamageCause damageCause;
    DamageType damageType;
    FullBreakSpec fullBreakSpec;

    std::string rawLine;

};

// Exposure structure for formatting summary file according to the damages from each exposure ID
struct Exposure
{
    int exposureID = 0;
    std::vector<DamageEntry> damages;
};


//--------------------------------------------------
// Chromosome damage summary
//--------------------------------------------------
struct ChromosomeDamageSummary // Important damage information to be summarized for each chromosome, used to format summary file output.
{
    int dnaStructure = 0;
    int chromosomeNumber = 0;
    int numBaseDamages = 0;
    int numSingleStrandBreaks = 0;
    int numDoubleStrandBreaks = 0;
};

//-------------------------------------------------
// Plotting double-strand/ single-strand break locations 
//-------------------------------------------------
struct DamageLocation		// Stores the chromosome position in base pairs where the damage occurs, converted to x and y coordinates for karyogram.
{
    ChromosomeID chromosomeID;
    ChromosomePosition chromosomePosition;
    int numSingleStrandBreaks;
};

} // namespace karyoke

#endif
