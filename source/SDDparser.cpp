#include <fstream>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>

#include "SDDparser.h"
#include "generalUtilities.h"
#include "SDDparserVersion.h"
#include "SDDutilities.h"




namespace sddparser
{

// Function to check if SDD file successfully loads, and if header and data blocks were successfully processed.
bool SDDparser::load(
    const std::string& filename)		// Takes as input the SDD file path
{
    std::ifstream file(filename);		// Create input file stream object file using the input filename.

    if(!file)					// Check if SDD file loads correctly.
        { 
	 std::cerr << "ERROR: SDD header parsing failed: " << filename << std::endl;
	 return false;
	}

    header = Header();				// Clear header and exposure objects prior to loading file
    exposures.clear();

    if(!parseHeader(file))			// Check if SDD header can be parsed
    {
    	std::cerr << "ERROR: SDD header parsing failed." << std::endl;
        return false;
    }

    if(!parseDamageEntries(file))		// Check if SDD data block can be parsed
    {
    	std::cerr << "ERROR: SDD data entry parsing failed.\n";
    	return false;
    }

    return true;				// If all checks pass, return true
}

// List all valid SDD header field names. Check to see if invalid SDD header field name has been encountered.
const std::unordered_set<std::string> validSDDHeaderFields =
{
    "sdd version",
    "software",
    "author",
    "simulation details",
    "source",
    "source type",
    "incident particles",
    "mean particle energy",
    "energy distribution",
    "particle fraction",
    "dose or fluence",
    "dose rate",
    "irradiation target",
    "volumes",
    "chromosome sizes",
    "dna density",
    "cell cycle phase",
    "dna structure",
    "in vitro / in vivo",
    "proliferation status",
    "microenvironment",
    "damage definition",
    "time",
    "damage and primary count",
    "data entries",
    "additional information"
};

// Function to store the values associated with the various SDD header fields, to be summarized later in printHeaderSummary
bool SDDparser::parseHeader(
    std::ifstream& file)							// Take as input SDD file
{
    std::string line;								// Create line variable

    while(std::getline(file,line))						// Loop through each line in the SDD file
    {
        line = trim(line);							// Remove all white spaces, convert line to a vector of strings

        if(line.empty())							// If line empty, check next line
        {
	    continue;
	}
        
        if(line.back()==';')							// If line ends in ';'
        {
	    line.pop_back();							// Remove last element of the line vector ';'
	}

        if(normalizeHeaderKey(line) == "***endofheader***")			// If end of header block reached, give success message
        {    std::cout << "SDD header parsed successfully.\n";
		
	     return true;
	}

        std::vector<std::string> tokens = split(line,',');			// Split header entries by commas, creating a vector of strings

        if(tokens.empty())			
        {
	    continue;
	}

        std::string original_key = tokens[0]; 					// Header field names before removing whitespaces and capitalizations
	std::string key = normalizeHeaderKey(original_key); 			// Header field names are now lower case and have no leading or trailing whitespaces

	// Check if SDD header field is in the above list of supported fields
	if(validSDDHeaderFields.find(key) == validSDDHeaderFields.end())
	{
    	    std::cerr
            << "ERROR: Unknown SDD header field: " << original_key << std::endl;
    	    return false;
	}

        std::vector<std::string> values(tokens.begin()+1, tokens.end()); 	// Check next header field 

        //----------------------------
        // Map SDD header fields
        //----------------------------

	// Fields with different numbering data types are converted to floats.
	// Fields with strings and ints/floats are kept as strings and handled later.

        if(key=="sdd version")
        {
            header.sdd_version = values.empty() ? "" : values[0]; 
        }

        else if(key=="software")
        {
            header.software = values.empty() ? "" : values[0];
        }

        else if(key=="author")
        {
            header.author = values.empty() ? "" : values[0];
        }

	else if(key=="simulation details")
	{
	    header.simulation_details = values.empty() ? "" : values[0]; 
	}

	else if(key=="source")
	{
	    header.source = values.empty() ? "" : values[0];
	}

	else if(key=="source type")
	{
	    if (values.empty())
	    {
		header.source_type = SDD_INT_FIELD_NOT_MEASURED;
	    }
	    else
	    {
		try
		{
	    	    header.source_type = std::stoi(values[0]);
		}
		catch (const std::exception& error)
		{
		    std::cerr << "ERROR: Invalid 'Source type' value in SDD header: " << values[0] << " (" << error.what() << ")\n";
            	    return false;
		}
	    }
	}

	else if(key=="incident particles")
	{
	    if (values.empty())
	    {
		header.incident_particles.clear();
	    }
	    else
	    {
		try
	    	{
	            header.incident_particles = parseIntList(values); 
	    	}
	    	catch (const std::exception& error)
	    	{
		    std::cerr << "ERROR: Invalid 'Incident particles' value in SDD header (" << error.what() << ")\n";
		    return false;
	    	}
	    }
	}

	else if(key=="mean particle energy")
	{
	    if (values.empty())
	    {
	    	header.mean_particle_energy.clear();
	    }
	    else
	    {
		try
		{
		    header.mean_particle_energy = parseDoubleList(values);
		}
		catch (const std::exception& error)
		{
		    std::cerr << "ERROR: Invalid 'Mean particle energy' value in SDD header (" << error.what() << ")\n";
		    return false;
		}
	    }
	}

	else if(key=="energy distribution")
	{
	    header.energy_distribution = values.empty() ? "" : values[0];
	}

	else if(key=="particle fraction")
	{
	    if (values.empty())
	    {
		header.particle_fraction.clear();
	    }
	    else
	    {
	    	try
	    	{
	    	    header.particle_fraction = parseDoubleList(values);
	    	}
	    	catch (const std::exception& error)
	    	{
		    std::cerr << "ERROR: Invalid 'Particle fraction' value in SDD header (" << error.what() << ")\n";
		    return false;
	    	}
	    }
	}

        else if(key=="dose or fluence")
        {
	    if (values.empty())
	    {
		header.dose_or_fluence.clear();
	    }
	    else
	    {
            	try
		{
		    header.dose_or_fluence = parseDoubleList(values);
		}
		catch (const std::exception& error)
		{
		    std::cerr << "ERROR: Invalid 'Dose or fluence' value in SDD header (" << error.what() << ")\n";
		    return false;
		}
	    }
	}

        else if(key=="dose rate")
        {
	    if (values.empty())
	    {
		header.dose_rate = SDD_DOUBLE_FIELD_NOT_MEASURED;
	    }
	    else
	    {
		try
		{
            	    header.dose_rate = std::stod(values[0]);
		}
		catch (const std::exception& error)
		{
		    std::cerr << "ERROR: Invalid 'Dose rate' value in SDD header (" << error.what() << ")\n";
		    return false;
		}
	    }
        }

	else if(key=="irradiation target")
	{
	    header.irradiation_target = values.empty() ? "" : values[0];
	}

	else if(key=="volumes")
	{
	    if (values.empty())
	    {
		header.volumes.clear();
	    }
	    else
	    {
		try
		{
		    header.volumes = parseDoubleList(values);
		}
		catch (const std::exception& error)
		{
		    std::cerr << "ERROR: Invalid 'Irradiation target' value in SDD header (" << error.what() << ")\n";
		    return false;
		}
	    }
	}

        else if(key=="chromosome sizes")		// REQUIRED ENTRY
        {
	    if (values.empty())
	    {
		std::cerr << "ERROR: 'Chromosome sizes' is a required SDD header field and cannot be left blank.\n";
		return false;
	    }
	    else
	    {
	        try
		{
		    header.chromosome_sizes = parseDoubleList(values);
		}
		catch (const std::exception& error)
		{
		    std::cerr << "ERROR: Invalid 'Chromosome sizes' value in SDD header (" << error.what() << ")\n";
		    return false;
		}
	    }
        }

	else if(key=="dna density")
	{
	    if (values.empty())
    	    {
        	header.DNA_density = SDD_DOUBLE_FIELD_NOT_MEASURED;
    	    }
    	    else
    	    {
        	try
        	{
            	    header.DNA_density = std::stod(values[0]);
        	}
        	catch (const std::exception& error)
        	{
           	     std::cerr << "ERROR: Invalid 'DNA density' value in SDD header: " << values[0] << " (" << error.what() << ")\n";
            	     return false;
        	}
    	    }
	}

	else if(key=="cell cycle phase")		// REQUIRED
	{
	    if (values.empty())
    	    {
        	std::cerr << "ERROR: 'Cell cycle phase' is a required SDD header field and cannot be blank.\n";
        	return false;
    	    }
	    else
	    {
    		try
    		{
        	    header.cell_cycle_phase = parseDoubleList(values);
    		}
    		catch (const std::exception& error)
    		{
       	 	    std::cerr << "ERROR: Invalid 'Cell cycle phase' value in SDD header (" << error.what() << ")\n";
        	    return false;
    		}
	    }
	}

	else if(key=="dna structure")
	{

	    if (values.empty())
	    {
		header.DNA_structure.clear();
	    }

	    try
    	    {
        	header.DNA_structure = parseIntList(values);
    	    }
    	    catch (const std::exception& error)
    	    {
        	std::cerr << "ERROR: Invalid 'DNA structure' value in SDD header (" << error.what() << ")\n";
        	return false;
    	    }

	}

	else if(key=="in vitro / in vivo")
	{

	    if (values.empty())
    	    {
        	header.in_vitro_or_in_vivo = SDD_INT_FIELD_NOT_MEASURED;
    	    }
    	    else
    	    {
        	try
        	{
            	    header.in_vitro_or_in_vivo = std::stoi(values[0]);
        	}
        	catch (const std::exception& error)
        	{
            	    std::cerr << "ERROR: Invalid 'In vitro / in vivo' value in SDD header: " << values[0] << " (" << error.what() << ")\n";
            	    return false;
        	}
    	    }
	}

	else if(key=="proliferation status")
	{
	   if (!values.empty())
    	   {
        	header.proliferation_status.push_back(values[0]);
    	   }
	}

	else if(key=="microenvironment")
	{
	    if (values.empty())
	    {
		header.microenvironment.clear();
	    }
	    else
	    {
		try
    		{
        	    header.microenvironment = parseDoubleList(values);
    		}
    		catch (const std::exception& error)
    		{
        	    std::cerr << "ERROR: Invalid 'Microenvironment' value in SDD header (" << error.what() << ")\n";
        	    return false;
    		}
	    }
	}

	else if(key=="damage definition")
	{
	    if (values.empty())
	    {
		header.damage_definition.clear();
	    }
	    else
	    {
		try
		{
		    header.damage_definition = parseDoubleList(values);
		}
		catch (const std::exception& error)
		{
		    std::cerr << "ERROR: Invalid 'Damage definition' value in SDD header (" << error.what() << ")\n";
		    return false;
		}
	    }
	}

	else if(key=="time")
	{
	    if (values.empty())
    	    {
        	header.time = SDD_DOUBLE_FIELD_NOT_MEASURED;
    	    }
    	    else
    	    {
        	try
        	{
            	    header.time = std::stod(values[0]);
        	}
        	catch (const std::exception& error)
        	{
            	    std::cerr << "ERROR: Invalid 'Time' value in SDD header: " << values[0] << " (" << error.what() << ")\n";
            	    return false;
        	}
    	    }
	}

	else if(key=="damage and primary count")
	{
	    if (values.empty())
	    {
		header.damage_and_primary_count.clear();
	    }
	    else
	    {
		try
		{
		    header.damage_and_primary_count = parseIntList(values);
		}
		catch (const std::exception& error)
		{
		    std::cerr << "ERROR: Invalid 'Damage and primary count' value in SDD header (" << error.what() << ")\n";
		    return false;
		}
	    }
	}

	else if(key=="data entries")			// REQUIRED TO INTERPRET SDD DATA BLOCK
	{
	    if (values.empty())
	    {
		std::cerr << "Data entries must be filled with 14 comma separated values of either 0 or 1.\n";
		return false;
	    }
	    else
	    {
		try
		{
		    header.data_entries = parseIntList(values);
		}
		catch (const std::exception& error)
		{
		    std::cerr << "ERROR: Invalid 'Data entries' value in SDD header (" << error.what() << ")\n";
		    return false;
		}
	    }
	}

	else if(key=="additional information")
	{
	    header.additional_information = values.empty() ? "" : values[0];
	}

        else
        {
            std::cerr << "Unknown SDD header field entry = " << key << "\n"; // If specified header name unknown
        }

    }

	std::cerr << "ERROR: End of SDD header not found.\n"; // If ***EndOfHeader*** not found, exit with error
	return false;
}


















// Function to store the different data field entries to be summarized in printExposureSummary
bool SDDparser::parseDamageEntries(
    std::ifstream& file)			// Takes as input SDD file
{

    std::string line;

    int currentExposureID = 0;			// Counter to record number of exposures stored in the SDD file.

    // Determine expected number of fields. 
    std::size_t expectedFields = 0;		// Used to align SDD header 'Data entries' field with the associated data fields. 

    for(int enabled : header.data_entries)	// Check if the entries in 'Data entries' header = 1, indicating that field is used in the data block.
    {						// Increment the number of expected fields if a 1 is encountered.
	if(enabled == 1)
	{
        expectedFields++;
	}
    }

    while(std::getline(file,line))		// Loop through the lines of the SDD file
    {
	line = trim(line);			// Remove whitespaces

        if(line.empty())			// If line is empty, check next line.
	{ 
           continue;
	}

        DamageEntry damage;			// Instantiate DamageEntry object called damage

        damage.rawLine = line;			// Use to store the raw string of the data entry

	//--------------------------------------------------------------------------------
        // Split the data entry into Fields 1-7 -- Will do Optional fields 8-14 in future.
        //--------------------------------------------------------------------------------

        std::vector<std::string> fields = split(line, ';');			// Split data entries by the field delimiter ';'

        if(fields.size() != expectedFields) 					// If statement to check if the SDD header 'Data entries' number of true fields equals the number of fields in the
	{			            					// data section. If not, exit with error.
            std::cerr << "ERROR: SDD data entry contains " << fields.size() << " fields, but the header specifies " << expectedFields << " fields.\n";
    	    return false;
	}

	std::size_t fieldIndex = 0;	    					// Variable to store the correct index and thus data field used in the data block according to the 'Data entries' header.

	//---------------------------------------------------------------------//
	// Validating which header fields are mandatory and which are optional //
	//---------------------------------------------------------------------//

	bool hasField1 = header.data_entries.size() > 0 && header.data_entries[0] == 1;

	bool hasField2 = header.data_entries.size() > 1 && header.data_entries[1] == 1;

	bool hasField3 = header.data_entries.size() > 2 && header.data_entries[2] == 1;

	bool hasField4 = header.data_entries.size() > 3 && header.data_entries[3] == 1;

	bool hasField5 = header.data_entries.size() > 4 && header.data_entries[4] == 1;

	bool hasField6 = header.data_entries.size() > 5 && header.data_entries[5] == 1;

	bool hasField7 = header.data_entries.size() > 6 && header.data_entries[6] == 1;

	if(!hasField1)								// Field 1 is mandatory
    	{
            std::cerr << "ERROR: Field 1 (Classification) is mandatory.\n";
            return false;
    	}

	if(!hasField2 && !(hasField3 && hasField4))				// Either Field 2 OR both Fields 3 and 4 must be present.
	{
    	    std::cerr << "ERROR: SDD data entries must contain either Field 2 (Spatial Position) or both "
        	      << "Field 3 (Chromosome IDs) and Field 4 (Chromosome Position).\n";
    	    return false;
	}

	if(!hasField6 && !hasField7)						// Either Field 6 or Field 7 must be present.
	{
    	    std::cerr << "ERROR: SDD data entries must contain either Field 6 (Damage Types) or Field 7 "
        	  << "(Full Break Specification) or both.\n";
    	    return false;
	}

	//--------------------------------------------------
        // Field 1: Classification
        //--------------------------------------------------

        if(hasField1) 								// Check if first 'Data entries' field is = 1.
        {

 	    if(fieldIndex >= fields.size())					// Check if mandatory Data Field 1 entry is missing.
    	    {
        	std::cerr << "ERROR: Missing SDD data entry Field 1.\n";
        	return false;
    	    }

	    std::vector<std::string> values = split(fields[fieldIndex], ',');

            if(values.size() >= 1 && !trim(values[0]).empty())
            {

	        int exposureMarker = std::stoi(trim(values[0]));		// This field marks the start of a new exposure
		damage.classification.exposureMarker = exposureMarker;

        	if(exposureMarker == 2 || currentExposureID == 0) 		// A value of 2 indicates the beginning of a new exposure. If SDD data field 1
		{						  		// does not start with a 2, force it to be considered a new exposure anyway.      	
            	    currentExposureID++;					// Exposure ID starts at 1 and increases each time exposureMarker = 2.

            	    Exposure newExposure;					// Instantiate the Exposure object called newExposure.
            	    newExposure.exposureID = currentExposureID;			// Associate the newExposure object's ID with the currentExposureID

            	    exposures.push_back(newExposure);				// Add newExposure object to the vector of exposures to be summarized in printExposureSummary
        	}
	    }

            if(values.size() >= 2 && !trim(values[1]).empty())			// Check if optional second entry 'event ID' is present, if so, store
            {
	        damage.classification.eventID = std::stoi(trim(values[1]));
	    }

	    fieldIndex++;							// After data Field 1 handled, increase the index to check next field and compare with 'Data entries' header field.
        }


	//--------------------------------------------------
        // Field 2: Spatial Position
        //--------------------------------------------------

        if(hasField2)								// Check if second data field present
        {

 	    if(fieldIndex >= fields.size())					// Check for mismatch between data field number and 'Data entries' fields
    	    {
        	std::cerr << "ERROR: Missing SDD data entry Field 2.\n";
        	return false;
    	    }

	    std::vector<std::string> fieldEntries;				// Outer vector to store the entire field, which can either be three comma separated values, or 3 groups of 3 comma separated values, each group delimited by '/'

	    std::vector<std::string> subfields = split(fields[fieldIndex], '/');// If subfields of Field 2 are present, split by '/', handle each subfield
	    if (subfields.size() == 3)						// Subfield size must be 3
	    {
		for (const std::string& subfield: subfields)		
		{
       	            std::vector<std::string> values = split(subfield, ',');
		    fieldEntries.insert(fieldEntries.end(), values.begin(), values.end()); // fieldEntries vector just becomes a vector of 9 values, removing the delimiters.
	        }

		if (fieldEntries.size() != 9)					// If subfields were present, there must be 9 values in Field 2.
		{
		    std::cerr << "ERROR: Incorrect number of entries for Data field 2.\n";
    		    return false;
		}

		else								// If correct number of elements, associate the values to their designated variable
		{
		damage.position.x = std::stod(trim(fieldEntries[0]));		// x, y, and z describe the centers of the damages
		damage.position.y = std::stod(trim(fieldEntries[1]));
		damage.position.z = std::stod(trim(fieldEntries[2]));
         	damage.position.x_max = std::stod(trim(fieldEntries[3])); 	// the next six variables describe the extent of the damages from the x, y, and z coordinates.
		damage.position.y_max = std::stod(trim(fieldEntries[4]));
		damage.position.z_max = std::stod(trim(fieldEntries[5]));
         	damage.position.x_min = std::stod(trim(fieldEntries[6]));
		damage.position.y_min = std::stod(trim(fieldEntries[7]));
		damage.position.z_min = std::stod(trim(fieldEntries[8]));
         	}

	    }

	    else if (subfields.size() == 1)					// If no subfields, the values are delimited by ','
	    {
		std::vector<std::string> values = split(fields[fieldIndex], ',');
		if (values.size() == 3)						// Make sure exactly three values are in the field if only on subfield
		{
		    damage.position.x = std::stod(trim(values[0]));		// Each value describes the center of the damages in the x, y, and z coordinates.
		    damage.position.y = std::stod(trim(values[1]));
		    damage.position.z = std::stod(trim(values[2]));
		}
		else								// If Field 2 does not contain either exactly 3 or 9 values, exit with error
		{
		    std::cerr << "ERROR: Incorrect number of entries for Data field 2.\n";
		    return false;
		}
            }

	    else								// If number of subfields is not exactly 1 or 3, exit with error.
	    {
		std::cerr << "ERROR: Incorrect number of entries for Data field 2.\n";
		return false;
	    }
	    
    	    fieldIndex++;							// Increment field index to compare with 'Data entries' header field.
        }

	//--------------------------------------------------
        // Field 3: Chromosome IDs
        //--------------------------------------------------

        if(hasField3)								// Is data Field 3 present?
        {

 	    if(fieldIndex >= fields.size()) 					// Make sure the Field Index matches the previous number of seen data fields.
            {
		std::cerr << "ERROR: Missing SDD data entry Field 3.\n";
        	return false;
    	    }

            std::vector<std::string> values = split(fields[fieldIndex], ',');	// Store field values

	    if(!trim(values[0]).empty() && values.size() == 4)			// This field should contain exactly four values
	    {
		damage.chromosomeID.dnaStructure = std::stoi(trim(values[0]));	// First value describes the DNA structure
	        damage.chromosomeID.chromosomeNumber = std::stoi(trim(values[1])); // Second value describes the chromosome number
	        damage.chromosomeID.chromatidNumber = std::stoi(trim(values[2]));  // Third value describes which chromatid is indicated (either 1 or 2)
	        damage.chromosomeID.chromosomeArm = std::stoi(trim(values[3]));	// Fourth value indicates which arm of the chromatid the damage is on (0 for short, 1 for long).
	    }
            else								// If number of entries in the field not equal to 4, exit with error.
	    {
		std::cerr << "ERROR: Incorrect number of entries in data Field 3.\n";
		return false;
	    }

	    fieldIndex++;							// Increment field index for comparison with 'Data entries' header field.
	}

	//--------------------------------------------------
        // Field 4: Chromosome Position
        //--------------------------------------------------

        if(hasField4)								// Check if data field 4 is present
        {

 	    if(fieldIndex >= fields.size()) 					// Make sure the field index matches the previous number of seen data fiels.
            {
		std::cerr << "ERROR: Missing SDD data entry Field 4.\n";
        	return false;
    	    }

            double value = std::stod(trim(fields[fieldIndex]));			// Field is a single value storing the location of the DNA damage

            damage.chromosomePosition.position = value;				// The value is either in base pairs or is fractions (if a number < 1 is specified)

            if(value < 1.0)
                damage.chromosomePosition.isFractional = true;
            else
                damage.chromosomePosition.isFractional = false;

	    fieldIndex++;							// Increment field index for comparison with 'Data entries' header field.
        }



 	//--------------------------------------------------
        // Field 5: Cause
        //--------------------------------------------------

        if(hasField5)								// check if field 5 present
        {

 	    if(fieldIndex >= fields.size()) 					// Make sure the field index matches the previous number of seen data fields.
            {
		std::cerr << "ERROR: Missing SDD data entry Field 5.\n";
        	return false;
    	    }
            
	    std::vector<std::string> values = split(fields[fieldIndex], ',');


	    if(values.size() != 3)						// Field 5 must contain exactly three values.
    	    {
                std::cerr << "ERROR: Field 5 must contain exactly 3 values.\n";
                return false;
    	    }

    	    damage.damageCause.cause = std::stoi(trim(values[0]));		// First entry stores the cause of damage from direct/indirect effects or charge migration

    	    damage.damageCause.numDirectDamages = std::stoi(trim(values[1]));	// Counter for the number of direct damages in the given damage block

    	    damage.damageCause.numIndirectDamages = std::stoi(trim(values[2])); // Counter for the number of indirect damages in the given damage block

    	    fieldIndex++;							// Increment field index for comparison with 'Data entries' header field.

	}


        //--------------------------------------------------
        // Field 6: Damage Types
        //--------------------------------------------------

        if(hasField6)								// Is field 6 present?
        {

 	    if(fieldIndex >= fields.size()) 					// Make sure the field index matches the previous number of seen data fields
            {
		std::cerr << "ERROR: Missing SDD data entry Field 6.\n";
        	return false;
    	    }

            std::vector<std::string> values = split(fields[fieldIndex], ',');

	    if(!values.empty() && values.size() == 3)				// If field is empty/does not contain exactly three values, exit with error
    	    {
		damage.damageType.numBaseDamages = std::stoi(trim(values[0]));	// First entry is the number of base damages in the site
		damage.damageType.numSingleBackboneBreaks = std::stoi(trim(values[1])); // Second entry is the number of single strand breaks in the site
                damage.damageType.presenceOfDoubleStrandBreaks = std::stoi(trim(values[2])); // Boolian flag to indicate if the single strand breaks consitute the presence of a DSB or not
	    }

	    else
	    {
		std::cerr << "ERROR: Incorrect number of entries in Field 6 - requires three comma separated values.\n";
	    	return false;
	    }

	    fieldIndex++;							// Increment field index for comparison with 'Data entries' header field.
	}


        //--------------------------------------------------
        // Field 7: Full Break Specification
        //--------------------------------------------------

        if(hasField7)								// Is field 7 present?
        {
 	    if(fieldIndex >= fields.size()) 					// Make sure the field index matches the previous number of seen data fields
            {
		std::cerr << "ERROR: Missing SDD data entry Field 7.\n";
        	return false;
    	    }

            std::vector<std::string> groups = split(fields[fieldIndex], '/');	// Field can have subgroups

            for(const auto& group : groups)					// Make sure all subgroups are non-empty
            {
                if(trim(group).empty())
                {
                    continue;
                }

		std::vector<std::string> values = split(group, ',');		// Each group consists of: STRAND, BASE, baseDamageType

		if(values.size() != 3)						// Each group must have three values
        	{
            	    std::cerr << "ERROR: Invalid Field 7 group: [" << group << "]\n";
            	    return false;
        	}

	        int strand = std::stoi(trim(values[0]));			// First value is the strand (numbered 1-4)

        	int base = std::stoi(trim(values[1]));				// Second strand is the base number from the break site specified in Field 4 from the p arm towards the q arm.

        	int baseDamageType = std::stoi(trim(values[2]));		// Damage type can be non-damaging (0), direct effects (1), indirect effects (2), or combined direct + indirect (3)

		// Validate STRAND
        	if(strand < 1 || strand > 4)					// Ensure the first value in the subgroup is always a number between 1 and 4
        	{
            	    std::cerr << "ERROR: Invalid Field 7 strand value: " << strand << "\n";

            	    return false;
        	}

        	// Validate base damage type
        	if(baseDamageType < 0 || baseDamageType > 3)			// Ensure damage type value is a valid option.
        	{
            	    std::cerr << "ERROR: Invalid Field 7 damage type: " << baseDamageType << "\n";

            	    return false;
        	}


		damage.fullBreakSpec.strand.push_back(strand);			// If strand values are valid, store them in the fullBreakSpec object.
        	damage.fullBreakSpec.base.push_back(base);			// If base values are valid, store them in the fullBreakSpec object.
        	damage.fullBreakSpec.baseDamageType.push_back(baseDamageType);  // If damage type values are valid, store them in the fullBreakSpec object

	   }
        
	    fieldIndex++;							// Make sure the field index matches the previous number of seen data fields.
	}


        //--------------------------------------------------
        // Add damage to the appropriate exposure
        //--------------------------------------------------

	if(currentExposureID <= 0 || currentExposureID > static_cast<int>(exposures.size()))
	{
    	    std::cerr << "ERROR: Invalid current exposure ID.\n";
   	    return false;
	}

	exposures[currentExposureID - 1].damages.push_back(damage);		// Exposure IDs start at 1, but indexing starts at 0, store damages to the exposure object.

    }

    return true;								// If no errors encountered previously, return true.
}













// -------------------------------------------------------------------------- //
// Functions to get header field / exposure fields and function to create 
// summaries of the header fields and exposure fields
// -------------------------------------------------------------------------- //


const Header& SDDparser::getHeader() const					// Getter function to return header fields
{
    return header;
}


const std::vector<Exposure>& SDDparser::getExposures() const 			// Getter function to return the exposure data fields
{
    return exposures;
}














// Function to print out the associated chromosome damages in a given exposure, mapping the exposure ID to the chromosome damages from that exposure.
std::map<int, ChromosomeDamageSummary> SDDparser::summarizeChromosomeDamage(
    const Exposure& exposure							// Access exposure data entries to summarize number of each damage type
    ) const
{
    std::map<int, ChromosomeDamageSummary> summary;				// initializing a summary object mapping the exposure id to the chromosome damages.

    for(const auto& damage : exposure.damages)
    {
        const auto& chromosome = damage.chromosomeID;				// For now we will summarize only the number of base damages, single strand breaks and double strand breaks in each chromosome for a given exposure.
        const auto& damageType = damage.damageType;

        if(chromosome.chromosomeNumber <= 0)					// Only strands that have an associated chromosome number (non-zero integer) will be valid.
        {
            continue;
        }

        int chromosomeNumber = chromosome.chromosomeNumber;

        ChromosomeDamageSummary& chromosomeSummary = summary[chromosomeNumber];	// Mapping the chromosome damage summary for each chromosome number to the given summary at that chromosome number.

        chromosomeSummary.dnaStructure = chromosome.dnaStructure;		

        chromosomeSummary.chromosomeNumber = chromosomeNumber;

        chromosomeSummary.numBaseDamages += damageType.numBaseDamages;		// Keep track of number of base damages for per exposure totals and totals over all exposures.

        chromosomeSummary.numSingleStrandBreaks += damageType.numSingleBackboneBreaks; // Keep track of number of single-strand breaks for per exposure totals and totals over all exposures.

        chromosomeSummary.numDoubleStrandBreaks += damageType.presenceOfDoubleStrandBreaks; // Keep track of boolian flag for presence of double strand breaks to be summed later for per exposure totals and totals over all exposures.
    }

    return summary;								// Summary returns number of base damages, single-strand breaks, and double strand breaks, per chromosome, per exposure, as well as per exposure totals and totals over all exposures.
}
















// Summary function to summarize the header fields only and output to a file, with suffix '_summary.txt' as specified in main.cpp. 
// Will be used in the larger printSummary function later with printExposureSummary.
void SDDparser::printHeaderSummary(
    std::ostream& output							// Takes as input the output path for the summary file
    ) const
{


	// ------------------- First section is the radiation information ---------------------- //
	output << "-----------------------     Incident Radiation Information"	
	       << "      -----------------------\n\n";

	output << "Radiation source: " << (header.source.empty() ? "N/A" : header.source) << "\n";		// Source summary - prints the string entry
	
	output << "Source type: ";
	if (header.source_type == SDD_INT_FIELD_NOT_MEASURED)
	{
	    output << "N/A";
	}
	else
	{
	   output << header.source_type;
	}
	output << " (" << sourceTypeMeaning(header.source_type) << ")\n";

	if (header.incident_particles.empty())					// Incident Particles summary
	{
	    output << "Incident Particles: N/A\n";
	}
	else
	{
	    output << "Incident particle fluence fractions: \n";		// Output fluence fraction associated with each listed incident particle
	    
	    if (header.particle_fraction.size() != header.incident_particles.size())
	    {
		output << "WARNING: 'Incident particles' (" << header.incident_particles.size()
	               << ") and 'Particle fraction' (" << header.particle_fraction.size()
	               << ") have mismatched lengths; showing only matched entries.\n";
	    }

	    const std::size_t matchedCount = std::min(header.incident_particles.size(), header.particle_fraction.size());

	    for (std::size_t i = 0; i < matchedCount; i++)
	    {
	        output << incidentParticlesMeaning(header.incident_particles[i]) << ": " << header.particle_fraction[i] << "\n";
	    }
	}

	output << "Mean particle energy: ";				// Mean particle energy summary
	if (header.mean_particle_energy.empty())
	{
	    output << "N/A";
	}
	else
	{
	    for (std::size_t i = 0; i < header.mean_particle_energy.size(); i++)
	    {
		output << header.mean_particle_energy[i];
		if (i + 1 < header.mean_particle_energy.size())
		{
		    output << ", ";
		}
	    }
	}
	output << "\n";

	output << "Dose or fluence specified: " << doseOrFluenceMeaning(header.dose_or_fluence) << "\n"; 	// Dose or fluence summary

	// ------------------ Second section is the target information ---------------------- //
	output << "-----------------------      Radiation Target Information" <<
		  "      -----------------------\n\n";

	output << "Radiation incident on: " << (header.irradiation_target.empty() ? "N/A" : header.irradiation_target) << "\n";// Irradiation target

	output << volumesMeaning(header.volumes) << "\n" ;			// Summary of target shapes and volumes

	output << "The number of chromosomes specified were: " << header.chromosome_sizes[0] << // Summarize number of chromosomes and DNA density
	". Subsequent field entries are the chromosome sizes in units of Mega base pairs (Mbp),with an average DNA density of ";
 	if (header.DNA_density == SDD_DOUBLE_FIELD_NOT_MEASURED)
	{
	    output << "N/A";
	}
	else
	{
	    output << header.DNA_density;
	}
	output << " Mbp per cubic micrometer.\n";

	output << "Cell Cycle Phase: " << cellCyclePhaseMeaning(header.cell_cycle_phase) << "\n" ;	// Summarize cell cycle phase

	output << "The DNA is structured as: ";
	if (header.DNA_structure.empty())
	{
	    output << "N/A";
	}
	else
	{
	   output << header.DNA_structure[0];
	}
	output << " (" << dnaStructureMeaning(header.DNA_structure) << ")\n";

	output << "The cell Proliferation status is: ";
	if (header.proliferation_status.empty())
	{
	    // proliferationStatusMeaning() indexes status[0] without checking
	    // it's non-empty, so it must not be called here at all when blank.
	    output << "N/A\n";
	}
	else
	{
	    output << header.proliferation_status[0] << " (" << proliferationStatusMeaning(header.proliferation_status) << ")" << "\n";
	}

	// ------------------- Third section is the DNA damage and exposure information ---------------------- //
	output << "\n-----------------------      DNA Damage Information" <<
		  "      -----------------------\n\n";

	output << damageDefinitionMeaning(header.damage_definition) << "\n"; // Summary of the Damage definition and damage and primary count header fields
	output << "The number of distinct damage lesions scored is ";
	if (header.damage_and_primary_count.size() < 2)
	{
	    output << "N/A";
	}
	else
	{
	    output << header.damage_and_primary_count[0] << ", as a result of " << header.damage_and_primary_count[1] << " primary particles simulated";
	}
	output << ".\n";


	output << microenvironmentMeaning(header.microenvironment) << "\n"; 	// Summary of microenvironment (temperature and oxygenation).

	if (header.time == SDD_DOUBLE_FIELD_NOT_MEASURED)			// timeMeaning() doesn't know about the sentinel and would otherwise
	{									// report it as an invalid entry, so bypass it entirely when absent.
	    output << "Time: N/A\n";
	}
	else
	{
	    output << timeMeaning(header.time) << "\n";
	}
	output << dataEntriesMeaning(header.data_entries) << "\n";		// Summary of Data field entries
}















// Summary function to summarize the data fields only and output to a file, with suffix '_summary.txt' as specified in main.cpp. 
// Will be used in the larger printSummary function later with printHeaderSummary.
void SDDparser::printExposureSummary(
    std::ostream& output							// Takes as input output file path for summary file
    ) const
{
	output << "----------------------------------------------------------------------------\n";	
	output << "-----------------------      Chromosome Damages      -----------------------\n";
	output << "----------------------------------------------------------------------------\n";

	output << "\nNumber of exposures: " << exposures.size() << "\n";	// Report total number of exposures in the SDD file.

	int overallBaseDamages = 0;						// Counter for the number of base damages for all exposures and all chromosomes
	int overallSingleStrandBreaks = 0;					// Counter for the number of single strand breaks for all exposures and all chromosomes
	int overallDoubleStrandBreaks = 0;					// Counter for the number of double strand breaks for all exposures and all chromosomes

    for(const auto& exposure : exposures)					// Counting the damages in each given exposure
    {
        output << "\n-------------------------------------------\n";
        output << "Exposure " << exposure.exposureID << "\n";
        output << "-------------------------------------------\n";

        output << "Number of damage entries: "					// Report the number of damage entries (rows) for a given exposure.
            << exposure.damages.size()
            << "\n";

        auto chromosomeSummary =
            summarizeChromosomeDamage(exposure);				// Call the above summarizeChromosomeDamage function for each exposure.

        if(chromosomeSummary.empty())						// Make sure there are damages to report before reporting.
        {
            output << "\nNo chromosome-associated damage entries found.\n";
            continue;
        }

        output << "\nDamage by chromosome:\n\n";				// Report the number of base damages, single strand breaks, and double strand breaks per chromosome.

	int totalBaseDamages = 0;						// per exposure totals
	int totalSingleStrandBreaks = 0;
	int totalDoubleStrandBreaks = 0;

        for(const auto& [chromosomeNumber, summary] :				// Sum the total number of different damages for each chromosome in a given exposure before reporting
            chromosomeSummary)
        {
	    totalBaseDamages += summary.numBaseDamages;
	    totalSingleStrandBreaks += summary.numSingleStrandBreaks;
	    totalDoubleStrandBreaks += summary.numDoubleStrandBreaks;

            output << "Chromosome " << chromosomeNumber << "\n";

            output << "  Base damages: " << summary.numBaseDamages << "\n";

            output << "  Single-strand breaks: " << summary.numSingleStrandBreaks << "\n";

            output << "  Double-strand breaks: " << summary.numDoubleStrandBreaks << "\n";

	}
										// Report the total number of base damages, single strand breaks and double strand breaks in a given exposure.
	output << "\nExposure totals:\n";

	output << "  Base damages: " << totalBaseDamages << "\n";

	output << "  Single-strand breaks: " << totalSingleStrandBreaks << "\n";

	output << "  Double-strand breaks: " << totalDoubleStrandBreaks << "\n";


        overallBaseDamages += totalBaseDamages;					// Sum exposure totals and add to the over all exposures totals
        overallSingleStrandBreaks += totalSingleStrandBreaks;
        overallDoubleStrandBreaks += totalDoubleStrandBreaks;

    }


    output << "\n===============================================\n";
    
    output << "Overall chromosome-associated damage:\n";			// Report the overall damages for all chromosomes over all exposures at the bottom of the summary file.

    output << "  Base damages: "
           << overallBaseDamages << "\n";

    output << "  Single-strand breaks: "
           << overallSingleStrandBreaks << "\n";

    output << "  Double-strand breaks: "
           << overallDoubleStrandBreaks << "\n";

    output << "=================================================\n";

}














// Function that outputs to a summary file the header and data field summary. User can have the option of commenting out one or the other if they
// do not want either the header or the data field summary.
void SDDparser::printSummary(
    std::ostream& output							// Takes as input output file summary path
    ) const
{
    output << "==========================================================================\n";
    output << "=======================      SDD FILE SUMMARY	   =======================\n";
    output << "==========================================================================\n";
    output << "Generated by " << SDDPARSER_VERSION_TEXT << "\n\n";

    printHeaderSummary(output);							// Outputs header summary to file

    output << "\n";

    printExposureSummary(output);						// Outputs data field summary to a file


}

} // namespace sddparser
