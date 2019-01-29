# SLMTransform
SLM-Transform: A Method for memory-efficient indexing of spectra for database search in LC-MS/MS proteomics

## Authors
Muhammad Haseeb and Fahad Saeed

# What do you need?
1. make
2. GCC v5.4.0 or later with C++11+ support
2. Windows: MinGW (x86_64-6.3.0-win32-seh-rt_v5-rev2)

# Pre SLM-Transform Steps
1. Digest the proteome database using Protein Digestion Simulator or OpenMS.
2. Remove the redundant peptide sequences in the digested database using DBToolkit.
3. Optional: The decoy database can be generated and appended to the target database using DBToolkit.
4. Convert the MS/MS data in (mzML/mzXML/MS2) format (MS2 preferably) using msconvert.exe

# Configure SLM-Transform & Sample Application
1. Configure the SLM-Transform parameters in /slmtransform/include/config.h: 
*Make sure that #define/undef WINDOWS is correct according to your Host OS else there will be Seg faults.*
2. Add the database, dataset and mods information in the ./slm.cpp (Sample Application) using the following format:

## Format

    STRING modconditions = "#1 <res> #2"; 
    /* #1 is maximum allowed modified residues per sequence, 
    <res> modified residues of 1st type, 
    #2 how many residues of this modification in peptide sequence. 
    
    For example: we want max 5 mod residues; and max STY 3 and max M 2 we will set modconditions as
    "5 STY 3 M 2"
    */
    
    <some code> 
    
    /* Database and Dataset files */
    STRING filename = "/path/to/digested/and/unduplicated/database.fasta"; // Path to processed database
    STRING querypath = "/path/to/query/dataset"; // Path to Dataset folder
    STRING patt = ".ms2/mzML/mzXML"; // Keep one, remove others
    
    <some code> 
    
    /* Initialize the vModInfo */
    vModInfo.num_vars = #3; // Types of modifications been specified (max 7) for above example, set this to: 2
    vModInfo.vmods_per_pep = #1; // #1 from the previous explanation, from above example, set this to: 5
    
    /* List of Mods Info */
    vModInfo.vmods[0].aa_per_peptide = 3// #2 from previous explanation, for above example, set this to: 3
    vModInfo.vmods[0].modMass = 79.97 * SCALE; // Mass of the specified modification, for above example: 79.97
    vModInfo.vmods[0].residues[0] = 'S'; // modified residues list (max 4 per mod type allowed), for above example S
    vModInfo.vmods[0].residues[1] = 'T'; // T
    vModInfo.vmods[0].residues[2] = 'Y'; // Y
    
    vModInfo.vmods[0].aa_per_peptide = 2; // #2 from previous explanation, for above example, set this to: 3
    vModInfo.vmods[0].modMass = 15.997 * SCALE; // Mass of the specified modification, for above example: 79.97
    vModInfo.vmods[0].residues[0] = 'M'; // modified residues list (max 4 per mod type allowed), for above example S

# Building SLM-Transform Sample Application
1. Open the Git Bash shell (Windows) or normal Terminal (Ubuntu).
1. Navigate to SLM-Transform home directory /slmtransform/
2. Execute the following command: `make`

# Running SLM-Transform Sample Application
1. Navigate to /slmtransform
2. ./SLMTransform.exe

# Please Note:
1. Work is being done to move SLM-Transform configuration to runtime so please bear with us.
2. Max digested peptide mass allowed: 5000Da.
3. The SLM-Transform only returns the fragment-ion filtered PSMs, which can be used in any way e.g. post filtered based on precursor masses or sequence tags, formally scored for PSMs, FDRed.
4. You are free to fork this repo and make any changes for your application (under the specified license), however, please conform to code format in: lbe_eclipse.xml

# Please cite our work
For queries or questions about SLM-Transform, please contact: {fsaeed, mhaseeb}@fiu.edu. Thank you.
