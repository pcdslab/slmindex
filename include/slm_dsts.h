/*
 * This file is part of SLM-Transform
 *  Copyright (C) 2019  Muhammad Haseeb, Fahad Saeed
 *  Florida International University, Miami, FL
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#ifndef INCLUDE_SLM_DSTS_H_
#define INCLUDE_SLM_DSTS_H_

#include "common.h"

/* Types of modifications allowed by SLM_Mods     */
#define MAX_MOD_TYPES                        7

/************************* Common DSTs ************************/
/* Add distribution policies */
typedef enum _DistPolicy
{
    _chunk,
} DistPolicy;

typedef struct _SLM_varAA
{
    AA     residues[4]   = ""; /* Modified AA residues in this modification - Upto 4 */
    UINT  modMass         = 0; /* Scaled mass of the modification                    */
    USHORT aa_per_peptide = 0; /* Allowed modified residues per peptide sequence     */

    _SLM_varAA& operator=(const _SLM_varAA& rhs)
    {
        /* Check for self assignment */
        if (this != &rhs)
        {
            this->aa_per_peptide = rhs.aa_per_peptide;
            this->modMass = rhs.modMass;

            this->residues[0] = rhs.residues[0];
            this->residues[1] = rhs.residues[1];
            this->residues[2] = rhs.residues[2];
            this->residues[3] = rhs.residues[3];
        }

        return *this;
    }

} SLM_varAA;

typedef struct _SLM_Mods
{
    USHORT       vmods_per_pep = 0; /* Total allowed modified residues per sequence */
    USHORT            num_vars = 0; /* Number of types of modifications added to index - Max: 7 */
    SLM_varAA vmods[MAX_MOD_TYPES]; /* Information for each modification added to index */

    /* Overload = operator */
    _SLM_Mods& operator=(const _SLM_Mods& rhs)
    {
        /* Check for self assignment */
        if (this != &rhs)
        {
            this->vmods_per_pep = rhs.vmods_per_pep;
            this->num_vars = rhs.num_vars;

            for (UINT i = 0; i < MIN(num_vars, 7); i++)
            {
                this->vmods[i] = rhs.vmods[i];
            }
        }

        return *this;
    }

} SLM_vMods;

typedef struct _pepSeq
{
    AA    *seqs; /* Stores peptide sequence, could store as strings as well */
    UINT   *idx; /* Stores sequence length */
    UINT    AAs; /* Total number of characters */
} PepSeqs;

typedef struct _modAA
{
    ULONG  sites = 0x0; /* maxlen(pep) = 60AA + 2 bits (termini mods)      */
    UINT  modNum = 0x0; /* 4 bits per mods num, Max 8 mods allowed per pep */

    /* Overload = operator */
    _modAA& operator=(const _modAA& rhs)
    {
        /* Check for self assignment */
        if (this != &rhs)
        {
            this->sites = rhs.sites;
            this->modNum = rhs.modNum;
        }

        return *this;
    }

} modAA;

typedef struct _pepEntry
{
    FLOAT  Mass; /* mass of peptide */
} pepEntry;

typedef struct _varEntry
{
    FLOAT  Mass; /* Mass of Peptide            */
    IDX   seqID; /* Normal Peptide Sequence ID */
    modAA sites; /* Modified AA information    */

    /* Overload = operator */
    _varEntry& operator=(const _varEntry& rhs)
    {
        /* Check for self assignment */
        if (this != &rhs)
        {
            this->Mass = rhs.Mass;
            this->seqID = rhs.seqID;
            this->sites.modNum = rhs.sites.modNum;
            this->sites.sites = rhs.sites.sites;
        }
        return *this;
    }

    /* Default constructor */
    _varEntry()
    {
        Mass = 0;
        seqID = 0;
        sites.modNum = 0;
        sites.sites = 0;
    }

} varEntry;

/************************* SLM Index DSTs ************************/
typedef struct _SLMchunk
{
    UINT    *iA = NULL; /* Ions Array (iA)   */
    UINT    *bA = NULL; /* Bucket Array (bA) */
    UCHAR   *sC = NULL; /* Scorecard (SC)    */
#ifdef FUTURE
    UCHAR *bits = NULL; /* Scorecard bits    */
#endif /* FUTURE */
} SLMchunk;

typedef struct _SLMindex
{
    UINT       nPepChunks =    0; /* Number of pep chunks         */
    UINT       nModChunks =    0; /* Number of mod chunks         */
    UINT          nChunks =    0; /* Total Number of chunks       */
    SLMchunk   *pepChunks = NULL; /* Chunks for normal peptides   */
    SLMchunk   *modChunks = NULL; /* Chunks for modified peptides */
} SLMindex;


/************************* SLM Query DSTs ************************/
typedef struct _Query
{
    /* Raw chosen fragments/peaks from the input MS/MS spectrum   */
    PEAK Peaks[QALEN] = {0};

/* HM: Enable if intensity information is also required in future */
#if (defined(REQUIRE_INTENSITY))
    INTENSITY Intensities[QUERYPK] = {0};
#endif /* (defined(REQUIRE_INTENSITY)) */
} Query;

#endif /* INCLUDE_SLM_DSTS_H_ */
