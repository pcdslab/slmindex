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

#include "dslim.h"

/* Sanity check for fragment ion mass tolerance */
#if (dF <= 0)
#error "ERROR: The fragment mass tolerance must be > 0"
#endif /* (dF <= 0) */

/* External Variables */
extern SLMindex          dslim; /* DSLIM Index       */

#ifdef FUTURE
extern pepEntry    *pepEntries; /* SLM Peptide Index */
extern PepSeqs          seqPep; /* Peptide sequences */
#ifdef VMODS
extern varEntry    *modEntries;
#endif /* VMODS */
#endif /* FUTURE */

/* Global Variables */
extern UINT chunksize;

/* FUNCTION: DSLIM_QuerySpectrum
 *
 * DESCRIPTION: Query the DSLIM for all query peaks
 *              and count the number of hits per chunk
 *
 * INPUT:
 * @QA     : Query Spectra Array
 * @len    : Number of spectra in the array
 * @Matches: Array to fill in the number of hits per chunk
 * @threads: Number of parallel threads to launch
 *
 * OUTPUT:
 * @status: Status of execution
 */
STATUS DSLIM_QuerySpectrum(UINT *QA, UINT len, UINT *Matches, UINT threads)
{
    STATUS status = SLM_SUCCESS;
    UINT chunk_number = 0;
    UINT *QAPtr = NULL;

#ifndef _OPENMP
    threads = 1;
#endif /* _OPENMP */

        /* Query each chunk in parallel */
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif /* _OPENMP */
        /* Process all the queries in the chunk */
        for (UINT queries = 0; queries < len; queries++)
        {
            /* Pointer to each query spectrum */
            QAPtr = QA + (queries * QALEN);

            UCHAR *SCPtr = dslim.pepChunks[omp_get_thread_num()].sC;

            /* Check if this chunk is the last chunk */
            UINT size = chunksize;

			UINT localMatches = 0;
							
            /* Query each chunk in parallel */
            for (chunk_number = 0; chunk_number < dslim.nChunks; chunk_number++)
            {
                UINT *bAPtr = dslim.pepChunks[chunk_number].bA;
                UINT *iAPtr = dslim.pepChunks[chunk_number].iA;

#ifdef FUTURE

                UCHAR *bits = dslim.pepChunks[chunk_number].bits;
                /* Calculate the size of bits */
                UINT bitsize = (size/BYTE) + 1;

#endif /* FUTURE */

                /* Query all fragments in each spectrum */
                for (UINT k = 0; k < QALEN; k++)
                {
                    /* Check for any zeros
                     * Zero = Trivial query */
                    if (QAPtr[k] < dF || QAPtr[k] > ((MAX_MASS * SCALE) - 1 - dF))
                    {
                        continue;
                    }

                    /* Locate iAPtr start and end */
                    UINT start = bAPtr[QAPtr[k] - dF];
                    UINT end = bAPtr[QAPtr[k] + 1 + dF];

                    /* Loop through located iAions */
                    for (UINT ion = start; ion < end; ion++)
                    {
                        /* Calculate parent peptide ID */
                        UINT ppid = (iAPtr[ion] / (2 * F));

#ifdef FUTURE
                        if (IS_INIT(bits[ppid], ppid))
                        {
                            SC[ppid]++;
                        }
                        else
                        {
                            INIT(bits[ppid], ppid);
                            SC[ppid] = 1;
                        }
#else
                        /* Update corresponding SC entry */
                        SCPtr[ppid] += 1;
#endif /* FUTURE */
                    }
                }

                /* Put the candidate peptides
                 * from each chunk
                 */
                for (UINT cntr = 0; cntr < size && localMatches < MAXMATCHES; cntr++)
                {
                    if (SCPtr[cntr] >= MIN_SHRD_PKS)
                    {
                        /* Avoid too many updates to the Matches */
                        Matches[(queries * MAXMATCHES)+ localMatches] = cntr;

                        /* Fill in the matches here or just count them */
                        localMatches++;

                    }
                }
            }
#ifdef FUTURE
                /* Clear the bitmask */
                std::memset(bits, 0x0, bitsize);
#else
                /* bitmask not active,
                 * reset the SC instead */
                std::memset(SCPtr, 0x0, size);
#endif /* FUTURE */
        }

    return status;
}
