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

/* Sanity check for CHUNKSIZE */
#if (CHUNKSIZE <= 0)
#error "ABORT: The macro CHUNKSIZE must be > 0"
#endif /* (CHUNKSIZE == 0) */

/* Enternal Variables */
extern ULONGLONG    pepCount;
extern ULONGLONG  totalCount;
extern PepSeqs        seqPep;

/* Global Variables */
SLMindex          dslim; /* DSLIM Index       */
pepEntry    *pepEntries; /* SLM Peptide Index */
UINT           *SpecArr; /* Spectra Array     */

#ifdef VMODS
varEntry    *modEntries; /* SLM Mods Index    */
#endif /* VMODS */

/* Global Variables */
UINT chunksize = CHUNKSIZE;
UINT lastchunksize = 0;
UINT nchunks   = 0;


/* FUNCTION: DSLIM_Construct
 *
 * DESCRIPTION: Construct DSLIM chunks
 *
 * INPUT:
 * @threads: Number of parallel threads to launch
 * @modInfo: DSLIM mods Information
 *
 * OUTPUT:
 * @status: status of execution
 */
STATUS DSLIM_Construct(UINT threads, SLM_vMods *modInfo)
{
    STATUS status = SLM_SUCCESS;

#ifdef VMODS
    /* Update gModInfo */
    status = UTILS_InitializeModInfo(modInfo);
#else
    LBE_UNUSED_PARAM(modInfo);
#endif /* VMODS */

#ifndef _OPENMP
    threads = 1;
#endif /* _OPENMP */

    if (modInfo == NULL)
    {
        status = ERR_INVLD_PTR;
    }

    if (status == SLM_SUCCESS)
    {
        /* Spectra Array (SA) */
        SpecArr = new UINT[chunksize * iSERIES * F];

        /* Check if Spectra Array has been allocated */
        if (SpecArr == NULL)
        {
            status = ERR_BAD_MEM_ALLOC;
        }
    }
    if (status == SLM_SUCCESS)
    {
        /* Allocate memory for SLMChunks and SPI*/
        status = DSLIM_AllocateMemory(chunksize, nchunks);

        /* Construct DSLIM.iA */
        if (status == SLM_SUCCESS)
        {
            /* Distributed SLM Ions Array construction */
            for (UINT chno = 0; chno < nchunks && status == SLM_SUCCESS; chno++)
            {
                /* Construct each DSLIM chunk in Parallel */
                status = DSLIM_ConstructChunk(threads, chno);

                /* Apply SLM-Transform on the chunk */
                if (status == SLM_SUCCESS)
                {
                    status = DSLIM_SLMTransform(threads, chno);
                }
            }
        }
    }

    if (status == SLM_SUCCESS)
    {
        /* Construct DSLIM.bA */
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static, 1)
#endif /* _OPENMP */
        for (UINT chunk_number = 0; chunk_number < nchunks; chunk_number++)
        {
            UINT *bAPtr = dslim.pepChunks[chunk_number].bA;

            /* Get size of the chunk */
            UINT csize = ((chunk_number == nchunks - 1) && (nchunks > 1)) ?
                            lastchunksize                                 :
                            chunksize;

            UINT count = bAPtr[0];

            /* Initialize first and last bA entries */
            bAPtr[0] = 0;
            bAPtr[(MAX_MASS * SCALE)] = (csize * iSERIES * F);

            for (UINT li = 1; li <= (MAX_MASS * SCALE); li++)
            {
                UINT tmpcount = bAPtr[li];
                bAPtr[li] = bAPtr[li - 1] + count;
                count = tmpcount;
#ifdef VALIDATE_SLM
                if (bAPtr[li] < bAPtr[li - 1])
                {
                    std::cout << chunk_number << " " << li << " " << bAPtr[li - 1] << " " << count << " " << bAPtr[li] << std::endl;
                    status = ERR_INVLD_SIZE;

                    while(true);
                }
#endif /* VALIDATE_SLM */
            }

            /* Check if all correctly done */
            if (bAPtr[(MAX_MASS * SCALE)] != (csize * iSERIES * F))
            {
                status = ERR_INVLD_SIZE;
            }
        }
    }

    /* Remove the temporary SpecArray (SA) */
    delete[] SpecArr;
    SpecArr = NULL;

    return status;
}

/*
 * FUNCTION: DSLIM_AllocateMemory
 *
 * DESCRIPTION: Allocate memory for DSLIM chunks
 *
 * INPUT:
 * @chsize: Chunk size
 * @Chunks: Number of chunks
 *
 * OUTPUT:
 * @status: Status of execution
 */
STATUS DSLIM_AllocateMemory(UINT chsize, UINT Chunks)
{
    STATUS status = SLM_SUCCESS;

    /* Initialize DSLIM pepChunks */
    dslim.pepChunks = new SLMchunk[Chunks];

    if (dslim.pepChunks != NULL)
    {
        dslim.nChunks = Chunks;

        /* Counter for Chunk size construction */
        INT totalpeps = (INT) totalCount;

        /* At every loop, check for remaining peps and status */
        for (UINT i = 0; i < Chunks && totalpeps > 0 && status == SLM_SUCCESS; i++)
        {
            /* Initialize direct hashing bA */
            dslim.pepChunks[i].bA = new UINT[(MAX_MASS * SCALE) + 1];

            if (dslim.pepChunks[i].bA != NULL)
            {
                /* Calculate the iA chunk size */
                INT size = ((INT)(totalpeps - chsize)) > 0 ? chsize : totalpeps;
                totalpeps -= size; // Update the counter

                /* Total Number of Ions = peps * #ion series * ions/ion series */
                dslim.pepChunks[i].iA = new UINT[(size * iSERIES * F)];

                if (dslim.pepChunks[i].iA == NULL)
                {
                    status = ERR_INVLD_MEMORY;
                }
            }
            else
            {
                status = ERR_INVLD_MEMORY;
            }
        }
    }
    else
    {
        status = ERR_INVLD_MEMORY;
    }

    /* Allocate memory for SPI and Spectra Array */
    if (status == SLM_SUCCESS)
    {
        pepEntries = new pepEntry[pepCount];

        if (pepEntries == NULL)
        {
            status = ERR_INVLD_MEMORY;
        }
    }

    return status;
}

/*
 * FUNCTION: DSLIM_ConstructChunk
 *
 * INPUT:
 * @threads:      Number of parallel threads
 * @chunk_number: Chunk Index
 *
 * OUTPUT:
 * @status: Status of execution
 */
STATUS DSLIM_ConstructChunk(UINT threads, UINT chunk_number)
{
    STATUS status = SLM_SUCCESS;

    /* Check if this chunk is the last chunk */
    BOOL lastChunk = (chunk_number == (nchunks - 1))? true: false;

#ifdef _OPENMP
    /* SA ptr for each thread */
    UINT *SAPtrs[threads] = {NULL};

    /* Temporary Array needed for Theoretical Spectra */
    UINT *Spectra = new UINT[threads * MAX_SEQ_LEN * iSERIES * MAXz];

    UINT *BAPtrs[threads] = {NULL};
    UINT *bA = new UINT[threads * SCALE * MAX_MASS];

    /* Initialize SAPtrs for each thread */
    if (Spectra != NULL && bA != NULL)
    {
        for (UINT i = 0; i < threads; i++)
        {
            SAPtrs[i] = Spectra + (i * MAX_SEQ_LEN * iSERIES * MAXz);
            BAPtrs[i] = bA + (i * SCALE * MAX_MASS);
        }

        /* Clear the bAPtrs and SAPtrs */
#pragma omp parallel for num_threads(threads) schedule (static)
        for (UINT i = 0; i < threads; i++)
        {
            std::memset(BAPtrs[i], 0x0, (SCALE * MAX_MASS * sizeof(UINT)));
            std::memset(SAPtrs[i], 0x0, (MAX_SEQ_LEN * iSERIES * MAXz * sizeof(UINT)));
        }
    }
    else
    {
        status = ERR_INVLD_MEMORY;
    }
#else

    UINT SAPtr[(MAX_SEQ_LEN * iSERIES * MAXz)] = {};
    std::memset(dslim.pepChunks[chunk_number].bA, 0x0, (SCALE * MAX_MASS * sizeof(UINT) + 1));

#endif /* _OPENMP */

    if (status == SLM_SUCCESS)
    {
        UINT start_idx = chunk_number * chunksize;
        UINT interval = chunksize;

        /* Check for last chunk */
        if (lastChunk == true && nchunks > 1)
        {
            interval = lastchunksize;
        }

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static, 1) /* schedule(dynamic, 1) */
#endif /* _OPENMP */
        for (UINT k = start_idx; k < (start_idx + interval); k++)
        {
            /* Filling point */
            UINT nfilled = (k - start_idx) * iSERIES * F;

#ifdef _OPENMP
            UINT *Spec = SAPtrs[omp_get_thread_num()];
            UINT *bAPtr  = BAPtrs[omp_get_thread_num()];
#else
            UINT *Spec = SAPtr;
            UINT *bAPtr = dslim.pepChunks[chunk_number].bA;
#endif /* OPENMP */

            /* Extract peptide Information */
            FLOAT pepMass = 0.0;
            CHAR *seq = NULL;
            INT len = 0;
            UINT pepID = LBE_RevDist(k);

#ifdef VMODS
            /* Check if pepID belongs to peps or mods */
            if (pepID >= pepCount)
            {
                /* Extract from Mods */
                varEntry *entry = modEntries + (pepID - pepCount);
                seq = &seqPep.seqs[seqPep.idx[entry->seqID]];
                len = (INT)seqPep.idx[entry->seqID + 1] - (INT)seqPep.idx[entry->seqID];

                /* Generate the Mod. Theoretical Spectrum */
                pepMass = UTILS_GenerateModSpectrum(seq, (UINT)len, Spec, entry->sites);

                /* Fill in the pepMass */
                entry->Mass = pepMass;
            }
            else
#endif /* VMODS */
            {
                /* Extract from Peps */
                pepEntry *entry = pepEntries + pepID;
                seq = &seqPep.seqs[seqPep.idx[pepID]];
                len = (INT)seqPep.idx[pepID + 1] - (INT)seqPep.idx[pepID];

                /* Generate the Theoretical Spectrum */
                pepMass = UTILS_GenerateSpectrum(seq, len, Spec);

                /* Fill in the pepMass */
                entry->Mass = pepMass;
            }

            /* If a legal peptide */
            if (pepMass >= MIN_MASS && pepMass <= MAX_MASS)
            {
                /* Calculate the length/2 of Spec */
                UINT half_len = ((iSERIES * MAX_SEQ_LEN * MAXz) / 2);

                /* Sort by ion Series and Mass */
                UTILS_Sort<UINT>(Spec, half_len, false);
                UTILS_Sort<UINT>((Spec + half_len), half_len, false);

                /* Choose SpecArr filling method */
                UINT nIons = (len - 1) * MAXz;
                UINT maxfill = F;
                UINT fill = 0;

                /* Fill the sorted ions */
                for (fill = 0; fill < nIons && fill < maxfill; fill++)
                {
                    /* Check if legal b-ion */
                    if (Spec[half_len - nIons + fill] >= (MAX_MASS * SCALE))
                    {
                        Spec[half_len - nIons + fill] = (MAX_MASS * SCALE) - 1;
                    }

                    SpecArr[nfilled + fill] = Spec[half_len - nIons + fill]; // Fill in the b-ion
                    bAPtr[SpecArr[nfilled + fill]]++; // Update the BA counter

                    /* Check for a legal y-ion */
                    if (Spec[(2 * half_len) - nIons + fill] >= (MAX_MASS * SCALE))
                    {
                        Spec[(2 * half_len) - nIons + fill] = (MAX_MASS * SCALE) - 1;
                    }

                    SpecArr[nfilled + F + fill] = Spec[(2 * half_len) - nIons + fill]; // Fill in the y-ion
                    bAPtr[SpecArr[nfilled + F + fill]]++; // Update the BA counter

                }

                /* Fill the rest with zeros */
                for (; fill < maxfill; fill++)
                {
                    SpecArr[nfilled + fill] = 0; // Fill in the b-ion
                    SpecArr[nfilled + F + fill] = 0; // Fill in the y-ion
                    bAPtr[0] += 2; // Update the BA counter
                }
            }

            /* Illegal peptide, fill in the container with zeros */
            else
            {
                /* Fill zeros for illegal peptides
                 * FIXME: Should not be filled into the chunk
                 *        and be removed from SPI as well
                 */
                std::memset(&SpecArr[nfilled], 0x0, sizeof(UINT) * iSERIES * F);
                bAPtr[0] += (iSERIES * F); // Update the BA counter
            }
        }

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule (static)
        for (UINT i = 0; i < (MAX_MASS * SCALE); i++)
        {
            /* Initialize to zero */
            dslim.pepChunks[chunk_number].bA[i] = 0;

            for (UINT j = 0; j < threads; j++)
            {
                dslim.pepChunks[chunk_number].bA[i] += BAPtrs[j][i];
            }
        }
#endif /* _OPENMP */

    }

#ifdef _OPENMP
    delete[] bA;
    bA = NULL;
#endif /* OPENMP */

    return status;
}

/*
 * FUNCTION: DSLIM_SLMTransform
 *
 * DESCRIPTION: Constructs SLIM Transform
 *
 * INPUT:
 * @threads     : Number of parallel threads
 * @chunk_number: Chunk Index
 *
 * OUTPUT:
 * @status: Status of execution
 */
STATUS DSLIM_SLMTransform(UINT threads, UINT chunk_number)
{
    STATUS status = SLM_SUCCESS;

    /* Check if this chunk is the last chunk */
    UINT size = ((chunk_number == nchunks - 1) && (nchunks > 1))?
                 lastchunksize                                  :
                 chunksize;

    UINT *iAPtr = dslim.pepChunks[chunk_number].iA;
    UINT iAsize = size * iSERIES * F;

#ifdef VALIDATE_SLM
    UINT *integ = new UINT[size * iSERIES * F];
    std::memcpy(integ, SpecArr, sizeof(UINT) * iAsize);
#endif /* VALIDATE_SLM */

    /* Construct DSLIM.iA */
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif /* _OPENMP */
    for (UINT k = 0; k < iAsize; k++)
    {
        iAPtr[k] = k;
    }

#ifdef _OPENMP
    /* Parallel Key Value Sort */
    KeyVal_Parallel<UINT>(SpecArr, iAPtr, iAsize, threads);
#else
    KeyVal_Serial<UINT>(SpecArr, iAPtr, iAsize);
#endif /* _OPENMP */

    /* Check integrity of SLM-Transform */
#ifdef VALIDATE_SLM
    BOOL integrity = true;

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif /* _OPENMP */
    for (UINT k = 1; k < iAsize; k++)
    {
        if (integ[iAPtr[k]] < integ[iAPtr[k - 1]])
        {
            integrity = false;
        }
    }

    if (integrity == false)
    {
        while (!integrity);
        status = ERR_INVLD_SORT;
    }
#endif /* VALIDATE_SLM */

    return status;
}

/*
 * FUNCTION: DSLIM_InitializeSC
 *
 * DESCRIPTION: Initialize Scorecard for DSLIM
 *
 * INPUT:
 * @threads: Number of parallel threads
 *
 * OUTPUT:
 * @status: Status of execution
 */
STATUS DSLIM_InitializeSC(UINT threads)
{
    STATUS status = SLM_SUCCESS;

    UINT chunk_number = 0;

#ifndef _OPENMP
    threads = 1;
#endif /* OPENMP */

    /* Loop through all the chunks */
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif /* _OPENMP */
    for (chunk_number = 0; chunk_number < nchunks; chunk_number++)
    {
        /* Check if this chunk is the last chunk */
        UINT size = chunksize;
#ifdef FUTURE
        /* Calculate the size of bits */
        UINT bitsize = (size/BYTE) + 1;
#endif /* FUTURE */

        /* Allocate memory for scorecard */
        UCHAR *SC = new UCHAR[size];
        dslim.pepChunks[chunk_number].sC = SC;

        if (SC != NULL)
        {
            std::memset(SC, 0x0, size);
#ifndef FUTURE
#else
            /* Allocate memory for bitmask */
            UCHAR *bits = new UCHAR[bitsize];

            if (bits != NULL)
            {
                /* Set the bitmask to zero */
                std::memset(bits, 0x0, bitsize);

                /* Set the pointers */
                dslim.pepChunks[chunk_number].bits = bits;
            }
            else
            {
                status = ERR_INVLD_MEMORY;
            }
#endif /* FUTURE */
        }
        else
        {
            status = ERR_INVLD_MEMORY;
        }
    }

    return status;
}

/*
 * FUNCTION: DSLIM_Analyze
 *
 * DESCRIPTION: Analyze the DSLIM distribution
 *
 * INPUT:
 * @threads: Number of parallel threads
 * @avg    : Pointer to DSLIM mean load
 * @std    : Pointer to DSLIM load distribution
 *
 * OUTPUT:
 * @status: Status of execution
 */
STATUS DSLIM_Analyze(UINT threads, DOUBLE &avg, DOUBLE &std)
{
    STATUS status = SLM_SUCCESS;
    DOUBLE sum_std = 0;
    ULONGLONG truecount = 1;

#ifdef _OPENMP
    DOUBLE stds[threads] = {};
    ULONGLONG count[threads] = {};
#endif /* _OPENMP */

    if (nchunks <= 1)
    {
        status = ERR_INVLD_SIZE;
    }

    if (status == SLM_SUCCESS)
    {
        DOUBLE *arr = new DOUBLE[(MAX_MASS * SCALE)+1];

        /* Trivial to include zeros or last value */
        arr[(MAX_MASS * SCALE)] = 0x0;
        arr[0] = 0x0;

        /* Compute the means */
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif /* _OPENMP */
        for (UINT i = 1; i < (MAX_MASS * SCALE) -1; i++)
        {
            arr[i] = 0x0;

            for (UINT j = 0; j < nchunks; j++)
            {
                arr[i] += (dslim.pepChunks[j].bA[i+1] - dslim.pepChunks[j].bA[i]);
            }

            arr[i] /= nchunks;

            if (arr[i] > 0)
            {
#ifdef _OPENMP
                count[omp_get_thread_num()]++;
#else
                truecount++;
#endif /* _OPENMP */
            }
        }

        /* Gather the counts to truecount */
#ifdef _OPENMP
        for (UINT kk = 0; kk < threads; kk++)
        {
            truecount += count[kk];
        }
#endif /* OPENMP */

        /* Compute the Standard Deviations */
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif /* _OPENMP */
        for (UINT i = 1; i < (MAX_MASS * SCALE) - 1; i++)
        {
            DOUBLE mean = arr[i];
            arr[i] = 0x0;

            if (mean > 0)
            {
                for (UINT j = 0; j < nchunks; j++)
                {
                    arr[i] += (((DOUBLE) dslim.pepChunks[j].bA[i+1]
                            - (DOUBLE) dslim.pepChunks[j].bA[i] - mean)
                            * ((DOUBLE) dslim.pepChunks[j].bA[i+1]
                            - (DOUBLE) dslim.pepChunks[j].bA[i] - mean));
                }

                arr[i] = sqrt(arr[i] / nchunks);
            }
        }

        /* Compute mean of stdevs */
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif /* _OPENMP */
        for (UINT i = 0; i < MAX_MASS * SCALE; i++)
        {
#ifdef _OPENMP
            stds[omp_get_thread_num()] += arr[i];
#else
            sum_std += arr[i];
#endif /* _OPENMP */
        }

/* Gather the counts to truecount */
#ifdef _OPENMP
        for (UINT kk = 0; kk < threads; kk++)
        {
            sum_std += stds[kk];
            stds[kk] = 0;
        }
#endif /* OPENMP */

        sum_std /= truecount;
        avg = sum_std;

        sum_std = 0;

        /* Compute stdev of stdevs */
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif /* _OPENMP */
        for (UINT i = 0; i < MAX_MASS * SCALE; i++)
        {
#ifdef _OPENMP
            stds[omp_get_thread_num()] += ((arr[i] - avg) * (arr[i] - avg));
#else
            sum_std += ((arr[i] - avg) * (arr[i] - avg));
#endif /* _OPENMP */
        }

        /* Gather the counts to truecount */
#ifdef _OPENMP
        for (UINT kk = 0; kk < threads; kk++)
        {
            sum_std += stds[kk];
        }
#endif /* OPENMP */

        sum_std /= truecount;
        sum_std = std::sqrt(sum_std);

        std = sum_std;

        delete[] arr;
        arr = NULL;
    }

    return status;
}

/*
 * FUNCTION: DSLIM_Deinitialize
 *
 * DESCRIPTION: Deallocate DSLIM memory
 *
 * INPUT: none
 *
 * OUTPUT:
 * @status: Status of execution
 */
STATUS DSLIM_Deinitialize(VOID)
{
    /* Deallocate all the DSLIM chunks */
    for (UINT chno = 0; chno < dslim.nChunks; chno++)
    {
        SLMchunk curr_chunk = dslim.pepChunks[chno];

        if (curr_chunk.bA != NULL)
        {
            delete[] curr_chunk.bA;
            curr_chunk.bA = NULL;
        }

        if (curr_chunk.iA != NULL)
        {
            delete[] curr_chunk.iA;
            curr_chunk.iA = NULL;
        }

        if (curr_chunk.sC != NULL)
        {
            delete[] curr_chunk.sC;
            curr_chunk.sC = NULL;
        }
    }

    if (dslim.pepChunks != NULL)
    {
        delete[] dslim.pepChunks;
    }

    /* Reset Global Variables */
    nchunks = 0;
    dslim.modChunks = 0;
    dslim.pepChunks = 0;
    dslim.nChunks = 0;
    chunksize = 0;
    lastchunksize = 0;

    return SLM_SUCCESS;
}
