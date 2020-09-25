# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from multiprocessing import Pool, cpu_count
from q2_types.feature_data import (DNAFASTAFormat, AlignedDNAIterator)
from skbio import DNA

CHUNK_LEN = 75


def worker(seq_obj):
    '''Degap a single sequence and place the result in a queue'''
    dg_seq = seq_obj.degap()
    return dg_seq


def degap_seqs(aligned_sequences:
               AlignedDNAIterator,
               min_length: int = 1) -> DNAFASTAFormat:
    result = DNAFASTAFormat()

    pool = Pool(4)
    with result.open() as out:
        for seq in aligned_sequences:
            seq_chunks = [seq[i:i + CHUNK_LEN] for i in
                          range(0, len(seq), CHUNK_LEN)]
            dg_seqs = pool.imap(worker, seq_chunks)
            dg_seq = DNA.concat(dg_seqs)
            dg_seq.metadata = seq.metadata
            if len(dg_seq) >= min_length:
                dg_seq.write(out)

    pool.close()
    pool.join()

    return result
