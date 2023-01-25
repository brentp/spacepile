import spacepile
print(spacepile)
import numpy as np
import pysam


# make matricies of at most this width. usually this will be ~ 10% longer than read-length
max_width = 9
# NOTE, these would likely be drawn from a bam file that is sorted by fragment (an option in samtools sort)
# these are available as aln.cigartuples in pysam.
cigs = [
        [(pysam.CMATCH, 4)],                                     # read 1. 4M     ACTG
        [(pysam.CMATCH, 2), (pysam.CDEL, 2), (pysam.CMATCH, 1)], # read 2. 2M2D1M  CTC
        [(pysam.CMATCH, 2), (pysam.CINS, 3), (pysam.CMATCH, 2)], # read 2. 2M2D1M ACGGGTG
        ]

posns = [0, 1, 0]

idxs = np.zeros((len(cigs), max_width), dtype=np.uint16)
spacepile.space(idxs, cigs, posns)

#  now we have the indexes into the original sequence or base-qualities or IPDs. etc.
assert np.array_equal(idxs,
                [[    0,    1,65534,65534,65534,    2,    3,65535,65535],
                 [65535,    0,65534,65534,65534,    1,65534,65534,    2],
                 [    0,    1,    2,    3,    4,    5,    6,65535,65535]])


# sequences likely retrieved from pysam's aln.query_sequence or aln.query_alingment_sequence
raw_sequences = ["ACTG", "CTC", "ACGGGTG"]
base_qs =       [[ 40, 50, 60, 70], [60, 60, 60], [45, 45, 45, 55, 55, 55, 60]]

sequences = np.zeros((len(cigs), max_width), dtype=np.int16)
bqs = np.zeros_like(sequences)
for i, s in enumerate(raw_sequences):
    sequences[i, :len(s)] = np.array(list(s), dtype='U1').view(np.int32)
    bqs[i, :len(s)] = np.array(base_qs[i], dtype=np.int32)


# NOTE we will put sequences and bqs in the same output matrix.
mat = np.zeros((2 * len(cigs), max_width), dtype=np.int16)
spacepile.translate(idxs, sequences, mat[:len(cigs)])

# NOTE how the columns generally match
assert np.array_equal(mat[:len(cigs)],
                   [[65,67,-2,-2,-2,84,71,-1,-1],
                    [-1,67,-2,-2,-2,84,-2,-2,67],
                    [65,67,71,71,71,84,71,-1,-1]])


spacepile.translate(idxs, bqs, mat[len(cigs):])
assert np.array_equal(mat[len(cigs):],
			[[40,50,-2,-2,-2,60,70,-1,-1],
			 [-1,60,-2,-2,-2,60,-2,-2,60],
			 [45,45,45,55,55,55,60,-1,-1]])


# now repeat for many sequences and send many mats to a learner.
