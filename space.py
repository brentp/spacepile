import spacepile
print(spacepile)
import numpy as np
import pysam


# make matricies of at most this width. usually this will be ~ 10% longer than read-length
max_width = 9
# NOTE, these would likely be drawn from a bam file that is sorted by fragment (an option in samtools sort) or grouped by UMI
# these are available as aln.cigartuples in pysam.
# here we show 3 reads for an example.
cigs = [
        [(pysam.CHARD_CLIP, 100), (pysam.CMATCH, 4)],                                     # read 1. 4M     ACTG
        [(pysam.CMATCH, 2), (pysam.CDEL, 2), (pysam.CMATCH, 1)], # read 2. 2M2D1M  CTC
        [(pysam.CMATCH, 2), (pysam.CINS, 3), (pysam.CMATCH, 2)], # read 2. 2M2D1M ACGGGTG
        ]

posns = [0, 1, 0] # the 2nd read starts 1 base after the other reads.

# idxs will be filled by `space` call.
idxs = np.zeros((len(cigs), max_width), dtype=np.uint16)
spacepile.space(idxs, cigs, posns)

#  now we have the indexes into the original sequence or base-qualities or IPDs. etc.
assert np.array_equal(idxs,
                [[    0,    1,65534,65534,65534,    2,    3,65535,65535],
                 [65535,    0,65534,65534,65534,    1,65534,65534,    2],
                 [    0,    1,    2,    3,    4,    5,    6,65535,65535]]), idxs


# sequences likely retrieved from pysam's aln.query_sequence or aln.query_alingment_sequence
raw_sequences = ["ACTG", "CTC", "ACGGGTG"]
base_qs =       [[ 40, 50, 60, 70], [60, 60, 60], [45, 45, 45, 55, 55, 55, 60]]

sequences = np.zeros((len(cigs), max_width), dtype=np.int16)
bqs = np.zeros_like(sequences)
for i, s in enumerate(raw_sequences):
    sequences[i, :len(s)] = np.array(list(s), dtype='U1').view(np.int32) # Q: is there a better way to do this in numpy?
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

#print(mat)
# now repeat for many sequences and send many mats to a learner.

# MIDNSHP=XB
cigs = [
        #[(0, 172-126), (4, 930)], # 1040S 172M 930S
        #[(0, 136-126), (2, 2), (0, 34), (5, 1767)],
        [(4, 1040), (0, 172), (4, 930)], # 1040S 172M 930S
        [(5, 205), (0, 136), (2, 2), (0, 34), (5, 1767)],
        [(5, 628), (0, 39), (2, 1), (0, 38), (1, 1), (0, 82), (5, 1354)],
        [(5, 1899), (0, 36), (2, 1), (0, 72), (1, 2), (0, 27), (5, 106)],
        [(5, 1487), (0, 40), (1, 1), (0, 99), (5, 515)]
        ]
posns = [x + 1040 for x in [0, 0, 8, 21, 33]][:len(cigs)]
print(posns)

max_width = 177
idxs = np.zeros((len(cigs), max_width), dtype=np.uint16)
spacepile.space(idxs, cigs, posns)
print(idxs)

cigs = [
        [(4, 799), (0, 16), (1, 2), (0, 23), (2, 1), (0, 90), (4, 12)], # 799S 16M 2I 23M 1D 90M 12 S
        [(5, 384), (0, 57), (2, 2), (0, 56), (5, 445)],                 # 384H 57M 2D 56M 445H
        ]
posns = [21271006, 21271021]
max_width = 177
idxs = np.zeros((len(cigs), max_width), dtype=np.uint16)
spacepile.space(idxs, cigs, posns)
print(idxs)

cigs = [
        [(5, 239), (0, 20), (2, 2), (0, 35), (5, 648)],
        [(5, 672), (0, 54), (2, 1), (0, 15), (1, 1), (0, 12), (5, 188)],
        [(5, 298), (0, 23), (1, 1), (0, 3), (2, 2), (0, 58), (1, 1), (0, 5), (1, 1), (0, 29), (5, 523)]]
posns = [ 0, 12, 27]
for max_width in range(10, 100, 10):
    idxs = np.zeros((len(cigs), max_width), dtype=np.uint16)
    spacepile.space(idxs, cigs, posns)
    print(idxs)
