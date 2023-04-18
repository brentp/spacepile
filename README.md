# spacepile (WIP)

`spacepile` efficiently creates (multi) alignments from replicate measurements of the same piece of DNA as
from **duplex** sequencing, single molecule consensus, pac-bio **CCS**, or nanopore [**cyclomics**](https://cyclomics.com/).
`spacepile` accepts alignments (cigars) of grouped sequences and _spaces_ them into matricies that can then be sent to any deep-learner.

Note that this is done nicely in [deepconsensus](https://github.com/google/deepconsensus) which is a major source
of inspiration for this work, but that is specific to PacBio and (likely) less efficient at extraction.

`spacepile` **will** also provide means to facilitate training (agnostic to the learner) by providing tools to
sample and augment reads and to separate training from test datasets.

Currently `spacepile` can be used by downloading and `pip install`ing the python wheel from the releases.

```Python
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
        [(pysam.CMATCH, 4)],                                     # read 1. 4M     ACTG
        [(pysam.CMATCH, 2), (pysam.CDEL, 2), (pysam.CMATCH, 1)], # read 2. 2M2D1M  CTC
        [(pysam.CMATCH, 2), (pysam.CINS, 3), (pysam.CMATCH, 2)], # read 3. 2M3I2M ACGGGTG
        ]

posns = [0, 1, 0] # the 2nd read starts 1 base after the other reads.

# idxs will be filled by `space` call.
idxs = np.zeros((len(cigs), max_width), dtype=np.uint32)
spacepile.space(idxs, cigs, posns)

#  now we have the indexes into the original sequence or base-qualities or IPDs. etc.
SPACE = np.iinfo(np.uint32) - 1
END = np.iinfo(np.uint32)
assert np.array_equal(idxs,
                [[  0,    1,SPACE,SPACE,SPACE,    2,    3,END  ,END],
                 [END,    0,SPACE,SPACE,SPACE,    1,SPACE,SPACE,  2],
                 [  0,    1,    2,    3,    4,    5,    6,END  ,END]])


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

# now repeat for many sequences and send many mats to a learner.
```

# Research

## steps

### training

#### read selection

Need enough unique reads from across genome to train. Given duplex-seq is often deep,
we want to sample wisely so we're not training on the same position 10,000 times.

#### calibration

[Paper covering this topic](https://arxiv.org/abs/1706.04599)

https://scikit-learn.org/stable/modules/calibration.html
https://arxiv.org/abs/2210.00045 ( from google 9/30/2022)

NOTE: this one seems best. Just get pytorch model into sklearn-like class.
https://github.com/dirichletcal/dirichlet_python

See this python module: https://github.com/classifier-calibration/PyCalib

The output from the model is not calibrated -- the error rates do not match the probabilities.
This is a known problem.
There is some work on this here:
https://github.com/gpleiss/temperature_scaling (needs logits, not softmax)

https://github.com/dirichletcal/experiments_dnn

But this may require further work because of the positional info where the ends of reads
often have lower quality.

### spacing

This implements the multi-alignment spacing which takes a set of reads and returns
something that is sent to the network. Currently this is implemented in python.
A rust version is in progress.

This involves clustering (by tag) and selection to max_depth so we can have a fixed
network size.

### consensus

Take a trained model and a bam and write a new fastq with well-calibrated base-qualities.

## build notes

```
docker run --rm -v $(pwd):/io ghcr.io/pyo3/maturin build --release
pip install target/wheels/spacepile-0.1.0-cp37-abi3-linux_x86_64.whl
```
