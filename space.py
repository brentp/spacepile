import spacepile
import numpy as np
import pysam

a = np.zeros((1, 5), dtype=np.uint16)
cigs = [[(int(pysam.CMATCH), 5)]]

posns = [0]
spacepile.space(a, cigs, posns)
print(a)
