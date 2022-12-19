+ use this for nd arrays: https://docs.rs/numpy/latest/numpy/
  https://github.com/PyO3/rust-numpy
  https://itnext.io/how-to-bind-python-numpy-with-rust-ndarray-2efa5717ed21

+ maturin: https://github.com/PyO3/maturin


# steps

## training

### read selection

Need enough unique reads from across genome to train. Given duplex-seq is often deep,
we want to sample wisely so we're not training on the same position 10,000 times.

### calibration

[Paper covering this topic](https://arxiv.org/abs/1706.04599)

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

## spacing

This implements the multi-alignment spacing which takes a set of reads and returns
something that is sent to the network. Currently this is implemented in python.
A rust version is in progress.

This involves clustering (by tag) and selection to max_depth so we can have a fixed
network size.

## consensus

Take a trained model and a bam and write a new fastq with well-calibrated base-qualities.
