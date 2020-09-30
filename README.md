# Microbial Evolution

To define the parameters for a simulation, create a `.py` file (eg, `sim.py`) using
 `sample_experiment.py` as a template.

To execute the simulation `sim.py`, run the following command
```
$ python3 sim.py <outfile>
```
which will save (pickle) the simulation results in a file called `<outfile>`.

To merge the replicates (ie, multiple runs with the same set of parameters) of a simulation, run the
 following command
```
$ python3 merge.py <dirname> <bins>
```
where `<dirname>` is the name of the directory that contains the `.pkl` files for the replicates
 and `<bins>` specifies the number of bins to use to compute the host/virus genotype and host
  mass distributions. The command will produce a `summary.pkl` file under `<dirname>` summarizing
   the results across all the replicates.
  
To generate visualizations, run the following command within the directory that contains the
 `summary.pkl` file.
```
$ python3 visualizations.py summary.pkl
```

## Software Dependencies

- Python 3
- NumPy
- Matplotlib
- Dill

## Contact Us

If you have any questions about the software, please email <swami.iyer@gmail.com>.