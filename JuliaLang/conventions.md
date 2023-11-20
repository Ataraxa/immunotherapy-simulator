# Conventions of file format

## CSV files for batch trajectories
These files should contain a matrix.
In row direction: different experiments
In col direction: time 

<u>Important:</u>
- The data should be raw, so not log'ed (i.e., in mm³) 
- The matrix should not be sliced in time dimension (i.e., do not discard the unused timepoints)