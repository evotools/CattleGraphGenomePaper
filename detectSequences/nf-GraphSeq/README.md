# Get non-reference sequences from a pg graph
## Introduction
The different script are run in a nextflow workflow to extract the sequences in the graph that are not present into any of the reference paths.

## Dependencies
The workflow will install all the dependencies using `environment.yml` file. Will then run all the processes in sequence, to identify, filter, combine and simplify the resulting regions.
Finally, it will perform a detection of novel genes in the regions.

## Notes
The workflow still needs to be tested extensively. However, in the legacy folder we provide the scripts used to perform the analysis, numbered by stage.