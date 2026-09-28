# COBRA
<img width="1000" height="400" alt="Logo" src="https://github.com/user-attachments/assets/b37f544f-85bb-4f51-8fc0-470f3a6bbcdf" />

**View contig read coverages and associated gene annotations in an easy, user-friendly, and publication-quality figure format.**

## Quick Start:

1. Open [COBRA](https://jlmaier12.github.io/COBRA/).
2. Load your [pileup](https://bbmap.org/tools/pileup) file.
3. Load your [gff](https://en.wikipedia.org/wiki/General_feature_format) file.
   
Note- These files tend to be large and will take a minute to load and process.

4. Select the reference name for the contig you'd like to view. Generally, it is a good idea to have some contigs of interest in mind before using COBRA.
5. (Optional) Choose the base pair range to subset on the contig. Gene annotations will only be visible on a <30 kbp subset.
6. (Optional) Type a specific gene annotation to highlight it on the subset plot.
7. (Optional) Put read coverage in log scale.

Example usage:

<img width="3809" height="1604" alt="image" src="https://github.com/user-attachments/assets/574a8bbe-7e58-4675-88bf-8862a42d4d8d" />

## Acknowledgements

COBRA relies heavily on the functions of many, many tools and packages, but I'd like to especially thank the developers of [gggenes](https://github.com/wilkox/gggenes) which is an **excellent** stand-alone tool for visualizing gene annotations. COBRA simply pairs the gggenes functionality with read coverage visualization making for a more comprehensive plot overall.
