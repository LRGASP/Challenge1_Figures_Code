# Challenge1_Figures_Code

Code to generate Challenge 1 figures for LRGASP paper.  The data is available from Sage Synapse. 

To generate the figures, obtain code from GitHub:

```
git clone git@github.com:LRGASP/Challenge1_Figures_Code.git
cd Challenge1_Figures_Code
```

Download the data files from:

```
https://cgl.gi.ucsc.edu/data/LRGASP/paper/Challenge1_Figures_Data.zip
```

Run the R programs to build figures into the `output` directory:

```
unzip -q Challenge1_Figures_Data.zip
Rscript Code_Figure_Challenge1.R
Rscript Code_SupplementaryFigures_Challenge1.R
```



