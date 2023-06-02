# Challenge1_Figures_Code

Code to generate Challenge 1 figures for LRGASP paper.  The data is available from Sage Synapse. 

To generate the figures, obtain code from GitHub:

```
git clone git@github.com:LRGASP/Challenge1_Figures_Code.git
cd Challenge1_Figures_Code
```

The data file is `Challenge1_Figures_Data.zip`, Synapse id  `syn51602848`.
Install synapseclient if necessary: `pip install synapseclient`

Download, extract, and run the R programs to build figures into the `output` directory:

```
synapse get syn51602848
unzip -q Challenge1_Figures_Data.zip
Rscript Code_Figure_Challenge1.R
Rscript Code_SupplementaryFigures_Challenge1.R
```



