# MVMA
Code, figures, and tables for paper "Multivariate meta-analysis with an increasing number of parameters" by Simina M. Boca, Ruth M. Pfeiffer, Joshua N. Sampson

Each directory corresponds to 1 or 2 figures or tables in the paper and provides the complete code to generate them.
Each directory also has a README.txt file with specific instructions. The directories corresponding to figures also have a subdirectory titled "figures" which contain the actual figures in the manuscript, which can be generated following the instructions in the README.txt files.

Figure 1, Table 3, Figure S2, and Figure S3 can be generated directly from the code provided and are not simulation-based. They are obtained by simply running the *.Rnw files in the respective directories.

Figures 2, S4, S5, S6, S7, and S8 are simulation based. The corresponding directories include a subdirectory called "simulation_code" which has both the R code and the code for running parallel simulations on a high-performance cluster. Note that these simulations can take several days to run. The make_Figure_2 directory also includes the combined simulation output in the "simResultsComb" subdirectory. This is not uploaded for the other figures, given storage limitations.
