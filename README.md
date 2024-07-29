# Extra quality control of Mad4hatter libraries and sequencing runs
The objective of this script is to visualize the performance of Mad4hatter libraries in a sequencing run.
# Summary
# 1. Inputs
The script needs two inputs for every sequencing run:
Input 1: Sample Number of Loci with 100 Reads or More per Pool. This table is in the QCplots.html file generated as an output of the Mad4hatter pipeline.
Input 2: Sequencing excel file of the run in xlsx format.
# 2. Parasitemia and libraries performance
This section will generate four different plots:
A) Density plot of parasitemias
B) Density plot of loci with 100 reads or more
C) Correlation between loci with 100 reads or more and parasitemia
D) Box plot of loci with 100 reads or more per library plate
# 3. Heatmaps
This section will generate a set of four heatmaps for every library plate. The heatmap represents the performance of every well in the experiment for P1A, P5, P2 and the sum of them all.
