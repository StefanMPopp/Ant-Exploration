Data & Code to "Ant colonies explore novel environments with increased activity and slower, curvier walks." Popp & Dornhaus 2023.
It takes in cleaned track data (from trackProcessing), processes them, analyzes them statistically, and produces publication-ready figures.

For just recreating the statistical analysis (note: it takes a few hours on a 16core machine and requires >16GB RAM due to the large data sets):

Open the R script
Change the fileDir to where your AntExploration folder sits.
Run the script
For just recreating the analysis data and figures in MATLAB:

Download the folder to your MATLAB path (usually .../Documents/MATLAB)
Open exploration.m
Change the one line of user input if wished (just running the figure code or recreating the analysis data, which takes a few minutes)
Run the rest of the script
In case of errors, make sure the required toolboxes are installed (indicated in the script description)
