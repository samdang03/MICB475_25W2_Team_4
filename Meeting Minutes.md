# February 5, 2026 Notes

## Variables 

Conductivity (Need binning to determine a threshold for high vs low), Exchange capacity (Also need binning to determine threshold), Soil_NH4, Soil_CH4

## Troubleshoot to verify predicted inverse relationship between conducitivty/exchange capacity and methane/ammonium production 

1.) Correlation graph:
See if there's a negative relationship between conducitvity and NH4, and test for signficance
conductivity vs NH4 
conductivity vs CH4 
exchange capacity vs NH4 
exchange capacity vs CH4
Conductivity vs exchange

2.) Diversity metrics:
First extract taxonomic info using qiime2, then look at high vs low comparisons of conductivity and exchange capacity 

3.) Plot taxa barplot:
To see if there are certain taxa more representative or not

4.) Relative abundance
After you have filtered out cyanobacteria
Get a boxplot that tells you whether one is higher or lower 
Hopefully see bacterial relatively more abundant in high or low conditions

5.) Functional analysis 
Not taught in course but may have to research how to do it on our own 
Module 19 teaches you 
Gives u predictive functional profile on your groups 

## Extra notes 
Start uploading your metadata on the shared server
MAKE SURE you email Hans before trimming data 

# Feburary 12, 2026 Notes

1.) Further analysis to determine the exact "threshold" (step behavior) in conductivity (as compared to net mineralization) needed

2.) Moving away from a high-throughput plot, focusing on two selected variables going forwards after ensuring consistency

3.) Net immobilization may be due to high NH4 levels in general and may not be entirely linked to cable bacteria

## Extra notes: Proposal question should be drafted


# Feburary 19, 2026 Notes
- Title should not be a paper title yet.. should be more open for what you are planning to do on your proposal eg. Investigation... or Exploring...

#March 5, 2026 Notes

## Proposal comments
- Not enough on the BACTERIA itself… want to focus on diversity/abundance
- Need more CONTEXT to appear in the INTRODUCTION…
- In the aims, we have these analysis to see what’s up with the bacteria
- For aim 1 (diversity analysis), be specific on what metric to do
- For Taxonomic analysis, do CORE MICROBIOME as well, not just ISA, because it shows you SHARED vs UNIQUE taxa 
- To do this, we need the phyloseq object which combines 4 things from qiime2 
- Eric made this already, she just needs to push it 
- NEW proposal is good as it focuses on the bacteria 

## Manuscript writing
- Plan to put the “Making sense of what we got” stuff into the DISUCSSION part of our manuscript
- The introduction should be an overview, a general summary of the data, as well as KEY findings of the paper 

## Presentation
- On someone else’s project, either on April 7 or April 9 (8-9:30AM)
- Need to submit your SLIDES by the FIFTH of april (sunday) 
- Have FOUR meetings until that day, excluding today 
- TWO more meetings from now to generate data and discuss them
- Thursday, March 12th and 19th will be generating plots and showing Hans
- On the 26th, we want to have all the data DONE, and think about the complete story to present on the 7th or 9th
- After the 26th, Go home and make the slide deck with the narrative by the SECOND of april 
Have 3 days to make edits based on suggestions 

## Data analysis
PLOT everything on R, including taxaonomic barplot 


# March 12, 2026 notes

## Diversity analysis: Currently run by Ruth and Rick, will probably be done by next week
- Run in R 
- Taxonomic analysis: Current run by Erica, files are in "results" on GitHub
- Trim phyla down to top 10, others clustered as "other"
- Possible further family runs under phylum on interest (Thermodesulfobacteriota) 
## Functional (pathway) analysis: 
- Medici 
- K-O recommended to run
- Looking for DNRA + Esox as the pathway of interest
- DESEQ heatmap + barplot generation after PiCrust2
## Cleaning Github
- How to subfolder
- .rm-r Folder, NOT rm folder password required. Add f after folder to force remove


