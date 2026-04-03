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

# March 19, 2026 notes
##Further tasks
- Functional analysis (Picrust) is being run by Hans now due to CPU issues
- Not sure what's wrong, possibly data
- Taxonomic analysis to Family level (Isolate desulfobacteriae family)
## Beta + Alpha Diversity results
- Alpha diversity: Heatmap gives significantly lower diversity in high conductivity levels
- Beta diversity: p 0.045 on beta diversity, likely due to outliers
## Graphing notes
- Format should map general theme/colors chosen
- Do not use default
- Theme_classic() is good
## Manuscript
- Introduction written by Sam
- Can start on rest of manuscript as needed

# March 26, 2026
## Results discussion
Alpha diversity: 
- Multiple dots at same heights
- Increased Faith's at low conductivity = high conductivity leads to reduced overall diversity, improved chances of desulfobacteria specifically
- Remove figure legend

Beta diversity (Weighted UniFrac): 
- Not much difference
- Unweighted is preferred for Faith's 
- Low conductivity at top, high conductivity below (REVERSE THAT)
- Prevalence = 0.5
- Desulfobacteriae appears in 0.2? 
- Abundance = 0.01
- Normal values
- No immediate interesting data. Possibly interesting to see what families specifically are seen in each
- Ignore shared ones
- Desulfobacteria not present? 
- Placeholder for now, details need to be finagled

Taxonomy analysis:
- No need for rarefication, needs removal
- Palettes need to be chosen differently (qualitative and distinct rather than quantitative)
- Only desulfobulbaceae vs others 
- % in phylum + within phylum of family stated rather than the % barplot
- Plot average abundance of all samples instead of per sample? 
- Box plot with relative abundance (Y%) 
- Current barplot can be changed to supplemental figure (after fixing)
- fix the y axis 
- Wilcoxon analysis

ISA: 
- Phylum, Family, Genus, Indicator Value, P Value are main important statistics for the paper
- Functional analysis: 
- Waiting for response from someone who can actually run it :( 

Manuscript
- Phyla names used often - family in italics or not? Italics seems to be neat 

Story
- No major difference other than alpha diversity indicating lower diversity in high conductivity environments
- Possibly due to something else, not likely desulfobacteriae given other data
- They exist, difficult to say whether it contributes to conductivity 
- ISA shows one comes up, but not sure if specific species 
- Possible for conductivity high/low threshold being too high? 
- Not much can be said about the correlation between conductivity and desulfobacteriae
- Taxa barplot graphing only family/phylum of interest as a continuous variable compared to conductivity? Is there a relationship there? 
- Spearman's test for p value (Continuous vs continuous)
- Also gives rho number (positive or negative correlation) 
- Will consider both linear and parabolic/other shape graphs
- Line of best fit + ribbon 
- x/sum x function
- deseq also does conversion

Slides
- Possible looking at binned vs unbinned? Is it more meaningful to have one or both

# April 2 Notes
## Presentation Slides 

### Introduction 
- Slide 4: Start off by introducing soil micoribota, then go into WHY electron transport is IMPORTANT… 
- Add “there’s DECREASED oxidative reduction” processes in the soil so increased dependence of anaerobic reduction
- Emphasize that there’s LESS electron availability DEEPER down in the soil so we can emphasize the importance of electron transport
- Slide 5: More description for the FIGURE TEXT wise
- Make sure to remove, in author et. al
- NO FIGURE CAPTIONS!!!
- Hans advice: Don’t explain more than ONE concept per slide 
Slide 7: 
- Ensure to “add in” for WETLANDS in hypothesis
- Say “MORE wll studied in freshwater sediments” 

### Results:
- Slide 11: “No significant difference in Desulfobulbaceae family abundance between conductivity categories” should be the TITLE
- Add DIVERSITY metrics FIRST, BEFORE taxonomic analysis 
- For alpha, directly say “HIGH conductivity” has LOWER faith’s diversity as title
- Explain the reasoning for how we bin the high vs low (the threshold values) 
- Combine Alpha with Beta diversity results and combine the titles
- Mention how conductivity epxlains 5% of the data at BEST, which is very little 
- ONLY use the LIST of family names BESIDES the venn diagram, REMOVE the donuts
- Get rid of the SHARED info 
- We see the “methyl” family in BOTH high and low conductivity groups sinc ewe’re LOOKING at a LOWER taxonomic LEVEL 
- Looking at the BIG picture it seems to change, but DOES it change when it comes to the specific family/species level 
- For ISA analysis, put the low and HIGH conductiviteis together
- COLORCODE them so that the HIGH is one color and LOW is another
- On slide 13, recheck the figure… why are there so many points?

### Order of figures
- Do alpha/beta DIVERSITY metrics FIRST
- Then go to desulfobaubalcae taxonomic graph and show it’s NOT signfiaccnt
- Do indicator species, though which we found the family of interest (Geobacterioceae)
- Even though it’s at LOW conductivity, which contradicts what we hypohteiszse… but it’s ok since we can still discuss it 
- Now look more specifically at Geobacterioceae in taxonomic analysis  
- Core microbiome can be in the supplemental 

## MANUSCRIPT notes:

### Figure 1: Alpha as panel 1a, beta as panel 1b 
- At the diversity level, the binning DOES result in differences 
- LOW conductivity group has HIGHER diveristy

### Figure 2: Core microbiome 
- Put list of microbes, since it has the Geobacterioceae which appears as a UNIQUE species for the LOW conductivity group

### Table 1: ISA 

### Figure 3: Taxonomic/Unbinned
- Panel a: BOX plot for desulfulbabaeae 
- Panel b: Continuous correlation plot for the UNBINNED version for desulfulbabaeae
- Panel c: Boxplot of relative abundance for Geobacterioceae
- Panel d: Continuous correlation plot for the UNBINNED version for desulfulbabaeae 

### 1 supplemental figure 1
- Panel A: Taxa barplot for Desulfulbabaeae 
- Panel B: Taxa barplot Geobacterioceae

## TO DO:
- Chance to resubmit presentation if needed
-Meet TUESDAY 3:30pm at LSI to finalize graph since team 5 is presenting OUR group on THURSDAY

