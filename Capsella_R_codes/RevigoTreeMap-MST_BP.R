# A treemap R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","freqInDbPercent","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0006631","fatty acid metabolic process",1.109,5.409,0.633,0.000,"fatty acid metabolic process"),
c("GO:0009691","cytokinin biosynthetic process",0.126,4.469,0.663,0.100,"fatty acid metabolic process"),
c("GO:0006099","tricarboxylic acid cycle",0.261,2.135,0.847,0.112,"fatty acid metabolic process"),
c("GO:0044275","cellular carbohydrate catabolic process",0.334,1.646,0.783,0.114,"fatty acid metabolic process"),
c("GO:1901141","regulation of lignin biosynthetic process",0.063,2.842,0.894,0.121,"fatty acid metabolic process"),
c("GO:0006629","lipid metabolic process",4.900,4.569,0.853,0.155,"fatty acid metabolic process"),
c("GO:0009308","amine metabolic process",0.518,2.996,0.873,0.181,"fatty acid metabolic process"),
c("GO:0016575","histone deacetylation",0.185,2.119,0.738,0.187,"fatty acid metabolic process"),
c("GO:0043455","regulation of secondary metabolic process",0.289,1.620,0.936,0.223,"fatty acid metabolic process"),
c("GO:0009699","phenylpropanoid biosynthetic process",0.383,2.049,0.783,0.288,"fatty acid metabolic process"),
c("GO:0030261","chromosome condensation",0.086,1.606,0.909,0.522,"fatty acid metabolic process"),
c("GO:0010817","regulation of hormone levels",1.222,1.795,0.861,0.542,"fatty acid metabolic process"),
c("GO:0044283","small molecule biosynthetic process",3.079,1.843,0.695,0.598,"fatty acid metabolic process"),
c("GO:0042761","very long-chain fatty acid biosynthetic process",0.095,1.994,0.659,0.659,"fatty acid metabolic process"),
c("GO:0000038","very long-chain fatty acid metabolic process",0.131,2.022,0.681,0.678,"fatty acid metabolic process"),
c("GO:0009555","pollen development",1.488,5.495,0.809,0.000,"pollen development"),
c("GO:0032989","cellular component morphogenesis",0.212,4.886,0.807,0.346,"pollen development"),
c("GO:0007568","aging",0.582,1.625,0.881,0.374,"pollen development"),
c("GO:0048229","gametophyte development",1.961,4.215,0.861,0.471,"pollen development"),
c("GO:0048869","cellular developmental process",3.426,1.817,0.825,0.479,"pollen development"),
c("GO:0009751","response to salicylic acid",0.685,1.844,0.994,0.000,"response to salicylic acid"),
c("GO:0050826","response to freezing",0.117,1.693,0.994,0.159,"response to salicylic acid"),
c("GO:0055114","(obsolete) oxidation-reduction process",0.483,2.234,1.000,0.000,"(obsolete) oxidation-reduction process"),
c("GO:0034504","protein localization to nucleus",0.289,1.945,0.964,0.033,"protein localization to nucleus"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
p2 <- treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = "REVIGO TreeMap",
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()

