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
revigo.data <- rbind(c("GO:0009926","auxin polar transport",0.302,2.889,0.730,0.000,"auxin polar transport"),
c("GO:1903338","regulation of cell wall organization or biogenesis",0.167,1.463,0.826,0.140,"auxin polar transport"),
c("GO:0051052","regulation of DNA metabolic process",0.388,1.304,0.745,0.183,"auxin polar transport"),
c("GO:0050794","regulation of cellular process",23.844,2.001,0.660,0.279,"auxin polar transport"),
c("GO:0007165","signal transduction",8.038,1.837,0.604,0.492,"auxin polar transport"),
c("GO:0006355","regulation of transcription, DNA-templated",10.869,1.779,0.584,0.533,"auxin polar transport"),
c("GO:0010817","regulation of hormone levels",1.222,1.978,0.735,0.588,"auxin polar transport"),
c("GO:0050789","regulation of biological process",26.558,1.470,0.683,0.618,"auxin polar transport"),
c("GO:0009889","regulation of biosynthetic process",13.268,1.508,0.613,0.696,"auxin polar transport"),
c("GO:0023052","signaling",8.205,1.816,1.000,0.000,"signaling"),
c("GO:0042542","response to hydrogen peroxide",0.316,1.795,0.897,0.000,"response to hydrogen peroxide"),
c("GO:0051716","cellular response to stimulus",13.795,1.809,0.861,0.241,"response to hydrogen peroxide"),
c("GO:1901698","response to nitrogen compound",1.352,1.530,0.886,0.346,"response to hydrogen peroxide"),
c("GO:0009605","response to external stimulus",8.178,1.330,0.903,0.416,"response to hydrogen peroxide"),
c("GO:0010200","response to chitin",0.591,1.719,0.866,0.453,"response to hydrogen peroxide"),
c("GO:0009733","response to auxin",1.785,1.353,0.874,0.486,"response to hydrogen peroxide"),
c("GO:1901700","response to oxygen-containing compound",7.046,1.347,0.862,0.522,"response to hydrogen peroxide"),
c("GO:0050896","response to stimulus",27.563,1.681,1.000,0.000,"response to stimulus"),
c("GO:0065007","biological regulation",29.718,1.327,1.000,0.000,"biological regulation"),
c("GO:0071897","DNA biosynthetic process",0.410,3.658,0.748,0.000,"DNA biosynthetic process"),
c("GO:0043543","protein acylation",0.410,1.280,0.859,0.210,"DNA biosynthetic process"),
c("GO:0016573","histone acetylation",0.239,1.406,0.767,0.307,"DNA biosynthetic process"),
c("GO:0019438","aromatic compound biosynthetic process",4.247,1.543,0.784,0.407,"DNA biosynthetic process"),
c("GO:0042401","cellular biogenic amine biosynthetic process",0.212,1.525,0.772,0.432,"DNA biosynthetic process"),
c("GO:1901362","organic cyclic compound biosynthetic process",4.761,1.342,0.796,0.536,"DNA biosynthetic process"),
c("GO:0018130","heterocycle biosynthetic process",3.904,1.373,0.783,0.548,"DNA biosynthetic process"),
c("GO:0034654","nucleobase-containing compound biosynthetic process",3.007,1.400,0.723,0.693,"DNA biosynthetic process"),
c("GO:0007154","cell communication",9.562,1.829,0.962,0.043,"cell communication"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
#pdf( file="revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
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

#dev.off()

