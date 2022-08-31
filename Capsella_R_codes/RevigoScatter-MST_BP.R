# A plotting R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# --------------------------------------------------------------------------
# If you don't have the ggplot2 package installed, uncomment the following line:
# install.packages( "ggplot2" );
library( ggplot2 );

# --------------------------------------------------------------------------
# If you don't have the scales package installed, uncomment the following line:
# install.packages( "scales" );
library( scales );

# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0006631","fatty acid metabolic process",1.109,-4.616,-4.792,2.393,-5.409,0.633,0.000),
c("GO:0009555","pollen development",1.488,5.984,-1.593,2.520,-5.495,0.809,0.000),
c("GO:0009751","response to salicylic acid",0.685,4.061,5.035,2.185,-1.844,0.994,0.000),
c("GO:0055114","(obsolete) oxidation-reduction process",0.483,1.377,6.386,2.000,-2.234,1.000,0.000),
c("GO:0034504","protein localization to nucleus",0.289,3.624,-7.317,1.813,-1.945,0.964,0.033),
c("GO:0009691","cytokinin biosynthetic process",0.126,-4.912,0.975,1.462,-4.469,0.663,0.100),
c("GO:0006099","tricarboxylic acid cycle",0.261,-1.299,-8.259,1.771,-2.135,0.847,0.112),
c("GO:0044275","cellular carbohydrate catabolic process",0.334,1.036,-7.328,1.875,-1.646,0.783,0.114),
c("GO:1901141","regulation of lignin biosynthetic process",0.063,-1.961,5.828,1.176,-2.842,0.894,0.121),
c("GO:0006629","lipid metabolic process",4.900,-1.653,-5.369,3.037,-4.569,0.853,0.155),
c("GO:0050826","response to freezing",0.117,5.409,4.022,1.431,-1.693,0.994,0.159),
c("GO:0009308","amine metabolic process",0.518,-6.258,2.272,2.064,-2.996,0.873,0.181),
c("GO:0016575","histone deacetylation",0.185,0.038,0.706,1.623,-2.119,0.738,0.187),
c("GO:0043455","regulation of secondary metabolic process",0.289,-3.096,5.528,1.813,-1.620,0.936,0.223),
c("GO:0009699","phenylpropanoid biosynthetic process",0.383,-6.435,-1.252,1.934,-2.049,0.783,0.288),
c("GO:0032989","cellular component morphogenesis",0.212,4.270,-1.027,1.681,-4.886,0.807,0.346),
c("GO:0007568","aging",0.582,5.545,-2.992,2.114,-1.625,0.881,0.374),
c("GO:0048229","gametophyte development",1.961,6.012,-2.139,2.639,-4.215,0.861,0.471),
c("GO:0048869","cellular developmental process",3.426,5.532,-2.312,2.881,-1.817,0.825,0.479),
c("GO:0030261","chromosome condensation",0.086,1.598,1.321,1.301,-1.606,0.909,0.522),
c("GO:0010817","regulation of hormone levels",1.222,-4.288,4.057,2.435,-1.795,0.861,0.542),
c("GO:0044283","small molecule biosynthetic process",3.079,-5.442,-2.964,2.835,-1.843,0.695,0.598),
c("GO:0042761","very long-chain fatty acid biosynthetic process",0.095,-5.553,-3.484,1.342,-1.994,0.659,0.659),
c("GO:0000038","very long-chain fatty acid metabolic process",0.131,-5.257,-4.917,1.477,-2.022,0.681,0.678));

one.data <- data.frame(revigo.data);
names(one.data) <- revigo.names;
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
one.data$plot_size <- as.numeric( as.character(one.data$plot_size) );
one.data$value <- as.numeric( as.character(one.data$value) );
one.data$frequency <- as.numeric( as.character(one.data$frequency) );
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );
#head(one.data);


# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below

p1 <- ggplot( data = one.data );
p1 <- p1 + geom_point( aes( plot_X, plot_Y, colour = value, size = plot_size), alpha = I(0.6) ) + scale_size_area();
p1 <- p1 + scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(one.data$value), 0) );
p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = plot_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) + scale_size_area();
p1 <- p1 + scale_size( range=c(5, 30)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
ex <- one.data [ one.data$dispensability < 0.15, ];
p1 <- p1 + geom_text( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 3 );
p1 <- p1 + labs (y = "semantic space x", x = "semantic space y");
p1 <- p1 + theme(legend.key = element_blank()) ;
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10);
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/10,max(one.data$plot_Y)+one.y_range/10);


# --------------------------------------------------------------------------
# Output the plot to screen

p1;

# Uncomment the line below to also save the plot to a file.
# The file type depends on the extension (default=pdf).

# ggsave("/path_to_your_file/revigo-plot.pdf");

