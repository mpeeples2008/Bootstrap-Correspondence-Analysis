# Bootstrap-Correspondence-Analysis
R-script for conducting and evaluating bootstrapped correspondence analyses following methods published by Peeples and Schachner (2012) in the Journal of Archaeological Science

This code relates to the following article:

Peeples, Matthew A. and Gregson Schachner
2012 Refining correspondence analysis-based ceramic seriation of regional data sets. Journal of Archaeological Science 39(8):2818-2827. 
http://dx.doi.org/10.1016/j.jas.2012.04.040

Online version of these instructions with visuals available at http://www.mattpeeples.net/caboot.html

Abstract
Ceramic seriation is one of the primary tools used by archaeologists to create chronologies, especially for surface survey data. In this article, we outline a bootstrapped approach to correspondence analysis-based seriation designed to help assess and improve the stability of relative orderings produced through such analyses. This procedure systematically identifies and removes small samples and sites with unusual samples, such as those with multiple components, which are not handled well by seriation and require secondary interpretation. Our approach combines data from multiple projects in the Zuni region of the American Southwest in order to gauge the effects of intraregional variation in ceramic distributions and reexamine prior interpretations of demographic change based on individual projects. This analysis highlights two previously unknown regional trends, including variation in the distribution of multicomponent sites and potential short-term depopulation of a significant part of the region. Based on this example, we further suggest that our method may be particularly useful in situations where little prior chronological knowledge is available.

Bootstrapped Correspondence Analysis:
This document provides a brief overview of the ca.boot.R script described in the article cited above. You can download a sample data file along with the script to follow along with this example. 

File Format:
This script is designed to use the *.csv (comma separated value) file format. 

Table Format:
Tables should be formatted with each of the samples/observations as rows and each of the variables to be included as columns. The first row of the spreadsheet should be a header that labels each of the columns. The first column should contain the name of each unit (i.e., level, unit, site, etc.). Row names may not be repeated. All of the remaining columns should contain numerical count data. This analysis will not work if there are missing data in any rows or columns, so samples with missing data should be removed before running the script. 

Requirements for Running the Script:
In order to run this script, you must install the R statistical platform (> version 2.13) and install the "ca" package.

Starting the Script:
The first step for running the script is to place the script file "ca.boot.R" and the data file "CA_data.csv" (or whatever you named your file) in the current working directory of R. To change the working directory, click on "File" in the R window and select "Change dir", then simply browse to the directory that you would like to use as the working directory. Next, to actually run the script, type the following line into the R command line:

source('ca.boot.R')

Running the Script:
After typing the command above into the command line, the console will automatically bring up a dialog box asking you to select the file on which you want to run this script. To follow along with the example here, choose the sample file CA_data.csv (right click and Save As to download). Once you have selected the file to use, a series of figures will appear in the console. To progress, simply click on the monitor window. 

1) The first figure to appear displays a plot of the first two dimensions for all data included in the original table. 

2) The next figure to appear shows the first two dimensions of the CA again, but with convex hulls representing the spread of 1,000 bootstrapped replicates for each site and type point.

3) Next, the console will display a histogram of the perimeter lengths for all site convex hulls. The dotted red line represents the mean perimeter length and the solid red line represents one standard deviation above the mean (by default the value for excluding points). 

4) Next, the console will display a histogram of the maximum dimension of all site convex hulls. The dotted red line represents the mean maximum dimension and the solid red line represents one standard deviation above the mean (by default the value for excluding points).

5) The next figure displays the first two dimensions of the CA for those sites that were retained through the site removal procedure. 

6) Finally, the console displays a plot of the first two dimensions of the CA on the reduced data set with convex hulls representing the spread of 1,000 bootstrapped replicates for each site and type point. By default the convex hulls on this final plot represent 90% convex hulls (i.e., they contain at least 90% of all bootstrapped replicates).
