
# The DBNN Software

Latest(2010) release of DBNN Software. may be downloaded from here.
The DBNN is a Bayesian classifier for machine learning applications. The software implementation works on Linux platforms. The latest version, a rebuild of the original with bug fixes and enhancements is available as autodbnn2.cpp in the zip file. To compile the code simply type in a terminal : c++ -O3 autodbnn2.cpp -o autodbnn2

## A brief description of DBNN training and testing procedure

The DBNN classifier requires training data arranged in a matrix-like format. Each row in the matrix should represent space delimited features of an example followed by its class label as the last column entry. All feature and class labels should be numeric. The algorithm requires two sets of parameter files. The first set defines the data and the second a set of input parameters that are used when the program executes. There are three parameter files relating to the data and they differ from the data file only in their extensions. So if the data file is named DBNN.dat, the parameter files should be DBNN.apf, DBNN.awf and DBNN.inf respectively. DBNN.awf is created during the training process. Both apf and inf files are to be created by the user. The inf file will be used in read only mode while apf file will be appended with additional data during the training process. All configuration file data are in text format. The format of the inf file is as follows: each parameter in the file should be specified in a new row. The first parameter should be the number of features in the training data, followed by the number of output classes defined in the data on the second line. The third line in the file is the allowed $\pm$ margin between the classes and is set to 0.45. This means anything within this margin from a class label will be treated as belonging to that class. This is significant only in regression problems as classification problems give the precise value for the class of the object. The next $m$ lines in the file are the labels assigned to the $m$ classes. They should be numeric, but need not be consecutive numbers. There are two more optional lines in the file. The first defines the lowest confidence value a prediction should have to qualify for classification. Any object that gets a confidence less than this value will be rejected by the classifier as not classifiable. This is a useful feature to increase the accuracy of the predictions if completeness is not the concern. The last entry in the file is the maximum number of bins into which a feature may be sliced. If the last two entries are missing, a default value is set by the classifier. But the other entries in the inf file are mandatory. Only one line of data is required in the apf file. For a data with $m$ features, the apf file should have $m$ entries in a single line delimited by spaces that defines the number of bins that we want to allocate to each of them.

It is mandatory that the training and test files be in the same format. So if DBNN_test.dat represents the test data, symbolic links with the same name to the inf, apf and awf train files is the simplest way to define them for the test data.

The two parameter files used to configure the runtime inputs to the classifier are to be labeled par0 and par1. These files are optional and has the command line inputs for the runtime settings of the classifier. If they are missing, the program will ask for the values of these parameters at run time. There is no need to edit par0 file and one can use the par0 file given in the distribution. The par1 file has two entries. The first line is the learning parameter that is typically set to 0.9 and the second line is the number of training epochs that defines how many times the classifier should try to update its weights before making a prediction. A typical value of 10 is usually appropriate.

Having created the training, test and parameter files, DBNN is executed by three simple commands executed in sequence.

autodbnn DBNN Out 0

autodbnn DBNN Out 1

autodbnn DBNN_test Out 2

Here autodbnn is the classifier program, 'DBNN' is the name of the data file without its extension .dat, 0 for the initial round, 1 the weight updating round and 2 the test round.

The output of the classifier is written to a .cmp file. 'Out' is a dummy label that allows different .cmp files to be created and stored from the same dataset for evaluations. In the above case, the command will produce the file DBNNOut.cmp and DBNN_testOut.cmp respectively for the training and test rounds. The .cmp file stores the first and second predictions along with the actual class labels and the respective confidence values the classifier has in each of its predictions. A few other intermediate files that are created and may be used for debugging are not discussed here.

 ## Photometric Classification of Quasars and other point sources in SDSS
 ![The animations show the 3D colour cube of about 6 million objects in our catalogue. The axis are the u-g, g-r and r-i colours derived from the SDSS u,g,r,i filters.](http://www.iucaa.in/~nspp/CC.gif) 
 ![The animation show the distribution of galaxies (black), stars (green) and quasars (blue) in sequence at different orientations of the cube.](http://www.iucaa.in/~nspp/CCsmooth.gif)

This page suppliments the "PHOTOMETRIC CLASSIFICATION OF QUASARS AND OTHER POINT SOURCES IN SDSS" by Abraham et.al, with DBNN software and the steps that were followed to produce the catalog of unresolved point sources. This content is licensed under the Creative Commons License.

A ready to run gzipped directory of the training and test data can be downloaded from here

The Revised Large Photometric Catalogue (July 2011 version) can be downloaded from here

The Our Predictions for Richard's Photometric Catalogue can be downloaded from here

### How to run the code and reproduce the results


You require Linux to run the software

Download the zipped directory and unzip it. It will create a directory named DBNN with a run subdirectory. Move to the run subdirectory and issue ./test from the Linux command prompt. It will first run autodbnn with option 0, that is the configuration mode, followed by option 1, the weights update mode and finally with option 2, the test mode. Option 0 and 1 are run on the training data 7.7K size while option 2 is run on the 90K test data - all the objects that have spectroscopic confirmation from the region. That is all.

The directory includes some script files that are used to format the output from DBNN to human readable form like confusion matrix

The subdirectory query has two sample queries, one to get training/test samples of spectroscopically confirmed object types and the other to get photometric data to prepare the catalog. We have used R-language to produce the data files required for DBNN. The R script with functions to do it are given in the DBNN directory. To try it, run R from the directory and press up arrow Key to see the commands. R is a very simple language to use and ten minutes on the R-manual that you may find on the web can help you run the commands. Note that the first command to run is the source(makedata0.R) that initiates the functions. Change the input file name in the makeU11Tdata with the name of your file of course! Some times, if the file is very large, R may just hang swapping memory. Use the makesub perl script to make your huge input file into small sub-files and then append the results files by setting apnd=TRUE in the subsequent R commands using each of the subsamples. DBNN can handle billions of rows and so combining the outputs is okay.

## Old Data and software

### DBNN for Star Galaxy Classification Problem

## How To

This page contains the DBNN software used in our paper on star galaxy classification. This software is licensed under the GNU Public License.

### ![Original DBNN Software - Click here](http://www.iucaa.in/~nspp/MySpace/MySpace/DBNN4SG.tar.gz).

### The ![Training and Test sets used by us for Qusar Classification are here](http://www.iucaa.in/~nspp/MySpace/MySpace/DBNN4SGData.tar.gz). Y

### ![The A&A paper may be found here](http://www.iucaa.in/~nspp/SG2.pdf).
