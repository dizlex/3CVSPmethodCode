# PevznerCode
# Gives fracture orientation for zero-offset 3C VSP data at a specified depth interval

This is a MATLAB code for the (Pevzner et al 2011) method: Estimation of azimuthal anisotropy from 3C VSP data using multicomponent S-wave velocity analysis.

(Pevzner et al 2011) refers to this paper: https://espace.curtin.edu.au/handle/20.500.11937/16761

REQUIREMENTS:
1) The data organized using an SEG-Y format in a .txt file with or without headers must be put in the same directory as the downloaded git files.
2) The curveFittingToolBox must be installed in MATLAB 2015a and forward.
2) Powerpoint 2013 or later version (used for visualizing the results in a single file)

HOW TO RUN THIS PROGRAM:

0) Zero, download the repository into a folder of your computer containing your data.
1) First, run the a_ReadDataFromTextFile
2) Second, create a folder named after the geological formations or depth intervals to be analyzed. 
3) Third, copy all the files in the present working directory into the folder created in the previous step.
4) Fourth, run the pevznerMain and modify the inputs as needed.
5) Lastly, wait to see the results.

OUTPUTS OF THIS PROGRAM:

The images of the coherency plots at the specified depth intervals.

A powerpoint presentation containing all the specified depth intervals for result visualization.

A matrix named 'FinalCoh.mat' containing all the coherency results at every depth.

A .txt file containing the parameters used and the log of the errors presented in each run.

DESCRIPTION OF THE FILES:
The rest of the scripts which contain the pevzner name are functions which are called by pevznerMain, except for exportToPPTX which is also called.
The last files are used to:
1) z_DrawStratColumn
The z_DrawStratColumn draws an stratigraphic column if the formation well tops are available in an excel file.
The syms toolbox must be installed in MATLAB in order for z_DrawStratColumn.m file to run. 
2) z_clipping_sorting
The z_clipping_sorting file is used to do trace clipping.
the z_clipping_sorting maybe used for sorting at a single keyword of a file that contains the normal keyword of SeismicUNIX
The z_clipping_sorting is recommended only if SeismicUNIX is not available in the machine.
3) z_checingEachLevel
The z_checingEachLevel file is used for the final adjustments for plotting better looking graphs at specified intervals.

NOTE:
1) All the .m files are heavily commented aiming to be self-explanatory and to be modified according to the user needs.
2) MATLAB shows the error at every line, so you can correct as needed.
3) SeismicUNIX is hard to learn for the novices.

3.1) If there are doubts regarding the SEGY and SU formats please refer to the following book:
Seismic Data Processing with Seismic Un*x A 2D Seismic Data Processing Primer by David Forel, Thomas Benz, and Wayne D. Pennington.
Published by the SEG this is a great resource for novices and beginners learning SeismicUNIX and understanding SEG-Y files
