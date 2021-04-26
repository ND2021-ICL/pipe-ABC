#######

This folder holds all files relevant to the ABC portion of the ND2021 project
All files used Python/ 3.8.5 and packagaes included in Anaconda
Files must be within the same folder as the data files to be able to run

It includes:

###

4 x python (.py) scripts for the four different GRN models that were parameterised using ABC:
- SMCABC_thermo_8ODEs.py: Thermodynamic mRNA regulatory functions with mRNA and protein incorporated into one equation (4 ODEs) - Cohen et al., 2014
- SMCABC_thermo_4ODEs.py: Thermodynamic mRNA regulatory functions with mRNA and protein in seprate equations (8 ODEs) - Rayon et al., 2020
- SMCABC_hill_8ODEs.py: Hill function mRNA regulatory functions with mRNA and protein incorporated into one equation (4 ODEs) - new model proposed by us
- SMCABC_hill_4ODEs.py: Hill function mRNA regulatory functions with mRNA and protein in seprate equations (8 ODEs) - new model proposed by us

The .py scripts use multiprocessing and can be used to fit the models to either the mouse or human transcript data. This is specified in the bash script

Multiprocessing doesn't work well in Jupyter notebooks so keep in mind if using part of this code. These scripts were all run using 8 cores on Imperial's High Performance Computer farms, to fit each model to both the mouse and human variance stabilied transformed expression data obtained from earlier parts of the project, resulting in 8 ABC runs being carried out. The main shell line to run the script would be as follows:

python SMCABC_DESeq_ruben_4ODEs_maxscore.py $NUMBER_OF_THREADS $True >> $OUTPUT_FILE 2>&1

- Set True to use mouse data, False to use human data
- Running the script will also save the final layer parameters to a csv file, and produce a plot of the trajectories against the experimental data for each layer, to be able to track the progress of the runs
- The python file to be run needs to be changed within the shell script (DEseq.com)

DEseq.com: The bash script used to run the python scripts on Imperial College London's server. Run using the folloing line:

qsub -q long -lnodes=1:ppn=8:msc3 DEseq.com true

- ppn = specifies number of nodes to use - set in command line as depends on what hardware available
- msc3 = specifies node to use
- DEseq.com = shell script to run
- true/ false = whether to run on mouse (true) or human (false) data

##

Results processing.ipynb: The Jupyter notebook used to analyse the results of the ABC parameterisation for EACH of the models above fit to mouse and human data, including:
- Visualising the trajectories created using the 1000 parameter particles accepted in the final layer of the SMC ABC
- Creating joint pairwise ditributions of the postarior paramters in a gridplot
- Calculating, visualising and scoring the models created using the parameter values at which its distribution peaked
- Calculating the credibility intervals of the posterior parameters
- Comparing the posterior distributions obtained for protein production and degradation in each model

This file takes input in from the following results csvs, which include the paarameter values of the final 1000 particles obtained from each ABC run in a csv (these were renamed from the csvs obtained from the .py script). The name includes the final threshold the ABC managed to reach within 3 days of running: 
- 'MOUSE4_riju_2021_03_11_1426_final_params_thresh0.3.csv'
- 'HUMAN4_riju_2021_03_11_1540_final_params_thresh0.2.csv'
- 'HUMAN4_ruben_2021_03_11_1459_final_params_thresh0.2.csv'
- 'MOUSE4_ruben_2021_03_11_1406_final_params_thresh0.3.csv'
- 'MOUSE8_riju_2021_03_13_0431_final_params_thresh0.4974474934820172.csv'
- 'HUMAN8_riju_2021_03_13_1710_final_params_thresh0.3.csv'
- 'MOUSE8_ruben_2021_03_13_1343_final_params_thresh0.45067636708696424.csv'
- 'HUMAN8_ruben_2021_03_13_1207_final_params_thresh0.35.csv'

NOTE: Where 'riju' is used, this refers to Hill function models, where 'ruben' is used, this refers to thermodynamic models (due to who the different models were provided by)

###

Simulating mutant patterning.ipynb: A Jupyter notebook used to model the spatial patterning as done by Cohen et al. (2014), as well as with the parameters obtained from parameterisation with the SMCABC scripts above.

###

Data files needed (not included)
- vst_transformed.csv: contains the experimental mouse and human reads
- The 8 csvs of results listed above 

#######

Email maiseydavidson99@gmail.com for any questions/ clarifications about the code!
