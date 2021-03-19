# Overview

This is a repository with the Python scripts to run the case study for selecting the pollution abatement activities for chemicals of concern and tracking chemical flows at the end-of-life stage. The data was obtained by means of data engineering using different publicly-available databases. The properties of chemicals were obtained using the GitHub repository [Properties_Scraper](https://github.com/jodhernandezbe/Properties_Scraper), while the PAU dataset using the repository [PAU4Chem](https://github.com/jodhernandezbe/PAU4Chem).

# Requirements:

This Pythos scripts were written using Python 3.x, Ubuntu 18.04, and Anaconda3. The following python packages are requiered to run the program:

1. [pandas](https://anaconda.org/conda-forge/pandas)
2. [graphviz](https://anaconda.org/conda-forge/graphviz)
3. [numpy](https://anaconda.org/conda-forge/numpy)
4. [scipy](https://anaconda.org/conda-forge/scipy)
5. [xlrd](https://anaconda.org/conda-forge/xlrd)
6. [xlutils](https://anaconda.org/anaconda/xlutils)
7. [python-graphviz](https://anaconda.org/conda-forge/python-graphviz)
8. [pomegranate](https://anaconda.org/anaconda/pomegranate)<sup>[1](#myfootnote1)</sup>

# Data-driven model: Bayesian Network (BN) for Pollution Abatement Unit (PAU) Identification

## Random variable names

| Name | Node| Name | Node |
| ------------- | ------------- | ------------- | ------------- |
| Is a byproduct? |	Node-1 | Waste flow	| Node-7 |
| Is a manufactured impurity?	| Node-2 | Chemical price	| Node-8 |
| Is a process impurity?	| Node-3 | Pollution abatement unit		| Node-9 |
| Type of waste	| Node-4 | CAPEX	| Node-10 |
| Chemical concentration	| Node-5 | OPEX	| Node-11 |
| Material efficiency	| Node-6 | Type of waste management	| Node-12 |

## Structure

<p align="center">
  <img src=https://github.com/jodhernandezbe/PAU_case_study/blob/master/bayesian_network/Bayesian_Network_PAU.png width="85%">
</p>

# Decision-making: Fuzzy Analytical Hierarchy Process (FAHP)

## Selection of PAU for a Concerning Chemical

<p align="center">
  <img src= https://github.com/jodhernandezbe/PAU_case_study/blob/master/fuzzy_analytical_hierarchy_process/FAHP_PAU.png width="100%">
</p>

## Sequence of PAUs for a Waste Stream<sup>[2](#myfootnote2)</sup>
<p align="center">
  <img src= https://github.com/jodhernandezbe/PAU_case_study/blob/master/fuzzy_analytical_hierarchy_process/FAHP_Seq.svg width="85%">
</p>

# Chemical Flow Tracking<sup>[3](#myfootnote3)</sup>

<p align="center">
 <img src= https://github.com/jodhernandezbe/PAU_case_study/blob/master/chemical_flow_analysis/pau_draws/PAU_sequence.svg width="95%">
</p>

# How to use

In order to use the code you need to dowload [TRI_releases.csv](https://drive.google.com/file/d/1sq3AzdCFMJ6Rh3Vx0dJau2YcMuu6T1nO/view?usp=sharing) and save the file in a folder which must be named tri_releases and located in the [chemical_flow_analysis](https://github.com/jodhernandezbe/PCU_case_study/tree/master/chemical_flow_analysis) folder. To run the case studies in the paper, navigate to the folder which contains main.py and write the following line on your terminal:

```
python main.py -CAS 67630 67561 7664417 107211 110543 108883 68122 75092 -Y 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 2004
```

The flag -CAS is for the CAS number for the chemicals used in the case study and -Y for the PAU dataset years used. The folder [Inputs.xls](https://github.com/jodhernandezbe/PAU_case_study/blob/master/Inputs.xls) has all the input information for the case study. The sheet ***Input*** contains the information for the streams, ***Chemical*** for the chemical substances, ***Industry_sector*** for the industry sectors, and ***Specifications*** for the guidelines about the factors in the Bayesian Network. 

When the following lines appear on your terminal, only press enter so that the program can continue running.

```
------------------------------------------------------------------------------------------------------------------------
Fill out the required information in the "Input.xls". Check the options in the sheet called "Specifications"

When you finish filling out the sheet "Input", please press enter
------------------------------------------------------------------------------------------------------------------------
```

# Output

The following are the program outputs:

1. In the folder [marginal](https://github.com/jodhernandezbe/PAU_case_study/tree/master/bayesian_network/probabilities/marginal), the files Marginal_probabilities_based_on_BN_for_CAS_in_stream_#.csv have the marginal proabilities under the BN. *CAS* indicates the chemical, while *#* the number for the stream that contains the chemical.
2. The file [PAU_selection_and_position_under_FAHP.csv](https://github.com/jodhernandezbe/PAU_case_study/tree/master/fuzzy_analytical_hierarchy_process/PAU_selection_and_position_under_FAHP.csv) has the PAU methods and sequences selected under the FAHP.
3. The file [Chemical_flow_tracking.csv](https://github.com/jodhernandezbe/PAU_case_study/blob/master/chemical_flow_analysis/Chemical_flow_tracking.csv) has the material flow analysis and allocation for the PAU methods and sequences. 
4. The folder [pau_draws](https://github.com/jodhernandezbe/PAU_case_study/tree/master/chemical_flow_analysis/pau_draws) has the draws that represent the flow allocation for the PAU methods and sequences.


# Disclaimer

The views expressed in this article are those of the authors and do not necessarily represent the views or policies of
the U.S. Environmental Protection Agency. Any mention of trade names, products, or services does not imply an endorsement by the U.S.
Government or the U.S. Environmental Protection Agency. The U.S. Environmental Protection Agency does not endorse any commercial products, service, or enterprises.

# Acknowledgement

This research was supported in by an appointment for Jose D. Hernandez-Betancur to the Research Participation
Program at the Center for Environmental Solutions and Emergency Response, Office of Research and Development,
U.S. Environmental Protection Agency, administered by the Oak Ridge Institute for Science and Education through an Interagency Agreement No. DW-89-92433001 between the U.S. Department of Energy and the U.S. Environmental Protection Agency.

-----------------------------------------------------------------------------------------------------------------------------

<a name="myfootnote1">1</a>: If you are using Anaconda distribution, you could have collisions between channels because pomegranate is only in the *anaconda* channel. Thus,if you have collision, activate your environment and run the following line:

```
conda config --set channel_priority false
```
<a name="myfootnote2">2</a>: NFPA: National Fire Protection Association

<a name="myfootnote3">3</a>: A03 is a code that represents Scrubber, while F71 Fume/Vapor incinerator. The list of codes for the PAU methos can be found in [Methods_TRI.csv](https://github.com/jodhernandezbe/PAU_case_study/blob/master/Methods_TRI.csv).




