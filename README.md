# Overview

This is a repository with the Python scripts to run the case study for selecting the pollution abatement activities for concerning chemicals and tracking chemical flows at the end-of-life stage. The data was obtained be means of data engineering using different publicly-available databases. The properties of chemicals were obained using the public GitHub repository "Properties_Scraper" <sup>[1](#myfootnote1)</sup>, while the PCU database using the private repository "PCUs" <sup>[2](#myfootnote2)</sup>.

# Bayesian Network (BN) for Pollution Control Unit (PCU) Selection

### Factor names

| Name | Node| Name | Node |
| ------------- | ------------- | ------------- | ------------- |
| Byproduct |	Node-1 | Waste flow	| Node-7 |
| Manufactured impurity	| Node-2 | Chemical price	| Node-8 |
| Process impurity	| Node-3 | Pollution control unit		| Node-9 |
| Type of waste	| Node-4 | Pollution abatement capital expenditure	| Node-10 |
| Concentration	| Node-5 | Pollution abatement operating cost	| Node-11 |
| Efficiency	| Node-6 | Type of waste management	| Node-12 |

### Structure

<p align="center">
  <img src=https://github.com/jodhernandezbe/PCU_case_study/blob/master/bayesian_network/bayesian_network_PCU.png width="85%">
</p>

# Fuzzy Analytical Hierarchy Process (FAHP)

## Selection of PCU for a Concerning Chemical

<p align="center">
  <img src= https://github.com/jodhernandezbe/PCU_case_study/blob/master/fuzzy_analytical_hierarchy_process/FAHP_PCU.png width="85%">
</p>

## Sequence of PCUs for a Waste Stream
<p align="center">
  <img src= https://github.com/jodhernandezbe/PCU_case_study/blob/master/fuzzy_analytical_hierarchy_process/FAHP_Seq.png width="85%">
</p>

## Disclaimer

The views expressed in this article are those of the authors and do not necessarily represent the views or policies of
the U.S. Environmental Protection Agency. Any mention of trade names, products, or services does not imply an endorsement by the U.S.
Government or the U.S. Environmental Protection Agency. The U.S. Environmental Protection Agency does not endorse any commercial products, service, or enterprises.

## Acknowledgement

This research was supported in by an appointment for Jose D. Hernandez-Betancur to the Research Participation
Program at the Center for Environmental Solutions and Emergency Response, Office of Research and Development,
U.S. Environmental Protection Agency, administered by the Oak Ridge Institute for Science and Education through an Interagency Agreement No. DW-89-92433001 between the U.S. Department of Energy and the U.S. Environmental Protection Agency.

-----------------------------------------------------------------------------------------------------------------------------
<a name="myfootnote1">1</a>: Properties_Scraper: https://github.com/jodhernandezbe/Properties_Scraper (Public).

<a name="myfootnote2">2</a>: PCUs: https://github.com/jodhernandezbe/PCUs (Private).

