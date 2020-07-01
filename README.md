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
  <img src=https://github.com/jodhernandezbe/PCU_case_study/blob/master/Bayesian_Network/Bayesian_Network_PCU.png width="85%">
</p>

# Fuzzy Analytical Hierarchy Process (FAHP)

## Selection of PCU for a Concerning Chemical

<p align="center">
  <img src= https://github.com/jodhernandezbe/PCU_case_study/blob/master/Fuzzy_Analytical_Hierarchy_Process/FAHP_PCU.png width="85%">
</p>

## Sequence of PCUs for a Waste Stream
<p align="center">
  <img src= https://github.com/jodhernandezbe/PCU_case_study/blob/master/Fuzzy_Analytical_Hierarchy_Process/FAHP_Seq.png width="85%">
</p>

-----------------------------------------------------------------------------------------------------------------------------
<a name="myfootnote1">1</a>: Properties_Scraper: https://github.com/jodhernandezbe/Properties_Scraper (Public).
<a name="myfootnote2">2</a>: PCUs: https://github.com/jodhernandezbe/PCUs (Private).

