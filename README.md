# MMNetReduce

## Description

This project aim to reduce the number of reactions in a metabolic network by applying Michaelis-Menten reduction.
The program take as input a SBML file and return a SBML file if there's a Michaelis-Menten reduction applied.

## Content

1. [Installation](#installation)
2. [Usage](#usage)
3. [Project structure](#project-structure)
4. [Data](#data)
5. [Scripts](#scripts)
6. [Authors](#authors)
7. [References](#references)


## Installation

To install the program, you need to have python3 installed on your computer.

1. Clone the repository
```bash
git clone https://github.com/ClementRasp/MMNetReduce.git
```

or you can download the zip file and extract it.

2. Install the requirements
```bash
pip install -r requirements.txt
```


## Usage

To run the program, there's mutliple options:

1. Run the program with windows for selecting files in input and the directory for output.
```bash
python3 MMNetReduce.py
```

2. Run the program with options for input and output files.
- input_file_or_directory: specify the input file, can be a file or a directory (in this case the program gonna take all `.xml` files in the directory)
- -o: specify the output directory (Optional)
```
python3 MMNetReduce.py [-o output_directory] input_file_or_directory_1 ... input_file_or_directory_*
```

## Project structure

- `MMNetReduce.py` : Main file of the program from the paper
- `MMNetReduce_2.py` : A general methode to compute Michaelis-Menten Reductions
- `useful.py` : Methods for manipulating files
- `merge_reaction.py` : Methods for reductiong biomodel irreversible reactions into reversible
- `requirements.txt` : File containing the required libraries


## Data

The data used for the project are manually curated biomodel from [EBI](https://www.ebi.ac.uk/biomodels/).

## Scripts

- `MMNetReduce.py` : contains the code for applying Michaelis-Menten reduction following the method described in the paper.
- `MMNetReduce_2.py` : contains the code for applying Michaelis-Menten reduction following the general method of computing the intermediate of the reactions by using the rates of the reactions.
- `merge_reaction.py` : contains the code for merging the irreversible reactions into reversible reactions.
- `useful.py` : contains functions used by `MMNetReduce.py` and `MMNetReduce_2.py`.



## Authors

- Manvel Gasparyan - [manvelgasparyan](https://github.com/manvelgasparyan)
- Clément Raspail - [ClementRasp](https://github.com/ClementRasp)

## References

- Nisha Ann Viswan, Alexandre Tribut, Manvel Gasparyan, Ovidiu Radulescu, Upinder S Bhalla. Hierarchical Optimization of Biochemical Networks. 2024. [hal-04593669](https://hal.science/hal-04593669)