# TwisTranscripT: Simulation of the Transcription-Supercoiling Coupling in Bacteria

This package aims at simulating the Transcription-Supercoiling Coupling (TSC) and its impact in gene expression in Bacteria. We used the Python programming language with the help of a Python-based scientific ecosystem for scientific computing. Equations, calibration and application of the model are described in [El Houdaigui et al., NAR 2019](https://doi.org/10.1093/nar/gkz300). 
Topoisomerase parameters may be adapted for eukaryotes, but not yet tested. 

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

The package is based on Python3. 

The following libraries are required in order for the script to run successfully:
* [Numpy](http://www.numpy.org/) - Fundamental package for scientific computing with Python.
* [Pandas](https://pandas.pydata.org/) - A library providing high-performance, easy-to-use data structures and data analysis tools.

The following libraries are required for the optional plotting scripts: 
* [matplotlib](https://matplotlib.org/) - A plotting library for the Python programming language and its numerical mathematics extension NumPy.
* [dnaplotlib](https://github.com/VoigtLab/dnaplotlib) - A library that enables highly customizable visualization of individual genetic constructs and libraries of design variants (only for visualization).

### Installation

The installation is simple in Linux/Ubuntu using `pip`. Make sure that Python3 is the default installed version of Python. 

1. Download the package & open the terminal in the script main directory.
2. Install the package using `pip`

```
pip install .
```

In MacOS, Windows or Linux, you can also download the package, include your installation directory into your PYTHONPATH
```
PYTHONPATH=$PYTHONPATH:\PATH_TO_INSTALLATION_DIR
```
Then either call the main function start_transcribing within your Python scripts
```
from TSC import start_transcribing
start_transcribing(...)
```
or use the executable script ```./start_simulation.py INI_file OUTPUT_dir``` in the main installation directory (make sure the script can be executed by the user).


## Inputs, Outputs & Parameters files :

### Input files :

#### The General Feature Format file (GFF):
Classical genome description file; Each line of the [GFF](https://www.ensembl.org/info/website/upload/gff.html) file contains 9 columns of data : *Seqname*, *Source*, *Feature*, *Start*, *End*, *Score*, *Strand*, *Frame* and *Attribute*.

This file is mainly used to extract the simulated genome length, given by the "region" feature. Gene positions are not used by the main script. 

```
##gff-version 3
#!gff-spec-version 1.20
#!processor NCBI annotwriter
##sequence-region chrom1genes 1 34000
chrom1genes RefSeq  region  1 34000 . + . ID=id0;Name=chrom1genes
chrom1genes RefSeq  gene  1201  13200 . + . ID=gene0;Name=gene0
```

#### Transcription Start Site file (TSS):
The TSS file describes the list of TSS in the simulated genome. For each of them, it contains 4 information data :
* *TUindex* : The transcription unit (TU) index.
* *TUorient* : The TU orientation ("+" for forward or "-" for reverse)
* *TSS_pos* : The TSS position in the chromosome.
* *TSS_strength* : The basal initiation rate (in $s^{-1}$).

The lists of TSSs and TTSs (below) together define the transcriptional landscape of the provided genome. In principle, the program first determines all possible transcripts with these start and stop sites  (for each TSS, transcription ends when a perfect terminator is reached in the given orientation). 
The use of TUs is optional: if all TSS and TTS (see below) have the same TU index, the program will generate a list of possible transcripts based on their positions and orientations in the genome. If different TU indexes are provided, then transcription can only occur between TSSs and TTSs of the same TU (take care that a perfect terminator is given for each TU, with meaningful position and orientation). In all cases, a list of possible transcripts with detailed information is generated (see outputs).

```
TUindex TUorient  TSS_pos TSS_strength
0 - 4250  0.001
1 + 5150  0.015
2 + 16150 0.02
3 - 22250 0.02
```

#### Transcription Termination Site file (TTS):
The TTS file also provides the *TUindex*, *TUorient*, *TTS_pos* and :
* *TTS_proba_off* : termination probability. A perfect terminator has a value of 1. For a value <1, stochastic readthrough can occur with a probability (1-p).
```
TUindex TUorient  TTS_pos TTS_proba_off
0 - 3150  1.
1 + 6250  .6
2 + 17250 .4
3 - 21150 1.
```

#### Fixed topological barrier file (BARR_FIX):
This file (default prot.dat) contains a list of fixed topological barrier positions, which typically represents topological insulator proteins such as H-NS. 

```
prot_name	prot_pos
hns 		2000
hns		25000
```

### Output files :

Structure of output files :

```
output
├── all_tr_info.csv
├── save_tr_times.csv
├── save_nbr_RNAPs_hooked.npz
├── all_res
│   ├── save_RNAPs_info.npz
│   ├── save_sigma_info.npz
│   └── save_tr_info.npz
└── resume_sim
    ├── resume_sim_Barr.npz
    ├── resume_sim_RNAPs.npz
    └── resume_sim_tr.npz


```
Main output files:
  * *all_tr_info.csv*: (Human readable csv file) List of possible transcripts in the simulation (based on the provided TSS/TTS files), with detailed information: TU index (as provided in TSS/TTS files), transcript ID, TSS (start) position, TSS strength (basal rate), strand (+1 or -1), end position, transcript strength (computed from the TSS basal rate, corrected by TTS probabilities in case of imperfect terminators), number of expressed RNAs for each transcript.
  * *save_tr_times.csv*: (Human readable csv file) RNAs expression times for each transcript ID (in seconds). 
  * *save_nbr_RNAPs_hooked.npz* : Number of elongating RNA Polymerases at regular time intervals along the simulation. 

Additional output files: 
  * *all_res* : contains files in which information is saved **at constant time intervals during the whole simulation**:
      * *save_RNAPs_info.npz*: Elongating RNA Polymerase IDs and their positions along the genome.
      * *save_sigma_info.npz*: Information related to supercoiling: size of each region (regions are delimited by protein barriers and bound RNAPs), supercoiling density (<a href="https://www.codecogs.com/eqnedit.php?latex=\sigma" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\sigma" title="\sigma" /></a>) in each of them, and genomic average value. 
      * *save_tr_info.npz*: Information about transcripts: number of each possible transcript (as described in the info file below) and instantaneous initiation rate (<a href="https://www.codecogs.com/eqnedit.php?latex=k_{on}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?k_{on}" title="k_{on}" /></a>, i.e. basal rate modulated by instantaneous supercoiling level at the promoter). 
  * *resume_sim*: Files saved **at the end of the simulation**, in order to resume it later (possibly with different parameters, e.g. to simulate gyrase inhibition).
      * *resume_sim_Barr.npz*: Information related to barriers: barrier positions and types, regions size and supercoiling density in each of them, and remaining timesteps of each elongating RNA Polymerase before reaching its termination site. 
      * *resume_sim_RNAPs.npz*: Elongating RNA Polymerases IDs and positions along the genome, unbound RNA Polymerases IDs. 
      * *resume_sim_tr.npz* :  Number of transcripts and TSS instantaneous initiation rates (<a href="https://www.codecogs.com/eqnedit.php?latex=k_{on}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?k_{on}" title="k_{on}" /></a>)

 NOTE : `.npz` files can be opened only by using [`numpy.load`](https://docs.scipy.org/doc/numpy-1.13.0/reference/generated/numpy.load.html) function. The different stored items are then available using dictionary keys. 

### Parameter file :

This configuration file of standard format (https://docs.python.org/3/library/configparser.html) contains simulation parameters separated into 5 sections. **In most cases, the user should edit only the first and last sections**:

1. **INPUTS: paths to the *GFF*, *TSS*, *TTS* and *BARR_FIX* files**

2. PROMOTER: parameters of the promoter response curve to local SC variations, based on a thermodynamic model of promoter opening: supercoiling activation factor, crossover threshold and width (see [El Houdaigui et al., NAR 2019](https://doi.org/10.1093/nar/gkz300)). 
3. TOPOISOMERASES: activity parameters for the two main topoisomerases. Default values are calibrated for bacteria (topoisomerase 1 and DNA gyrase), and may be adapted to eukaryotes: activity constant (in $s^{-1}$), crossover threshold and width. 
4. GLOBAL: space and time discretisation units. Default values: 60~nt and 2~s.

5. **SIMULATION: main simulation parameters**
    * *RNAPs_genSC*: Supercoiling generation rate by elongating RNAPs (default 0.2 coils per unit length)
    * *SIGMA_0*: Initial uniform supercoiling density (default -0.06)
    * *RNAPS_NB*: The number of RNA Polymerases (default 3)
    * *SIM_TIME*: Simulation time (in seconds).
    * *OUTPUT_STEP*: Time interval for information output (in seconds).
    * *GYRASE_CONC*: Gyrase concentration (in micromoles/L)
    * *TOPO_CONC*: Topoisomerase I concentration (micromoles/L)

Example of parameter file:
```
[INPUTS]
gff = test.gff
tss = TSS.dat
tts = TTS.dat
barr_fix = prot.dat

[PROMOTER]
m = 2.20
sigma_t = -0.042
epsilon = 0.005

[GLOBAL]
delta_x = 60
delta_t = 2

[SIMULATION]
rnaps_gensc = 0.2
sigma_0 = -0.05
rnaps_nb = 10
output_step = 100
gyrase_conc = 0.1
topo_conc = 0.05
sim_time = 10000

[TOPOISOMERASES]
gyrase_cte = 0.01
topo_cte = 0.005
k_gyrase = 50
x0_gyrase = 0.016
k_topo = 80
x0_topo = -0.04

```



## Code Example

You can test the example provided by typing the following command from the main directory :

```
# Execute the script by providing the parameter file and the input directory (the output directory path is optional)
# python start_simulation.py path/to/the/params.ini [output/path]
python start_simulation.py example_1/params.ini
```

NOTE : The input directory should contain the *GFF*, *TSS*, *TTS*  and the *Protein barrier* file.

If the simulation ended successfully, you'll get the output files as described [above](#output-files).

The plotting package TSC_plotting.py contains some useful functions for result analysis and plotting, but the latter usually require a specific code. 

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.

#### Contact
* [sam.meyer@insa-lyon.fr](sam.meyer@insa-lyon.fr)
* [bilal.elhoudaigui@gmail.com](bilal.elhoudaigui@gmail.com).
