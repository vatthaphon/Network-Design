# Network Design

<p align="justify">Design a network of E and I LIF neurons.</p>

<p align="justify">This repository contains the source code for designing a network of E and I Leaky-Integrate-and-Fire neurons that generates desired periodic dynamics.</p>

## Table of Contents  
[1. Introduction](#Introduction)  
[2. Installation](#Installation)  
[3. Usage](#Usage)  
[4. File Structure](#FileStructure)  
[5. License](#License)  
          
## Introduction<a name="Introduction"/>
<p align="justify">This repository provides the implementation of the methods and algorithms to design a network. 
The main goal of this study is to design a network of pulse-coupled delayed E and I LIF neurons with desired periodic spiking.</p>

## Installation<a name="Installation"/>
Matlab 2018.

## Usage<a name="Usage"/>
Run Network_analysis/networkAnalysis.m to create a network with predefined spiking times.

## File Structure<a name="FileStructure"/>
```plaintext
├── Network_analysis/            # Design the network.
├── Get_weights/                 # Generate synaptic weights between the LIF neurons using quadprog of Matlab.  
├── Gen_delays/                  # Generate synaptic delays between the LIF neurons.
├── Conn_matrix/                 # Generate a matrix specifying if two neurons are connected.
├── Extract_pdf_from_the_paper/  # Generate spiking time of the E and I LIF neurons.
├── Sim_network/                 # Generate dynamics the network that we designed.  
├── README.md  
└── LICENSE
```

## License<a name="License"/>
This project is licensed under the MIT License - see the LICENSE file for details.
