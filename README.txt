This program was created for my Master's thesis "Simulation of the one-locus K-type ancestral selection graph with parent-independent mutation".

It simulates ancestral selection graphs and aggregates the results, which will be output on the commandline.

If you have this readme, then you should also have a copy of my Master's thesis if any questions arise.

Usage: simulate_asgs.py [-h] -input REPEATS -reproductive_advantages [S ...] -mutation_rate T -time R -mutation_probabilities [V ...] [-dirichlet]

Options:
  -h, --help            show this help message and exit
  -input, -i REPEATS    The number of graphs simulated.
  -reproductive_advantages, -sigmas, -s [S ...]
                        A list of the reproductive advantages of the different types in the diffusionlimit, starting with type 0 and ending with type K-2. s_i needs to be higher than s_i+1. K-1 will be automatically added as    
                        s_(K-1)=0. Simply write the numbers separated by spaces e.g: "-s 0.3 0.2 0.1"
  -mutation_rate, -theta, -t T
                        The rate with which an individual mutates in the diffusionlimit.
  -time, -R R           The length of time graph will be simulated for.
  -mutation_probabilities, -nus, -v [V ...]
                        A list of the probabilities with which an individual mutates into a specific type, if a mutation occurs, starting with type 0 and ending with type K-1. The sum of these values needs to be 1. Simply       
                        write the numbers separated by spaces e.g: "-v 0.2 0.3 0.4 0.1"
  -dirichlet, -d        If this command is added a Dirichlet-distribution is used to determine potential ancestors instead of a Wright-distribution. This may be useful for certain parameter configurations, as the rejection      
                        sampling process used for the Wright-distribution may become extremely unlikely to produce results. The type distribution on the ancestral line should still provide a very close approximation of the      
                        ancestral type distribution.


Example: to simulate 1000 graphs of length 100 with reproductive advantages sigma=(5,4,3,2,1,0), mutation rate theta=6 and mutation probabilities nu=(0.01, 0.04, 0.05, 0.1, 0.3, 0.5) use the command:
simulate_asgs.py -i 1000 -R 100 -s 5 4 3 2 1 -t 6 -v 0.01 0.04 0.05 0.1 0.3 0.5 
If i wanted to use a dirichlet distribution to determine the ancestor's types, I would append the "-dirichlet" parameter.

How to set up Python:::::::::::::::::::::::::::::::::::::::

Python should come with a linux installation.
If you are on Mac or Linux, download the newest version from https://www.python.org/downloads/. 

This program was written with Python 3.13 in case future version don't work, try switching to it.

ON WINDOWS:::::::::::::::::::::::::::::::::::::::::::::

Next, i heavily recommend setting up a virtual environment.
To do that, navigate, via the command line, to the directory the program is located in and type "python -m venv ." into the command line.
Next, download the required modules into this new environment by typing "python -m pip install numpy" and "python -m pip install scipy" into the command line.
Now you can use "python simulate_asgs.py [parameters]" to start the program.

ON LINUX:::::::::::::::::::::::::::::::::::::::::::

Setting up the virtual environment:
On Linux you might need to install venv before using it. Navigate to the folder you want to work in (with the simulat_asgs.py) and type in "sudo apt install python3.12-venv". You may need to run "sudo apt update" beforehand.
Then proceed with "python3 -m venv .".
pip might also need to be installed with "sudo apt install python3-pip"
Then use pip to install numpy and scipy with "python3 -m pip install numpy --break-system-packages" and "python3 -m pip install scipy --break-system-packages"
Now you can use "python3 simulate_asgs.py [parameters]" to start the program.
