This program was created for my Master's thesis "Simulation of the one-locus K-type ancestral selection graph with parent-independent mutation".

It simulates ancestral selection graphs and aggregates the results, which will be output on the commandline.

If you have this readme, then you should also have a copy of my Master's thesis if any questions arise.

For people not familiar with Python:
Before executing the program, you may need to install numpy and scipy, which are 3rd party modules.
I recommend creating a virtual environment before installing them. To do that, navigate (through the command line) to the directory in which simulate_asgs.py is located and type "python -m venv".
Then, to install the modules, type "python -m pip install scipy" and "python -m pip install numpy" into the command line before using the program.

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
python simulate_asgs.py -i 1000 -R 100 -s 5 4 3 2 1 -t 6 -v 0.01 0.04 0.05 0.1 0.3 0.5 
