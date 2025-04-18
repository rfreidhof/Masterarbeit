This program was created for my Master's thesis "Simulation of the one-locus K-type ancestral selection graph with parent-independent mutation".

It simulates ancestral selection graphs and aggregates the results, which will be output on the commandline.

If you have this readme, then you should also have a copy of my Master's thesis if any questions arise.

Before executing the program, you may need to install numpy and scipy, which are 3rd party modules.
The pip functionality should come with your python installation.
Simply type "python -m pip install scipy" and "python -m pip install numpy" into your command line before using the program.

Usage: simulate_asgs.py [-h] [-input REPEATS] -reproduction_rates [S ...] -mutation_rates T -time R -mutation_probabilities [V ...] [-dirichlet]

options:
  -h, --help            show this help message and exit
  -input, -i REPEATS    The number of graphs simulated. Default = 1000
  -reproduction_rates, -sigmas, -s [S ...]
                        A list of the reproductive advantages of the different types in the diffusionlimit, starting with type 0 and ending with type K-2. s_i needs to be higher than s_i+1. K-1 will be automatically added as    
                        s_(K-1)=0. Simply write the numbers separated by spaces e.g: "-s 0.3 0.2 0.1"
  -mutation_rates, -theta, -t T
                        The rate with which an individual mutates in the diffusionlimit.
  -time, -R R           The length of time graph will be simulated for.
  -mutation_probabilities, -nus, -v [V ...]
                        A list of the probabilities with which an individual mutates into a specific type, if a mutation occurs, starting with type 0 and ending with type K-1. The sum of these values needs to be 1. Simply       
                        write the numbers separated by spaces e.g: "-v 0.2 0.3 0.4 0.1"
  -dirichlet, -d        If this command is added a Dirichlet-distribution is used to determine potential ancestors instead of a Wright-distribution. This may be useful for certain parameter configurations, as the rejection      
                        sampling process used for the Wright-distribution may become extremely unlikely to produce results. The type distribution on the ancestral line should still provide a very close approximation of the      
                        ancestral type distribution.