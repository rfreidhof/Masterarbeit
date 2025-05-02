"""
This program is a simulation of the ancestral selection graph with a single mutation locus and any amount of types with different reproductive advantages.
It was created for my master's thesis "Simulation of the one-locus K-type ancestral selection graph with parent-independent mutation" at Universität Bielefeld.
It, the file used to create the graphics used within and the licenses for modules used in this software can be found at https://github.com/rfreidhof/Masterarbeit. - Richard C. Freidhof

This software is licensed under the following license:

Copyright (c) 2025 Richard Casparus Freidhof.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above
   copyright notice, this list of conditions and the following
   disclaimer in the documentation and/or other materials provided
   with the distribution.

3. Neither commercial redistribution nor commercial use of this software, or any
   software derived from it, are permitted without specific
   prior written permission by the copyright holder.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""




import sys
import math
import numpy.random as random
from argparse import ArgumentParser
import copy

rng = random.default_rng()

#Helpfunctions------------------------------------------------------------------------

def add_L_r(L_r_list, L_r):
    while len(L_r_list) <= L_r:
        L_r_list.append(0)
    L_r_list[L_r] += 1

def sum_lists_of_different_lengths(list1, list2):
    while len(list1) < len(list2):
        list1.append(0)
    while len(list2) < len(list1):
        list2.append(0)
    return [sum(x) for x in zip(list1, list2)]


#Calculations -------------------------------------------------------------
from scipy.integrate import quad

def normalizing_constant(sigma0, theta, v):
    C = quad(lambda x, sigma0, theta, v : x**(theta*v[0]-1)*(1-x)**(theta*v[1]-1)*math.exp(sigma0*x), 0, 1, args=(sigma0, theta, v))[0]
    return C

def wright_expected_value(sigma0, theta, v):
    y = quad(lambda x, sigma0, theta, v : x**(theta*v[0])*(1-x)**(theta*v[1]-1)*math.exp(sigma0*x) , 0, 1, args=(sigma0, theta, v))[0]
    C = normalizing_constant(sigma0, theta, v)
    return y/C


def h_integral(sigma0, theta, v, n):
    y = quad(lambda x, sigma0, theta, v, n : x**(theta*v[0])*(1-x)**(theta*v[1]+n-1)*math.exp(sigma0*x), 0, 1, args=(sigma0, theta, v, n))[0]
    C = normalizing_constant(sigma0, theta, v)
    return y/C


def b(n, theta, sigma0, v):
    y = quad(lambda x, sigma0, theta, v, n : x**(theta*v[0]-1)*(1-x)**(theta*v[1]+n-1)*math.exp(sigma0*x), 0, 1, args=(sigma0, theta, v, n))[0]
    C = normalizing_constant(sigma0, theta, v)
    return y/C

def increment_a(n, theta, sigma0, nu1, a_n_minus_1, a_n_minus_2): #A non-recursive method, greatly reducing computation time.
    return ((n+sigma0+theta)*a_n_minus_1-sigma0*a_n_minus_2)/(n+theta*nu1)

def a(n, theta, sigma0, v):
    a_n_minus_2 = 1
    a_n_minus_1 = (sigma0/(1+theta*v[1]))*((b(2, theta, sigma0, v)-b(3, theta, sigma0, v))/(b(1, theta, sigma0, v)- b(2, theta, sigma0, v)))
    if n == 0:
        return a_n_minus_2
    elif n== 1:
        return a_n_minus_1
    else:
        for i in range(2, n+1):
            a_n = increment_a(i, theta, sigma0, v[1], a_n_minus_1, a_n_minus_2)
            a_n_minus_2 = a_n_minus_1
            a_n_minus_1 = a_n
        return a_n

def altbigp(theta, sigma0, v):
    a_n_minus_2 = 1
    a_n_minus_1 = (sigma0/(1+theta*v[1]))*((b(2, theta, sigma0, v)-b(3, theta, sigma0, v))/(b(1, theta, sigma0, v)- b(2, theta, sigma0, v)))
    sum = wright_expected_value(sigma0, theta, v)+a_n_minus_1*h_integral(sigma0, theta, v, 1) #The values for n=0 and n=1
    n = 2
    while True:
        a_n = increment_a(n, theta, sigma0, v[1], a_n_minus_1, a_n_minus_2)
        new_v = a_n*h_integral(sigma0, theta, v, n)
        a_n_minus_2 = a_n_minus_1
        a_n_minus_1 = a_n
        n += 1
        sum += new_v
        if new_v < 0.001:
            return sum

def L_inf_dist(n, sigma0):
    return (sigma0**n)/(math.factorial(n)*(math.exp(sigma0)-1))

#For mutation rate:

def qA10 (theta, sigma0, v):
    sum1 = 0
    sum2 = 0
    n = 1
    while True:
        new_sum1 = a(n-1, theta, sigma0, v)*b(n, theta, sigma0, v)
        new_sum2 = (a(n-1, theta, sigma0, v)-a(n, theta, sigma0, v))*b(n, theta, sigma0, v)
        sum1 += new_sum1
        sum2 += new_sum2
        n += 1
        if new_sum1 < 0.0000001 and new_sum2 < 0.0000001:
            return theta*v[0]*(sum1/sum2)

def qA01 (theta, sigma0, v):
    sum1 = 0
    sum2 = 0
    n = 1
    while True:
        new_sum1 = (a(n-1, theta, sigma0, v)-a(n, theta, sigma0, v))*(b(n-1, theta, sigma0, v)-b(n, theta, sigma0, v))
        new_sum2 = a(n-1, theta, sigma0, v)*(b(n-1, theta, sigma0, v)-b(n, theta, sigma0, v))
        sum1 += new_sum1
        sum2 += new_sum2
        n += 1
        if new_sum1 < 0.0000001 and new_sum2 < 0.0000001:
            return theta*v[1]*(sum1/sum2)




#Classes----------------------------------------------------------------------------------------------------------
class Individual:
    def __init__(self, typ, ancestor_type):
        self.typ = typ
        self.ancestor_type = ancestor_type

class Event:
    def __init__(self, r, target):
        self.r = r #time from the present in backwards time.
        self.target = target

class BranchingEvent(Event):
    def __init__(self, r, origin, target, typ):
        super().__init__(r, target)
        self.origin = origin
        self.typ = typ
        self.is_silent = True #Knowing whether a branching event is silent or not means one does not need to employ a backtracking algorithm in reverse time.

    def apply_event_forward(self, individuals): #Applies the event in forward time.
        if individuals[self.origin].typ <= self.typ:
            individuals[self.target] = copy.deepcopy(individuals[self.origin])
            self.is_silent = False
        del individuals[self.origin] #Either way, the origin individual is not considered going forward.

class CoalescingEvent(Event):
    def __init__(self, r, origin, target):
        super().__init__(r, target)
        self.origin = origin
        self.is_silent = False #Coalescing events are never silent.

    def apply_event_forward(self, individuals): #Applies the event in forward time.
        individuals.insert(self.target, Individual(None, None))#Creates an empty individual at the target location, so the target and origin index line up at the time of the event
        individuals[self.target] = copy.deepcopy(individuals[self.origin])

class MutationEvent(Event):
    def __init__(self, r, target, typ):
        super().__init__(r, target)
        self.origin = target
        self.typ = typ
        self.previousType = None
        self.is_silent = True

    def apply_event_forward(self, individuals): #Applies the event in forward time. This should always happen before applying an event in backwards time.
        self.previousType = individuals[self.target].typ # Save the target's previous type for use in backwards time.
        if self.previousType != self.typ:
            individuals[self.target].typ = self.typ
            self.is_silent = False

#Functions implementing algorithms ---------------------------------------------------------------------------------------------

def generate_event(r, L_r, s, v, E):
    e = rng.random()
    target = rng.integers(low=0, high=L_r) #low inclusive, high exclusive. Thus this can be used with the "individuals" list
    if e <= (L_r * s[0])/E: #Branching
        p = rng.random()
        for i in range(len(s)-1): #Last entry is of s_(K-1), which is the coalescing event
            if p >= (s[i+1])/s[0]:
                return 1, BranchingEvent(r=r, origin=L_r, target=target, typ=i)

    elif e < ((L_r * s[0])+(L_r*(L_r-1)))/E: #coalescence
        origin = rng.integers(low=0, high=L_r)
        while target == origin:
            origin = rng.integers(low=0, high=L_r) #The possibility of selecting the same individual twice is already accounted for in the probability, so we reroll until they're not the same.
        return -1, CoalescingEvent(r=r, origin = origin, target = target)

    else: #mutation
        m = rng.random()
        for i in range(len(v)):
            if m <= sum(v[0:i+1]):
                return 0, MutationEvent(r=r, target=target, typ=i)


def generate_graph(s, th, R, v): #Generates an ancestral selection graph backwards in time
    L_r = 1
    E = L_r * ((L_r-1)+th+s[0])#Rate at which events happen. Needs to be recalculated each cycle due to L_r changing
    tau = rng.exponential(scale=1/E, size=None) #Scale is the inverse of the rate, or the expected average.
    r = tau
    eventlist = []

    while r < R:
        addition, event = generate_event(r, L_r, s, v, E)
        eventlist.append(event)

        if addition != 0: #An if block is roughly 1% faster than recalculating every cycle.
            L_r += addition
            E = L_r * ((L_r-1)+th+s[0])

        tau = rng.exponential(scale=1/E, size=None)
        r += tau
    #Could return R-(r-tau) to get time from last event to R, but that wasn't needed in the program
    return L_r, eventlist #Returns the time from R (or t=0) to the last (first) event, number of individuals and the graph in backwards time.

def generate_potential_ancestors(L_R, theta, vs, sigmas, dirichlet):
    ancestors = []
    a = [] #positive concentration parameters
    for i in range(len(vs)):
        a.append(theta*vs[i])

    if dirichlet:
        V = rng.dirichlet(a) # Dirichlet sample
    else:
        counter = 0
        while True:
            V = rng.dirichlet(a) # Dirichlet sample
            U = rng.random()
            max = -sigmas[0]
            for i in range(len(sigmas)):
                max += sigmas[i]*V[i]
            counter += 1
            if counter == 1000000:
                print("More than 1 million attempts to generate a suitable distribution of ancestral types have been made without success.\n It may be the case that the chosen parameters are particularly ill suited to the rejection sampling approach.\n You can rerun the program with the -dirichlet command to instead use a dirichlet distribution.")
            if U <= math.exp(max):
                break
    for i in range(L_R):
        ancestor_type_chance = rng.random()
        for j in range(len(V)):
            if ancestor_type_chance <= sum(V[0:j+1]):
                ancestors.append(Individual(typ=j, ancestor_type =j))
                break


    return ancestors


def type_graph(individuals, eventlist, K, R):
    asgen_type_distribution = [0]*K
    last_r = R
    L_R = len(individuals)

    last_mut_rs = [R]*L_R #Times since the specific bloodline had a mut event, starts at R for ease of use despite not knowing when exactly the last mut was #First index individual
    mutation_to = [0]*K
    mutation_wait = [0]*K

    reproductions_from_type = [0]*K
    coalescence_from_type = [0]*K
    last_rep_rs = [R]*L_R
    reproduction_wait = [0]*K

    for event in reversed(eventlist): #List is iterated over in forward time
        for individual in individuals:
            asgen_type_distribution[individual.typ] += last_r - event.r
        last_r = event.r

        if type(event) == BranchingEvent:
            origin_type_before_branching = individuals[event.origin].typ
            target_type_before_branching = individuals[event.target].typ

        event.apply_event_forward(individuals) #Removing silent events from the list is more time intensive than simply checking for silent events later. Maybe if the typed graph was used for more things it'd make sense.

        #The entire following block is used to determine the actual mutation rate (without silent events) in the graph. It has been added late in development and may be a bit of a hack-job.
        match event:
            case MutationEvent():
                if not event.is_silent:
                    mutation_to[event.typ] += 1

                    for i in range(K):
                        if i != event.previousType:
                            mutation_wait[i] += last_mut_rs[event.target]-event.r
                    last_mut_rs[event.target] = event.r

            case BranchingEvent():
                if not event.is_silent:
                    for i in range(K):
                        if i != target_type_before_branching:
                            mutation_wait[i] += last_mut_rs[event.target]-event.r #This is needed so that the wait times of replaced individuals is not lost.
                    last_mut_rs[event.target] = last_mut_rs[event.origin]
                    reproductions_from_type[origin_type_before_branching] += 1

                else:
                    for i in range(len(mutation_wait)):
                        if i != origin_type_before_branching:
                            mutation_wait[i] += last_mut_rs[event.origin]-event.r #This is needed so that the wait times of deleted individuals is not lost.
                del last_mut_rs[event.origin]

            case CoalescingEvent():
                last_mut_rs.insert(event.target, event.r)
                coalescence_from_type[individuals[event.origin].typ] += 1/(len(individuals))

    for i in range(K): #Add wait times since the last event.
        if i != individuals[0].typ:
            mutation_wait[i] += last_mut_rs[0]


    for individual in individuals:#Finally, add the time from the last event
        asgen_type_distribution[individual.typ] += last_r

    ancestor_type, present_type = individuals[0].ancestor_type, individuals[0].typ

    return asgen_type_distribution, ancestor_type, present_type, mutation_to, mutation_wait, reproductions_from_type, coalescence_from_type


def run_simulation(K, s, th, R, v, dirichlet):
    L_R, eventlist = generate_graph(s, th, R, v)
    potential_ancestor_types = [0]*K

    ancestors = generate_potential_ancestors(L_R, th, v, s, dirichlet)
    for ancestor in ancestors:
        potential_ancestor_types[ancestor.typ] += 1



    asgen_type_distribution, ancestor_type, present_type, mutation_to, mutation_wait, reproductions_from_type, coalescene_from_type = type_graph(ancestors, eventlist, K, R)
    current_line = 0 #The current index of the ancestral line



    current_type = present_type
    L_r = 1
    last_r = 0
    #Ancestral Line status is added to the list AFTER the events are applied, meaning the events, if read in forward time, CHANGE the status given in their array.

    time_in_type = [0] * K
    reproduction_events = [0] * K
    L_r_at_AL_reproduction = [0]*2

    reproduction_to_AL = [0]*K
    reproduction_from_AL = 0

    mutation_to_on_AL = [0]*K
    mutation_wait_on_AL = [0]*K
    later_mutation_r_on_AL = 0

    reproductions_from_type_on_AL = [0]*K
    coalescence_from_type_on_AL = [0]*K
    coalescence_to_type_on_AL = [0]*K




    for event in eventlist:
        time_in_type[current_type] += (event.r - last_r) #The ancestral line spends the time since the last event in its current type.
        last_r = event.r

        match event:
            case MutationEvent():
                if event.target == current_line and not event.is_silent:
                    current_type = event.previousType
                    mutation_to_on_AL[event.typ] += 1
                    for i in range(len(mutation_wait_on_AL)):
                        if i != event.typ: #Add the waiting time to the next event to all other types that were not chosen
                            mutation_wait_on_AL[i] += event.r - later_mutation_r_on_AL #How long it is to the next mutation of this type on the AL.
                    later_mutation_r_on_AL = event.r

            case BranchingEvent():
                if not event.is_silent:#Only count loud events
                    reproduction_events[event.typ] += 1/L_r
                L_r +=1
                if event.target == current_line and not event.is_silent: #If the event targets the current_line and isn't silent:
                    add_L_r(L_r_at_AL_reproduction, L_r)
                    current_line = event.origin
                    reproduction_to_AL[event.typ] += 1
                    reproductions_from_type_on_AL[current_type] += 1
                #Branchingevents never occur FROM the ancestral line.


            case CoalescingEvent():
                reproduction_events[K-1] += 1/L_r
                L_r -=1 #This doesn't change the course of the action.

                if event.origin == current_line: #If the event comes from the current line
                    reproduction_from_AL += 1
                    coalescence_from_type_on_AL[current_type] += 1
                elif event.target == current_line : #If the event targets the current_line.
                    current_line = event.origin
                    reproduction_to_AL[K-1] += 1
                    coalescence_to_type_on_AL[current_type] += 1
                if event.target < current_line: #If a coalescing event removes an individual above the current line, the current lines' position needs to be adjusted.
                    current_line -= 1


    time_in_type[current_type] += R-last_r #finally, add the time from the last (in backwards time) event to R

    for i in range(len(mutation_wait_on_AL)): #Also add the waiting time to the first event.
        if i != current_type:
            mutation_wait_on_AL[i] += R-later_mutation_r_on_AL

    return time_in_type, ancestor_type, present_type, asgen_type_distribution, reproduction_events, potential_ancestor_types, L_r_at_AL_reproduction, L_R, mutation_to_on_AL, mutation_wait_on_AL, reproduction_to_AL, reproduction_from_AL, mutation_to, mutation_wait, reproductions_from_type, reproductions_from_type_on_AL, coalescene_from_type, coalescence_from_type_on_AL, coalescence_to_type_on_AL



#Main -----------------------------------------------------------------------

def main(s, th, R, v, repeats, dirichlet):
    K = len(s)
    if K != len(v):
        sys.exit("The list of mutationchances (-v) needs to be exactly one element longer than the input list of reproductive advantages (s).")
    if math.fsum(v) != 1:
        sys.exit("The sum of the list of mutationchances (-v) needs to be exactly 1.")
    for i in range(len(s)-1):
        if s[i]<=s[i+1]:
            sys.exit("Reproductive advantages (-s) values need to be ordered, such that s[i]>s[i+1] for 0 <= i <= K-2")
    if th < 0:
        sys.exit("Mutation rate (-t) has to be non-negative.")
    if any(x < 0 for x in s):
        sys.exit("All reproductive advantages (elements of -s) have to be non-negative.")
    if any(x < 0 for x in v):
        sys.exit("All mutation probabilities (elements of -v) have to be non-negative.")
    if R <= 0:
        sys.exit("The length of the graph (-R) needs to be positive.")
    if repeats < 1:
        sys.exit("The number of repeats (-i) needs to be an integer larger than 0.")

    time_in_type_total  = [0] * K
    ancestral_type_distribution = [0] * K

    mutation_to_on_AL_total = [0]*K
    mutation_wait_on_AL_total = [0]*K
    mutation_to_total = [0]*K
    mutation_wait_total = [0]*K

    stationary_type_distribution = [0] * K
    asgen_type_distribution_total = [0]*K
    potential_ancestor_types_total = [0] * K

    reproduction_total = [0] * K
    reproduction_to_AL_total = [0]*K
    reproduction_from_AL_total = 0

    L_r_at_reproduction_total = [0]*2
    L_R_dist = []

    reproductions_from_type_on_AL_total = [0]*K
    reproductions_from_type_total = [0]*K

    asgen_time_spent_total= [0]*K

    coalescence_from_type_total = [0]*K
    coalescence_from_type_on_AL_total = [0]*K
    coalescence_to_type_on_AL_total = [0]*K

    for i in range(repeats):
        print("Repeat number", i)
        time_in_type, ancestor_type, present_type, asgen_type_distribution, reproduction_events, potential_ancestor_types, L_r_at_AL_reproduction, L_R, mutation_to_on_AL, mutation_wait_on_AL, reproduction_to_AL, reproduction_from_AL, mutation_to, mutation_wait, reproductions_from_type, reproductions_from_type_on_AL, coalescence_from_type, coalescence_from_type_on_AL, coalescence_to_type_on_AL = run_simulation(K, s, th, R, v, dirichlet)

        time_in_type_total = [sum(x) for x in zip(time_in_type_total, time_in_type)]

        asgen_time_spent_total = [sum(x) for x in zip(asgen_time_spent_total, asgen_type_distribution)]
        asgen_type_distribution = [x/sum(asgen_type_distribution) for x in asgen_type_distribution]
        asgen_type_distribution_total = [sum(x) for x in zip(asgen_type_distribution_total, asgen_type_distribution)]
        potential_ancestor_types_total = [sum(x) for x in zip(potential_ancestor_types_total, potential_ancestor_types)]
        reproduction_total = [sum(x) for x in zip(reproduction_total, reproduction_events)]
        L_r_at_reproduction_total = sum_lists_of_different_lengths(L_r_at_reproduction_total, L_r_at_AL_reproduction)
        add_L_r(L_R_dist, L_R)

        mutation_to_on_AL_total = [sum(x) for x in zip(mutation_to_on_AL_total, mutation_to_on_AL)]
        mutation_wait_on_AL_total = [sum(x) for x in zip(mutation_wait_on_AL_total, mutation_wait_on_AL)]

        mutation_to_total = [sum(x) for x in zip(mutation_to_total, mutation_to)]
        mutation_wait_total = [sum(x) for x in zip(mutation_wait_total, mutation_wait)]

        reproduction_to_AL_total = [sum(x) for x in zip(reproduction_to_AL_total, reproduction_to_AL)]
        reproductions_from_type_on_AL_total = [sum(x) for x in zip(reproductions_from_type_on_AL_total, reproductions_from_type_on_AL)]
        reproductions_from_type_total = [sum(x) for x in zip(reproductions_from_type_total, reproductions_from_type)]

        coalescence_from_type_on_AL_total = [sum(x) for x in zip(coalescence_from_type_on_AL_total, coalescence_from_type_on_AL)]
        coalescence_from_type_total = [sum(x) for x in zip(coalescence_from_type_total, coalescence_from_type)]
        coalescence_to_type_on_AL_total = [sum(x) for x in zip(coalescence_to_type_on_AL_total, coalescence_to_type_on_AL)]


        stationary_type_distribution[present_type] += 1
        ancestral_type_distribution[ancestor_type] += 1
        reproduction_from_AL_total += reproduction_from_AL

    print("\nOutput:----------------------------------------------------------------------------------------------\n")
    print("Type distribution on the ancestral line:::::::::::")
    for i in range(K):
        print("Ancestral line spent time in type", i, ":", (time_in_type_total[i]/sum(time_in_type_total))*100,"%.") #There'll be negligible floating point errors with all of hese values.

    print("Ancestral type distribution:::::::::::")
    if K ==2:
        print("Calculated ATD:", altbigp(th, s[0], v))
    print("Simulated ATD:")
    for i in range(K):
        print("Type", i, ":", (ancestral_type_distribution[i]/sum(ancestral_type_distribution))*100,"%.")

    print("\nTime-averaged type distribution of the entire graph:::::::::::")
    for i in range(K):
        print("Graph spent time in type", i, ":", [x/repeats for x in asgen_type_distribution_total][i] * 100,"%.")

    print("Stationary type distribution:::::::::::")
    for i in range(K):
        print("Type", i, ":", (stationary_type_distribution[i]/sum(stationary_type_distribution))*100,"%.")
    if K ==2:
        print("Expected Wright-distribution:", wright_expected_value(s[0], th, v))


    print("Type distribution of potential ancestors:::::::::::::::")
    for i in range(K):
        print("Type", i, ":", (potential_ancestor_types_total[i]/sum(potential_ancestor_types_total))*100,"%.")

    print("\nMutations rate per type on the ancestral line:::::::::::::")
    for i in range(K):
        print("To type", i, ":", (mutation_to_on_AL_total[i]/mutation_wait_on_AL_total[i]))

    print("Mutations rate per type in the entire graph:::::::::::::")
    for i in range(K):
        print("To type", i, ":", (mutation_to_total[i]/mutation_wait_total[i]))

    print("\nReproduction events per time unit on the ancestral line. Doesn't include silent ones. Type K-1 is coalescence. Keep in mind that type K-1 is counted to and from the ancestral line separately. To gain the a statistic comparable to the total graph, divide the sum of incoming and outgoing coalescence events by 2, or double the amount of type K-1 reproduction events in the graph. See \"Simulation of the one-locus K-type ancestral selection graph with parent-independent mutation\", 2025 for more details. ::::::::::::")
    print("Average reproduction events per time unit from the ancestral line (Only counts Type K-1 due to a non-silent selective reproductive events being defined as going to, instead of coming from the ancestral line):", reproduction_from_AL_total/(R*repeats))

    print("Average rate of non-silent reproduction events targeting the ancestral line:", sum(reproduction_to_AL_total)/(R*repeats))
    for i in range(K):
        print("Type", i, "event rate:", (reproduction_to_AL_total[i]/(R*repeats)))

    print("Average non-silent reproduction events per time unit per individual in the graph, doesn't include silent ones:")
    for i in range(K):
        print("Type", i, "event rate:", (reproduction_total[i]/(R*repeats)))

    reproduction_rates_on_AL = [reproductions_from_type_on_AL_total[i]/time_in_type_total[i] for i in range(K)]
    print("\nRate of non-silent selective reproduction events while the graph/the ancestral line are in a specific type:::::::::::::::::::")
    print("Average rate of non-silent selective reproduction events on the ancestral line while it is in a certain type:")
    for i in range(K):
        print("While in type", i, ":", reproduction_rates_on_AL[i])

    reproduction_rates = [reproductions_from_type_total[i]/asgen_time_spent_total[i] for i in range(K)]
    print("Average rate of non-silent selective reproduction events in the graph per individual while that individual is in a certain type:")
    for i in range(K):
        print("While in type", i, ":", reproduction_rates[i])

    coalescence_rates = [coalescence_from_type_total[i]/(R*repeats*[x/repeats for x in asgen_type_distribution_total][i]) for i in range(K)]
    coalescence_rates_from_AL = [coalescence_from_type_on_AL_total[i]/time_in_type_total[i] for i in range(K)]
    coalescence_rates_to_AL = [coalescence_to_type_on_AL_total[i]/time_in_type_total[i] for i in range(K)]
    print("\nRate of coalescing events while the graph/the ancestral line are in a specific type::::::::::")
    print("Average rate of coalescence events originating from the ancestral line while it is in a certain type:")
    for i in range(K):
        print("While in Type", i, ":", coalescence_rates_from_AL[i])

    print("Average rate of coalescence events targeting the ancestral line while it is in a certain type:")
    for i in range(K):
        print("While in type", i, ":", coalescence_rates_to_AL[i])

    print("Average rate of coalescence events in the graph per individual while that individual is in a certain type:")
    for i in range(K):
        print("While in type", i, ":", coalescence_rates[i])


    reproduction_L_r_dist = [round(x/sum(L_r_at_reproduction_total), 6) for x in L_r_at_reproduction_total]
    print("\nLine allocations. All of them start at L_r=1:::::::::::::::::::::::::::::::::::::::::::::::::::")
    print("L_r distribution immediately before selective reproduction events on the ancestral line, starting with L_r=1\t:", reproduction_L_r_dist[1:])
    analyt_L_r = []
    for i in range(len(L_r_at_reproduction_total)):
        analyt_L_r.append(round(L_inf_dist(i, s[0]), 6))
    print("Analytical distribution\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t:", analyt_L_r[1:])
    print("Absolute difference between L_r before selective events and analytical distribution\t\t\t\t\t\t\t\t:", [round(abs(ai-bi),6) for ai ,bi in zip(reproduction_L_r_dist,analyt_L_r)][1:])
    print("L_R distribution\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t:", [round(x/sum(L_R_dist),6) for x in L_R_dist][1:])

    mutation_rates_on_AL = [mutation_to_on_AL_total[i]/mutation_wait_on_AL_total[i] for i in range(len(mutation_to_on_AL_total))]
    mutation_rates = [mutation_to_total[i]/mutation_wait_total[i]for i in range(len(mutation_to_total))]

    return asgen_type_distribution_total, time_in_type_total, ancestral_type_distribution, stationary_type_distribution, mutation_rates_on_AL, mutation_rates, reproduction_to_AL_total, reproduction_from_AL_total, reproduction_total, reproduction_rates_on_AL, reproduction_rates, coalescence_rates_to_AL, coalescence_rates_from_AL, coalescence_rates




def parse_args():#I considered simply sorting the list of σs to make it easier on the user, but then they'd need to give the νs in the same order anyway, at which point they could take the time to sort their inputs.
    parser=ArgumentParser(description="This program simulates ancestral selection graphs with a single mutation locus and any number of types with different reproductive advantages.")
    parser.add_argument("-input", "-i", dest="repeats",  required=True, type=int, help="The number of graphs simulated.")
    parser.add_argument("-reproductive_advantages", "-sigmas", "-s", dest="s", nargs="*", type=float, required=True, help="A list of the reproductive advantages of the different types in the diffusionlimit, starting with type 0 and ending with type K-2. s_i needs to be higher than s_i+1. K-1 will be automatically added as s_(K-1)=0. Simply write the numbers separated by spaces e.g: \"-s 0.3 0.2 0.1\"")
    parser.add_argument("-mutation_rate", "-theta", "-t", dest="t", type=float, required=True, help="The rate with which an individual mutates in the diffusionlimit.")
    parser.add_argument("-time", "-R", dest="R", type=float, required=True, help="The length of time graph will be simulated for.")
    parser.add_argument("-mutation_probabilities", "-nus", "-v", dest="v", nargs="*", type=float, required=True, help="A list of the probabilities with which an individual mutates into a specific type, if a mutation occurs, starting with type 0 and ending with type K-1. The sum of these values needs to be 1. Simply write the numbers separated by spaces e.g: \"-v 0.2 0.3 0.4 0.1\"")
    parser.add_argument("-dirichlet", "-d", dest="dirichlet",action='store_true', required=False, help="If this command is added a Dirichlet-distribution is used to determine potential ancestors instead of a Wright-distribution. This may be useful for certain parameter configurations, as the rejection sampling process used for the Wright-distribution may become extremely unlikely to produce results. The type distribution on the ancestral line should still provide a very close approximation of the ancestral type distribution.")
    args=parser.parse_args()
    return args




if __name__ == "__main__":

    inputs = parse_args()
    sigmas, theta, R, nus, repeats = inputs.s, inputs.t, inputs.R, inputs.v, int(inputs.repeats)
    sigmas.append(0) #Add the reproductive advantage of s_(K-1), which is 0.
    main(s=sigmas, th=theta, R=R, v=nus, repeats=repeats, dirichlet=inputs.dirichlet)




