#THis script finds all dark particles and their daughter particles and puts them in their
#respective lists, what we need to do now is use the particles that are stored as lists of strings
#which is basically their line in the hepmc file, the 7th element in each particle list should
#be the energy of that particle, we need to create histograms of particles vs energy

import sys
import matplotlib.pyplot as plt
import numpy as np
from Particle import Particle

dark_particles = []
daughter_particles = []
vertex_number = None
event_num = None

#this functions finds potential daughter particles for dark particles
#if found add it to daughter particles list
def check_daughters(file_iter, num):
    global event_num, vertex_number
    for _ in range(num):
        try:
            line = next(file_iter)
        except StopIteration:
            break

        data = line.split()
        if data[2] == '21':
            daughter_particles.append(Particle(data,vertex_number,event_num))
        elif data[2] == '1023':
            dark_particles.append(Particle(data,vertex_number,event_num))
            vertex_number = data[11]

#This funcion finds every 1023 electron which is what we are lloking for
#and adds it the dark particles list, when it finds a 1023 particle
#it assign vertex to looks for "daughter" particles
def find_electron(file_iter, data):
    global vertex_number, event_num

    if data[0] == 'E':
        event_num = data[1]

    elif data[0] == 'P' and data[2] == '1023':
        dark_particles.append(Particle(data,vertex_number,event_num))
        vertex_number = data[11]


    elif data[0] == 'V' and data[1] == vertex_number:
        num_daughters = int(data[8])
        check_daughters(file_iter, num_daughters)

def plot_histograms():
    #plots the 1023 particles
    data = [np.float64(x.data[7]) for x in dark_particles]

    plt.hist(data, bins=10)  # bins is optional
    plt.xlabel("Energy of Particles")
    plt.ylabel("Frequency")
    plt.title("Histogram of 1023 Particles")
    plt.show()

    #plots the daughter particles
    data = [np.float64(x.data[7]) for x in daughter_particles]

    plt.hist(data, bins=10)  # bins is optional
    plt.xlabel("Energy of Particles")
    plt.ylabel("Frequency")
    plt.title("Histogram of daughter Particles")
    plt.show()

#opens the hep mc file and reads line by line
with open(sys.argv[1]) as f:
    file_iter = iter(f)

    for line in file_iter:
        data = line.split()
        find_electron(file_iter, data)

plot_histograms()

#for particle in dark_particles:
#    if float(particle[7]) < 0:
#        print(particle)
#for particle in daughter_particles:
    #print(f"E: {particle.event_num} V: {particle.vertex_num} {particle.data[7]}")
