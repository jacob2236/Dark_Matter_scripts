#THis script finds all dark particles and their daughter particles and puts them in their
#respective lists, what we need to do now is use the particles that are stored as lists of strings
#which is basically their line in the hepmc file, the 7th element in each particle list should
#be the energy of that particle, we need to create histograms of particles vs energy

import sys

dark_particles = []
daughter_particles = []
vertex_number = None

#this functions finds potential daughter particles for dark particles
#if found add it to daughter particles list
def check_daughters(file_iter, num):
    for _ in range(num):
        try:
            line = next(file_iter)
        except StopIteration:
            break

        data = line.split()
        if data[2] == '11' or data[2] == '-11':
            daughter_particles.append(data)
        elif data[2] == '1023':
            dark_particles.append(data)
            vertex_number = data[11]
            break

#This funcion finds every 1023 electron which is what we are lloking for
#and adds it the dark particles list, when it finds a 1023 particle
#it assign vertex to looks for "daughter" particles
def find_electron(file_iter, data):
    global vertex_number

    if data[0] == 'P' and data[2] == '1023':
        dark_particles.append(data)
        vertex_number = data[11]


    elif data[0] == 'V' and data[1] == vertex_number:
        num_daughters = int(data[8])
        check_daughters(file_iter, num_daughters)

#opens the hep mc file and reads line by line
with open(sys.argv[1]) as f:
    file_iter = iter(f)

    for line in file_iter:
        data = line.split()
        find_electron(file_iter, data)

for particle in dark_particles:
    print(particle)
for particle in daughter_particles:
    print(particle)
