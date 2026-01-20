import sys

def find_electron(data):
    if data[0] == 'P' and (data[2] == '11' or data[2] == '-11'):
        print(data)
    #Right now this is just finding electrons going line by line in the hepmc file and printing
    # those lines
    #But we need to look into printing the energy pietro wants.
    #

with open(sys.argv[1]) as f:
    for line in f: find_electron(line.split())

