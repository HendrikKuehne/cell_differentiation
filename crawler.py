import numpy as np
import os

outfile = open("simulation_output/metadata.txt",mode="w")

# going into every folder in "simulation_results" and getting the parameters of the simulation
for path in os.listdir("simulation_output"):
    if not os.path.isfile("simulation_output/" + path):
        # fetching the simulation parameters
        file = open("simulation_output/" + path + "/abstract.txt","r")
        for line in file:
            line = line.replace(" ","")
            line = line.replace("\n","")
            line = line.split("=")

            if(line[0] == "Nmax"):
                Nmax = int(line[1])
            elif(line[0] == "theta"):
                n_genes = len(line[1].split(","))
            elif(line[0] == "tdiv"):
                tdiv = int(line[1])                 # ill-named: This is the number of updates to the expresion levels that are carried out, not
                                                    # the time in seconds between cell divisions
            elif(line[0] == "dt"):
                dt = float(line[1])
        file.close()

        # writing to outfile
        outfile.write(path + " :\n    Nmax = {}\n    tdiv = {}\n    dt = {}\n".format(Nmax,tdiv,dt))

outfile.close()