"""
Module designed to read specific keywords from IMEX data set to perform an upscaling

NECESSARY DATA:

-> Mesh
-> Porosity
-> Permeability
"""

grid_token = '*GRID'
di_token = '*DI'
dj_token = '*DJ'
dk_token = '*DK'

tokens = [grid_token, di_token, dj_token, dk_token]
words = []

with open('../imex_datasets/teste.dat', 'r') as dataset:

    i = 0
    lines = dataset.readlines()

    for token in tokens:
        for line in lines:
             if line.find(token) is 0:
                 correct_line = line
                 break
        words.append(correct_line.split(" "))
        i =+ i

    nx, ny, nz = words[0][2], words[0][3], words[0][4]
    lx, ly, lz = words[1][2], words[2][2], words[3][2]
