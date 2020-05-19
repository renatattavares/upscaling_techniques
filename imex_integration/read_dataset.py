"""
Module to interpret IMEX dataset of structured meshes and generate a h5m mesh file, legible to IMPRESS
"""
import numpy as np
from imex_integration.interpreter import Interpreter
from imex_integration.mesh_constructor import MeshConstructor

def read_dataset(dataset_file):

    with open(dataset_file, 'r') as dataset:
        lines = dataset.readlines()

    number_elements, length_elements = Interpreter.read_mesh_data(lines)
    porosity_data = Interpreter.read_porosity(lines)
    permeability_data = Interpreter.read_permeability(lines)
    mesh_file = 'mesh/generated_mesh.h5m'
    generator = MeshConstructor(number_elements, length_elements, mesh_file)

    return porosity_data, permeability_data, number_elements, length_elements
