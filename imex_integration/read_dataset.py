"""
Module to interpret IMEX dataset of structured meshes and generate a h5m mesh file, legible to IMPRESS
"""
import numpy as np
from pymoab import core, types, rng
from imex_integration.interpreter import Interpreter
from imex_integration.mesh_constructor import MeshConstructor

def read_dataset(dataset_file):

    number_elements, length_elements = Interpreter.read_mesh_data(dataset_file)
    generator = MeshConstructor(number_elements, length_elements)
    print('\n##### Mesh file created #####')
