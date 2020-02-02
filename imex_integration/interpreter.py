import numpy as np

class Interpreter:

    def read_mesh_data(lines):

        print('\n##### Read mesh data from dataset #####\n')
        grid_token = '*GRID'
        di_token = '*DI'
        dj_token = '*DJ'
        dk_token = '*DK'
        directions = np.array(['x','y','z'])
        lines = lines

        tokens = [grid_token, di_token, dj_token, dk_token]
        words = []

        i = 0

        for token in tokens:
            for line in lines:
                if line.find(token) is 0:
                    correct_line = line
                    break
            words.append(correct_line.split(" "))
            i =+ 1

        number_elements_x, number_elements_y, number_elements_z = int(words[0][2]), int(words[0][3]), int(words[0][4])
        lenght_elements_x, lenght_elements_y, lenght_elements_z = int(words[1][2]), int(words[2][2]), int(words[3][2])

        number_elements = np.array([number_elements_x, number_elements_y, number_elements_z])
        length_elements = np.array([lenght_elements_x, lenght_elements_y, lenght_elements_z])

        for i,j in zip(directions, number_elements):
            print('Number of elements in {} direction: {}'.format(i,j))

        for i,j in zip(directions, length_elements):
            print('Length of elements in {} direction: {}'.format(i,j))

        return number_elements, length_elements

    def read_porosity(dataset):

        print('\n##### Read porosity data from dataset #####')
        porosity_token = '*POR'
        lines = lines

        for line in lines:
            if line.find(porosity_token) is 0:
                correct_line = line
                break

    def read_permeability():
        # print('\n##### Read porosity data from dataset #####')
        # permeability_token = '*PERM'
        #
        # with open(dataset_file) as dataset:
        #     i = 0
        #     lines = dataset.readlines()
        #     for line in lines:
        #         if line.find(permeability_token) is 0:
        #             correct_line = line
        #             break
