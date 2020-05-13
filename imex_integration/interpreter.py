import numpy as np

class Interpreter:

    def read_mesh_data(lines):

        print('\n##### Read mesh data from dataset #####\n')
        grid_token = '*GRID *CART'
        di_token = '*DI'
        dj_token = '*DJ'
        dk_token = '*DK'
        directions = np.array(['x','y','z'])
        lines = lines

        tokens = [grid_token, di_token, dj_token, dk_token]
        grid_info = []
        info = []
        start = False

        for line in lines:
            if grid_token in line:
                start = True
            if start is True:
                for token in tokens:
                    if token in line:
                        grid_info.append(line.split('**')[0])
            if dk_token in line:
                break

        for line in grid_info:
            for letter in line:
                if letter.isdigit() is True:
                    pos = line.index(letter)
                    info.append(line[pos:].strip())
                    break

        #number_elements_x, number_elements_y, number_elements_z = info[0].split(' ')
        number_elements = [int (i) for i in info[0].split(' ')]
        lenght_elements_x = int(info[1])
        lenght_elements_y = int(info[2])
        lenght_elements_z = int(info[3])

        #number_elements = np.array([number_elements_x, number_elements_y, number_elements_z])
        number_elements = np.array(number_elements)
        length_elements = np.array([lenght_elements_x, lenght_elements_y, lenght_elements_z])

        for i,j in zip(directions, number_elements):
            print('Number of elements in {} direction: {}'.format(i,j))

        for i,j in zip(directions, length_elements):
            print('Length of elements in {} direction: {}'.format(i,j))

        return number_elements, length_elements

    def read_porosity(lines):

        print('\n##### Read porosity data from dataset #####')
        porosity_token = '*POR'
        porosity_lines = int(1122000/5)
        lines = lines
        i = 0

        for line in lines:
            if line.find(porosity_token) == 0:
                porosity_data = lines[(i+1):int(i+1+porosity_lines)]
                break
            i += 1

        porosity_array = np.zeros((porosity_lines,5))

        for i in range(porosity_lines):
            porosity_array[i] = np.fromstring(porosity_data[i], sep = " ")

        porosity = porosity_array.flatten()

        return porosity

    def read_permeability(lines):
        print('\n##### Read permeability data from dataset #####')
        permeability_x_token = '*PERMI'
        permeability_y_token = '*PERMJ'
        permeability_z_token = '*PERMK'
        permeability_data = []
        lines = lines

        found = False

        for line in lines:
            if permeability_x_token in line:
                indexx = lines.index(line)

            if permeability_y_token in line:
                indexy = lines.index(line)

            if permeability_z_token in line:
                indexz = lines.index(line)

            if '*MODEL *OILWATER' in line:
                end = lines.index(line)

        permi = lines[indexx+1:indexy-1]
        permj = lines[indexy+1:indexz-1]
        permk = lines[indexz+1:end-1]

        flat_perm = []

        for line in permi:
            l = line.split()
            if l != '\n':
                for value in l:
                    flat_perm.append(float(value))

        permi_array = np.array(flat_perm, dtype = 'float')
        flat_perm = []

        for line in permj:
            l = line.split()
            if l != '\n':
                for value in l:
                    flat_perm.append(float(value))

        permj_array = np.array(flat_perm, dtype = 'float')
        flat_perm = []

        for line in permk:
            l = line.split()
            if l != '\n':
                for value in l:
                    flat_perm.append(float(value))

        permk_array = np.array(flat_perm, dtype = 'float')

        permeability = np.zeros((1122000,3))

        permeability[:,0] = permi_array
        permeability[:,1] = permj_array
        permeability[:,2] = permk_array

        return permeability
