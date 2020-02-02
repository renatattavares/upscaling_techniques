import numpy as np

class Interpreter:

  def read_mesh_data(dataset_file):

      print('\n##### Read mesh data from dataset #####')
      grid_token = '*GRID'
      di_token = '*DI'
      dj_token = '*DJ'
      dk_token = '*DK'

      tokens = [grid_token, di_token, dj_token, dk_token]
      words = []

      with open(dataset_file) as dataset:

          i = 0
          lines = dataset.readlines()

          for token in tokens:
              for line in lines:
                   if line.find(token) is 0:
                       correct_line = line
                       break
              words.append(correct_line.split(" "))
              i =+ i

          number_elements_x, number_elements_y, number_elements_z = int(words[0][2]), int(words[0][3]), int(words[0][4])
          lenght_elements_x, lenght_elements_y, lenght_elements_z = int(words[1][2]), int(words[2][2]), int(words[3][2])

          number_elements = np.array([number_elements_x, number_elements_y, number_elements_z])
          length_elements = np.array([lenght_elements_x, lenght_elements_y, lenght_elements_z])

          print('Number of elements in x direction: {}'.format(number_elements_x))
          print('Number of elements in y direction: {}'.format(number_elements_y))
          print('Number of elements in z direction: {}'.format(number_elements_z))

          print('\nLength of blocks in x direction: {}'.format(lenght_elements_x))
          print('Length of blocks in y direction: {}'.format(lenght_elements_y))
          print('Length of blocks in z direction: {}'.format(lenght_elements_z))

          return number_elements, length_elements

    def read_porosity():
        pass

    def read_permeability():
        pass
