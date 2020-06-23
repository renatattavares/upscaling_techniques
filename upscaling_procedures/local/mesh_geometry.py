import numpy as np

class MeshGeometry:

    def check_directions_and_coarse_face(self):

        direction = self.directions_dictionary.get(self.direction)

        if np.array_equal(direction, self.x) is True:
            self.perpendicular_direction_1 = self.y
            self.perpendicular_direction_2 = self.z
            self.number_faces_coarse_face = int((self.number_elements_y_direction/self.ny)*(self.number_elements_z_direction/self.nz))

        elif np.array_equal(direction, self.y) is True:
            self.perpendicular_direction_1 = self.x
            self.perpendicular_direction_2 = self.z
            self.number_faces_coarse_face = int((self.number_elements_x_direction/self.nx)*(self.number_elements_z_direction/self.nz))

        elif np.array_equal(direction, self.z) is True:
            self.perpendicular_direction_1 = self.x
            self.perpendicular_direction_2 = self.y
            self.number_faces_coarse_face = int((self.number_elements_x_direction/self.ny)*(self.number_elements_y_direction/self.ny))

    def areas(self):
        pass
