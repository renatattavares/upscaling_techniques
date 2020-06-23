import numpy as np

class MeshGeometry:

    def check_coarse_face(self):

        self.number_faces_coarse_face = np.array([], dtype = int)

        self.number_faces_coarse_face = np.append(self.number_faces_coarse_face,  int((self.number_elements[1]/self.ny)*(self.number_elements[2]/self.nz)))

        self.number_faces_coarse_face = np.append(self.number_faces_coarse_face, int((self.number_elements[0]/self.nx)*(self.number_elements[2]/self.nz)))

        self.number_faces_coarse_face = np.append(self.number_faces_coarse_face, int((self.number_elements[0]/self.ny)*(self.number_elements[1]/self.ny)))


    def center_distance_walls(self):

        self.center_distance_walls_x = np.array([])
        self.center_distance_walls_y = np.array([])
        self.center_distance_walls_z = np.array([])

        for i in range(self.number_coarse_volumes):
            min_coord = self.coarse.elements[i].volumes.center[:].min(axis = 0)
            max_coord = self.coarse.elements[i].volumes.center[:].max(axis = 0)
            self.center_distance_walls_x = np.append(self.center_distance_walls_x, (max_coord[0] - min_coord[0]))
            self.center_distance_walls_y = np.append(self.center_distance_walls_x, (max_coord[1] - min_coord[1]))
            self.center_distance_walls_z = np.append(self.center_distance_walls_x, (max_coord[2] - min_coord[2]))

    def areas(self):

        self.areas = np.array([])

        self.areas = np.append(self.areas, self.length_elements[1]*self.length_elements[2])
        self.areas = np.append(self.areas, self.length_elements[0]*self.length_elements[2])
        self.areas = np.append(self.areas, self.length_elements[0]*self.length_elements[1])
             
