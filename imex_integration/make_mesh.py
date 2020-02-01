from mesh_constructor import MeshConstructor

####################### Informações de entrada da malha ########################

nx = 60 # Número de elementos na direção x
ny = 220 # Número de elementos na direção y
nz = 85 # Número de elementos na direção Z
dx, dy, dz= 20.0, 10.0, 2.0 # Tamanho dos elementos nas direções x e y
dim = 2
num_elements = nx*ny*nz
################################################################################

# Inicializando a malha
MeshConstructor(nx,ny,nz,dx,dy,dz,num_elements)
print("Mesh created")
