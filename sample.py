import numpy as np
from scipy.spatial.transform import Rotation

class Sample:
    def __init__(self, position, orientation, cell_parameters, q_range,
                 cell_orientations, cell_colors=None, q_probe = 1,
                 sample_scale = 1):
        '''
        position: lab frame x,y,z coordiantes
        orientation: sample
        '''
        self.name = 'sample'
        self.position = np.asarray(position)
        self.orient_array = orientation
        self.cell_parameters = cell_parameters
        self.sample_scale = sample_scale
        self.a, self.b, self.c = cell_parameters
        self.v = self.a*self.b*self.c
        self.q_range = q_range
        self.q_probe = q_probe
        self.cell_orientations = [Rotation.from_euler('ZXZ', co, degrees = True)
                                  for co in cell_orientations]
        self.get_reciprocal_lattice_vectors()
        self.generate_reciprocal_lattice_to_q()
        self.orient_reciprocal_lattices_to_sample()
        if cell_colors == None:
            self.cell_colors = ['black']*len(cell_orientations)
        else:
            self.cell_colors = cell_colors
        self.update()

    def __repr__(self):
        return self.name
    def get_view_in_lab_frame(self, view_axis):
        return self.rot.apply(view_axis)
    def update(self):
        self.rot = Rotation.from_euler('ZXZ', self.orient_array, degrees = True)
        self.lab_frame_axes = self.rot.as_matrix()
        self.orient_reciprocal_lattices_to_lab()
    def get_reciprocal_lattice_vectors(self):
        a, b, c = np.array((self.a,0,0)), np.array((0,self.b,0)), np.array((0,0,self.c))
        v = self.v
        astar = np.cross(b,c)/v
        bstar = np.cross(c,a)/v
        cstar = np.cross(a,b)/v
        self.astar = astar[0]
        self.bstar = bstar[1]
        self.cstar = cstar[2]
    def generate_reciprocal_lattice_to_q(self):
        arange = np.ceil(self.q_range / self.astar)
        brange = np.ceil(self.q_range / self.bstar)
        crange = np.ceil(self.q_range / self.cstar)
        astar_inds = np.arange(-arange, arange+1)
        bstar_inds = np.arange(-brange, brange+1)
        cstar_inds = np.arange(-crange, crange+1)
        span = int((2*arange+1)*(2*brange+1)*(2*crange+1))
        astars = astar_inds*self.astar
        bstars = bstar_inds*self.bstar
        cstars = cstar_inds*self.cstar
        self.reciprocal_lattice = np.asarray(np.meshgrid(astars,bstars,cstars)).reshape((3, span))
        self.rl_qs = np.linalg.norm(self.reciprocal_lattice, axis = 0)
        q_mask = self.rl_qs < self.q_range
        self.rl_qs = self.rl_qs[q_mask]
        self.reciprocal_lattice = self.reciprocal_lattice[:,q_mask]
        self.rl_alphas = np.where(np.abs(self.rl_qs - self.q_probe) < 0.1, 5, 0.25)
        self.q_of_interest =  np.where(np.abs(self.rl_qs - self.q_probe)<0.1, True, False)
    def orient_reciprocal_lattices_to_sample(self):
        self.sample_space_rlatts = [co.apply(self.reciprocal_lattice.T).T for co in self.cell_orientations]
    def orient_reciprocal_lattices_to_lab(self):
        self.lab_space_rlatts = np.asarray([self.rot.apply(rl.T).T for rl in self.sample_space_rlatts])