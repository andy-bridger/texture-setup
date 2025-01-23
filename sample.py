import numpy as np
from scipy.spatial.transform import Rotation
import alphashape

class SampleObject:
    def __init__(self, position, smp, mesh = None, mesh_scale = 0.1):
        self.position = np.asarray(position)
        self.smp = smp
        self.mesh = mesh
        self.original_mesh_vectors = mesh.vectors
        self.scale = mesh_scale
        self.primitive = False

    def get_mesh_vectors(self, scale = 1):
        return np.array(list(
            map(Rotation.from_matrix(self.smp.get_rot()).apply, (self.original_mesh_vectors.copy()+ self.smp.get_translation()) * self.scale * scale)))

    def get_internal_path_length(self, ki, kd):
        # assume the beam is incident on 0,0,0
        mv = self.get_mesh_vectors()[:,0,:]
        mv_norm = mv/np.linalg.norm(mv, axis = 1)[:,None]
        vert_ki = np.argmax(np.dot(-ki, mv_norm.T)) #closest point to perfect alignment with ki
        vert_kd = np.argmax(np.dot(kd, mv_norm.T))  # closest point to perfect alignment with kd
        p_ki, p_kd = mv[vert_ki], mv[vert_kd]
        return np.linalg.norm(p_ki) + np.linalg.norm(p_kd), p_ki, p_kd

    def set_smp(self,smp):
        self.smp = smp

class Sample(SampleObject):
    def __init__(self, position, goniometer, cell_parameters, q_range,
                 cell_orientations, cell_colors=None, q_probe = 1,
                 sample_scale = 1, mesh = None):
        '''
        position: lab frame x,y,z coordiantes
        orientation: sample
        '''
        super().__init__(position, goniometer, mesh)
        self.smp = goniometer
        if type(self.mesh) == type(None):
            self.primitive = True
        self.name = 'sample'
        self.position = np.asarray(position)
        self.goniometer = goniometer
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
        self.rot = self.goniometer.rot
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
    def adjust_alphas(self, mask):
        self.rl_alphas = np.where(mask, 5, 0.25)
    def get_q_shell(self, q_probe, thresh = 0.1):
        self.q_of_interest =  np.where(np.abs(self.rl_qs - q_probe)<thresh, True, False)
        return self.q_of_interest
    def orient_reciprocal_lattices_to_sample(self):
        self.sample_space_rlatts = [co.apply(self.reciprocal_lattice.T).T for co in self.cell_orientations]
    def orient_reciprocal_lattices_to_lab(self):
        self.lab_space_rlatts = np.asarray([self.rot.apply(rl.T).T for rl in self.sample_space_rlatts])
