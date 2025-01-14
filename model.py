import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation
from goniometer import Goniometer
from detector import Detector
from source import Source
from sample import Sample
from helper_funcs import *



class Model():
    def __init__(self, source, detectors, sample, goniometer, ewald_radii, pole_radius,
                 sample_view_axis, detector_colors = None, ewald_steps = 2, q_probe = 1):
        self.source = source
        self.detectors = detectors
        self.sample = sample
        self.goniometer = goniometer
        self.r_dict = {'recip_sample': 1,
                                   'lab_sample': 1,
                                   'k_vecs': 1,
                                   'gonio_r': 1,
                                   'gonio_v': 1,
                                    'q_probe': q_probe}
        self.ewald_radii = ewald_radii
        self.ewald_radius = ewald_radii[0]
        self.q_range = np.linspace(ewald_radii[0], ewald_radii[1], ewald_steps)
        self.pole_radius = pole_radius
        self.sample_view_axis = np.asarray(sample_view_axis)
        self.no_pfi = True
        if detector_colors == None:
            self.detector_colors = ['black']*len(detectors)
        else:
            self.detector_colors = detector_colors
        self.ki_raw = self.sample.position - self.source.position
        self.ki_raw_scale = np.linalg.norm(self.ki_raw)
        self.ki = self.ki_raw/self.ki_raw_scale
        for det in self.detectors:
            det.ki = self.ki
            det.radii = ewald_radii
            det.r_steps = ewald_steps
        self.get_detector_Ks()
        self.get_probe_ind()
        self.setup_qrange_rl()
        self.setup_pole()

        self.show_lab_K_vecs = True
        self.lab_K_vecs = []
        self.setup_inplane_pole_figure()

        self.update()

    def get_probe_ind(self):
        self.probe_ind = np.searchsorted(self.q_range, self.r_dict['q_probe'], 'left')

    def update(self):
        self.pole_view_axis = self.sample.get_view_in_lab_frame(self.sample_view_axis)
        self.north_pole = self.pole_view_axis/np.linalg.norm(self.pole_view_axis)
        self.south_pole = -self.north_pole

        self.calc_pole()
        self.get_detector_signal()
        self.get_probe_ind()
        self.get_inplane_pole_figure()

    def get_sphere_cart_array(self, r = 1, res = 100, offset=(0,0,0)):
        u = np.linspace(0, 2 * np.pi, res)
        v = np.linspace(0, np.pi, res)
        x = r*np.outer(np.cos(u), np.sin(v)) + offset[0]
        y = r*np.outer(np.sin(u), np.sin(v))+ offset[1]
        z = r*np.outer(np.ones(np.size(u)), np.cos(v))+ offset[2]
        return np.concatenate((x[None,:],y[None,:],z[None,:]), axis = 0)


    def cuboid_data(self,center, size, rotation_matrix):
        """
        Generate the X, Y, Z coordinates for a cuboid.

        Parameters:
            center: tuple of 3 floats
                Coordinates of the cuboid's center (Px, Py, Pz).
            size: tuple of 3 floats
                Dimensions of the cuboid (a, b, c).
            rotation_matrix: 3x3 numpy array
                Rotation matrix to orient the cuboid.

        Returns:
            X, Y, Z: 2D numpy arrays
                Coordinates for the cuboid's surface.
        """
        # Extract cuboid size
        a, b, c = size
        # Generate corner points of the cuboid in local coordinates
        x = np.array([-0.5, 0.5]) * a
        y = np.array([-0.5, 0.5]) * b
        z = np.array([-0.5, 0.5]) * c
        # Create a meshgrid for the cuboid
        x, y, z = np.meshgrid(x, y, z, indexing="ij")

        # Flatten and combine into an array of points
        points = np.array([x.flatten(), y.flatten(), z.flatten()]).T

        # Apply the rotation
        rotated_points = points @ rotation_matrix.T

        # Translate to the center
        translated_points = rotated_points + np.array(center)

        # Reshape into 2D arrays for plotting
        X = translated_points[:, 0].reshape(2, 2, 2)
        Y = translated_points[:, 1].reshape(2, 2, 2)
        Z = translated_points[:, 2].reshape(2, 2, 2)

        return X, Y, Z

    def get_Ks(self):
        norm = lambda x : x/np.linalg.norm(x)
        ki = norm(self.sample.position - self.source.position)
        Ks = []
        for det in self.detectors:
            K = norm(norm(det.position - self.sample.position) - ki)
            Ks.append(K)
            det.K = K
        return Ks

    def get_detector_Ks(self):
        self.Ks = self.get_Ks()
        pKs = []
        for det in self.detectors:
            det_pKs = det.get_polychromatic_ks(self.q_range, self.ki)

            pKs.append(det_pKs)
        self.pKs = np.asarray(pKs)

    def setup_qrange_rl(self):
        self.qs_of_interest = []
        for iq, q_probe in enumerate(self.q_range):
            qs_of_interest = []
            for det in self.detectors:
                qs_of_interest.append(self.sample.get_q_shell(np.linalg.norm(det.pKs[iq]),0.2))
            self.qs_of_interest.append(np.asarray(qs_of_interest).sum(axis = 0)>0)


    def get_detector_signal(self, thresh = 0.1):
        q_space_readouts = []
        for iq, q_probe in enumerate(self.q_range):

            q_of_interest = self.qs_of_interest[iq]
            # detector, grain, coord, reflection  det qrange coord
            pKs = self.pKs[:,iq, :]
            disp_vecs = self.sample.lab_space_rlatts[None, :,:,q_of_interest] - pKs[:,None,:,None]
            dists = np.linalg.norm(disp_vecs, axis = 2)
            relevant_dists = np.where(dists < thresh, np.exp(-500*dists**2), 0)
            q_space_readouts.append(relevant_dists.sum(axis = (1,2)))
        q_space_readouts = np.asarray(q_space_readouts)

        # placeholder to keep other functions working
        self.detector_readout = q_space_readouts.T
        for idet, det in enumerate(self.detectors):
            det.readout = self.detector_readout[idet]


    def setup_pole(self):
        self.pole_cart_array = self.get_sphere_cart_array(r=self.pole_radius, offset=(0, 0, 0))

    def calc_pole(self):
        self.eq_cart = self.orient_to_pole(self.equator(), self.pole_view_axis)
        self.pole_Ks = self.get_pole_K_intercepts()
        self.pole_K_projs = np.asarray([self.get_orientated_stereographic_projection(pK) for pK in self.pole_Ks])
        self.pole_K_projs_vec = np.asarray([self.pole_K_projs[i] - K for i,K in enumerate(self.pole_Ks)])

    def calc_ewald(self, er):
        return self.get_sphere_cart_array(r = er, offset = -self.ki*er)

    def get_orientated_stereographic_projection(self, point):
        spoint = self.orient_to_pole(point, self.north_pole)
        oppole = np.array((0,0,-np.sign(spoint[2])*self.pole_radius))
        svec = spoint - oppole
        u = -oppole[2]/svec[2]
        proj_point = (u*svec) + oppole
        return self.orient_to_pole(proj_point, self.north_pole)

    def equator(self, r = None, res = 100, offset=(0,0,0)):
        if r == None:
            r = self.pole_radius
        u = np.linspace(0, 2 * np.pi, res)
        x = r*np.cos(u) + offset[0]
        y = r*np.sin(u) + offset[1]
        z = np.zeros_like(x) + offset[2]
        return np.concatenate((x[None,:],y[None,:],z[None,:]), axis = 0)

    def orient_to_pole(self, eq, view_axis):
        bisect = view_axis + ((np.array((0,0,1)) - view_axis)/2)
        bisect = (bisect/np.linalg.norm(bisect))*-np.pi
        rot = Rotation.from_rotvec(bisect)
        return rot.apply(eq.T).T

    def orient_to_pole_new(self, eq, view_axis):
        rot = self.sample.rot
        return rot.apply(eq.T).T


    def get_pole_K_intercepts(self):
        Gs = []
        for K in self.Ks:
            roots = sphere_line_intercept(np.zeros(3), K, self.sample.position, self.pole_radius)
            Gs.append(roots[roots!=0]*K)
        return np.asarray(Gs)

    def setup_inplane_pole_figure(self):
        self.pole_view_axis = self.sample.get_view_in_lab_frame(self.sample_view_axis)
        self.north_pole = self.pole_view_axis / np.linalg.norm(self.pole_view_axis)
        self.south_pole = -self.north_pole

        self.calc_pole()
        self.get_detector_signal()
        self.pole_figure_points = self.sample.rot.apply(self.pole_K_projs, inverse=True)
        pfi = []
        for di, det in enumerate(self.detectors):#just read the probe channel for now
            dr = det.readout[self.probe_ind]
            npfi = np.concatenate((self.pole_figure_points[di], np.array((dr,))))
            pfi.append(npfi)
        pfi = np.asarray(pfi)
        self.pole_figure_intensities = pfi[np.argsort(pfi[:,3])]


    def get_inplane_pole_figure(self):
        self.pole_figure_points = self.sample.rot.apply(self.pole_K_projs, inverse=True)
        current_is = self.pole_figure_intensities[:,3]
        for di, det in enumerate(self.detectors):#just read the probe channel for now
            dr = det.readout[self.probe_ind]
            pos = current_is.searchsorted(dr)
            npfi = np.concatenate((self.pole_figure_points[di], np.array((dr,))))
            self.pole_figure_intensities = np.concatenate((self.pole_figure_intensities[:pos,:], npfi[None,:], self.pole_figure_intensities[pos:,:]))



    def get_cube(self, size = (0.1, 0.1, 0.1),
                  offset = (0,0,0), use_gonio = True):
        if use_gonio == True:
            cube_data = self.cuboid_data(offset, size, self.sample.lab_frame_axes)
        else:
            cube_data = self.cuboid_data(offset, size, np.eye(3,3))
        return cube_data

    def get_detector_vectors(self):
        return [det.position - self.sample.position for det in self.detectors]

    def get_lab_frame(self):
        source_repr = self.get_cube(offset = self.source.position, size=[2,1,1], use_gonio= False)
        source_col = ['grey']*3
        det_reprs = []
        det_cols = []
        for idet, det in enumerate(self.detectors):
            det_reprs.append(self.get_cube(offset=det.position, size=[1, 1, 2], use_gonio=False))
            det_cols.append([self.detector_colors[idet]] * 3)
        return source_repr, source_col, det_reprs, det_cols


    def get_sample(self, ratio = 0.1):
        cube = self.get_cube(offset = (0,0, 0), size = self.ewald_radius*np.ones(3)*ratio)
        return cube

    ## under here needs to go at some point:

    def show_ewald(self):
        fig = self.fig
        ax = self.ax
        self.calc_ewald()
        self.plot_sphere(fig,ax, self.es_cart_array)



    def show_pole(self):
        fig = self.fig
        ax = self.ax
        self.calc_pole()
        self.plot_sphere(fig,ax, self.pole_cart_array)
        self.plot_line(fig,ax, self.eq_cart)
        for i, K in enumerate(self.pole_Ks):
            ax.quiver(0,0,0, K[0], K[1], K[2], arrow_length_ratio = 0.05, color = self.detector_colors[i])
            proj_vec = self.pole_K_projs[i]-K
            ax.quiver(K[0], K[1], K[2],
                      proj_vec[0], proj_vec[1], proj_vec[2],
                      arrow_length_ratio = 0.05, color = self.detector_colors[i])



    def plot_pole_figure(self, chosen_ax = None, use_intensities = False):
        if chosen_ax == None:
            chosen_ax = self.ax_new

        if use_intensities == False:
            if np.where(self.sample_view_axis ==1)[0] == 0:
                chosen_ax.scatter(self.pole_figure_points[:,1], self.pole_figure_points[:,2], c = self.detector_colors)
            if np.where(self.sample_view_axis ==1)[0] == 1:
                chosen_ax.scatter(self.pole_figure_points[:,0], self.pole_figure_points[:,2], c = self.detector_colors)
            if np.where(self.sample_view_axis ==1)[0] == 2:
                chosen_ax.scatter(self.pole_figure_points[:,0], self.pole_figure_points[:,1], c = self.detector_colors)
        else:
            if np.where(self.sample_view_axis ==1)[0] == 0:
                chosen_ax.scatter(self.pole_figure_intensities[:,1], self.pole_figure_intensities[:,2],
                                    c = self.pole_figure_intensities[:,3])
            if np.where(self.sample_view_axis ==1)[0] == 1:
                chosen_ax.scatter(self.pole_figure_intensities[:,0], self.pole_figure_intensities[:,2],
                                    c = self.pole_figure_intensities[:,3])
            if np.where(self.sample_view_axis ==1)[0] == 2:
                chosen_ax.scatter(self.pole_figure_intensities[:,0], self.pole_figure_intensities[:,1],
                                    c = self.pole_figure_intensities[:,3])

        eq = self.equator()


        chosen_ax.plot(eq[0], eq[1], c = 'grey')

    def plot_lab_K_vecs(self):
        lab_K_vecs = []
        for idet, det in enumerate(self.detectors):
            k_det = (self.detectors[idet].position - self.sample.position)
            k_det_quiver_kwargs = {'color': 'grey', 'arrow_length_ratio': 0.1, 'alpha': 0.25}
            for i, arg in enumerate(('X', 'Y', 'Z')):
                k_det_quiver_kwargs[arg] = self.sample.position[i]
            for i, arg in enumerate(('U', 'V', 'W')):
                k_det_quiver_kwargs[arg] = k_det[i]
            self.lab_ax.quiver(**k_det_quiver_kwargs)

            k_quiver_kwargs = {'color': self.detector_colors[idet], 'arrow_length_ratio': 0.1, 'alpha': 0.25}
            for i, arg in enumerate(('X', 'Y', 'Z')):
                k_quiver_kwargs[arg] = self.sample.position[i]
            for i, arg in enumerate(('U', 'V', 'W')):
                k_quiver_kwargs[arg] = self.get_Ks()[idet][i]*np.linalg.norm(k_det)
            self.lab_K_vecs.append(self.lab_ax.quiver(**k_quiver_kwargs))
        return lab_K_vecs







    def plot_all(self):
        self.update()
        for surf in self.lab_sample:
            surf.remove()
        self.calc_ewald()
        self.show_ewald()
        self.recip_sample = self.show_sample(ax = self.ax)
        self.lab_sample = self.show_sample(ax=self.lab_ax, ratio=self.sample.sample_scale)
        self.show_pole()
        self.plot_pole_figure(self.View.pole_proj_ax, )
        self.plot_pole_figure(chosen_ax=self.calc_pf_ax, use_intensities=True)