from PIL.ImageOps import scale
from helper_funcs import *


class Presenter():
    def __init__(self, Model, View):
        self.Model = Model
        self.View = View
        self.init_view()
        self.pole_figure_artists =[]
    def setup_view(self):
        self.View.setup_view()
    def toggle_visible(self, artist):
        artist.set_visible(not artist._visible)
    def plot_probe_ewald(self):
        probe_sphere_data = self.Model.calc_ewald(self.Model.q_probe)
        self.p_ewald_sphere_artist = self.View.plot_sphere(self.View.recip_ax, probe_sphere_data, c='gold')
    def plot_ewald(self):
        inner_sphere_data = self.Model.calc_ewald(self.Model.ewald_radii[0])
        outer_sphere_data = self.Model.calc_ewald(self.Model.ewald_radii[1])
        self.i_ewald_sphere_artist = self.View.plot_sphere(self.View.recip_ax, inner_sphere_data, alpha = 0.1)
        self.o_ewald_sphere_artist = self.View.plot_sphere(self.View.recip_ax, outer_sphere_data, alpha = 0.1)
        self.plot_probe_ewald()
    def plot_recip_lattices(self):
        self.Model.sample.adjust_alphas(self.Model.qs_of_interest[self.Model.probe_ind])
        self.r_latt_artists = self.View.show_reciprocal_lattices(self.View.recip_ax,
                                                                 self.Model.sample.lab_space_rlatts,
                                                                 self.Model.sample.cell_colors,
                                                                 self.Model.sample.rl_alphas)
    def plot_recip_sample(self):
        if self.Model.sample.primitive:
            self.recip_sample_artist = self.View.plot_cube(self.View.recip_ax,
                                                           self.Model.get_sample(ratio = 0.05*self.Model.r_dict['recip_sample']),
                                                           ('red', 'blue', 'green'))
        else:
            self.recip_sample_artist = self.View.plot_mesh(self.View.recip_ax, self.Model.sample.get_mesh_vectors(scale = 0.1))

    def plot_ewald_detector_Ks(self):
        self.ewald_detector_vector_artists = self.View.plot_ewald_detector_vectors(self.View.recip_ax,
                                                                             self.Model.pKs[:,self.Model.probe_ind],
                                                                             self.Model.detector_colors)
    def plot_lab_sample(self):
        if self.Model.sample.primitive:
            self.lab_sample_artist = self.View.plot_cube(self.View.lab_ax,
                                                         self.Model.get_sample(ratio=2*self.Model.r_dict['lab_sample']),
                                                         ('red', 'blue', 'green'))
        else:
            # only upon plotting do we want to recalc the mesh
            self.lab_sample_artist = self.View.plot_mesh(self.View.lab_ax, self.Model.sample.get_mesh_vectors(self.Model.r_dict['lab_sample']))

    def update_lab_sample(self):
        self.remove_artist_set(self.lab_sample_artist)
        self.plot_lab_sample()

    def plot_lab_beam_paths(self):
        self.det_path_artists = self.View.plot_lab_detector_paths(self.View.lab_ax, self.Model.sample.position,
                                                                  self.Model.get_detector_vectors())
        self.source_path_artists = self.View.plot_lab_source_path(self.View.lab_ax, self.Model.source.position,
                                                                  self.Model.ki_raw)
    def plot_lab_Ks(self):
        self.lab_K_vec_artists = self.View.plot_lab_K_vecs(self.View.lab_ax,
                                                           self.Model.sample.position,
                                                           self.Model.Ks,
                                                           self.Model.detector_colors,
                                                           scale = self.Model.ki_raw_scale*self.Model.r_dict['k_vecs']/4)
    def plot_lab_components(self):
        source_repr, source_col, det_reprs, det_cols = self.Model.get_lab_frame()
        self.source_artist = self.View.plot_cube(self.View.lab_ax, source_repr, source_col)
        self.det_artists = [self.View.plot_cube(self.View.lab_ax,
                                                det_repr,
                                                det_cols[i]) for i, det_repr in enumerate(det_reprs)]
    def plot_pole_sphere(self):
        self.pole_sphere_artist = self.View.plot_sphere(self.View.recip_ax, self.Model.pole_cart_array)
        self.pole_sphere_eq_artist = self.View.plot_line(self.View.recip_ax,
                                                         self.Model.pole_view_mat@equator(r=1))
        self.pole_sphere_vector_artists = self.View.plot_pole_figure_vectors(self.View.recip_ax,
                                                                             self.Model.pole_Ks,
                                                                             self.Model.detector_colors,
                                                                             self.Model.pole_K_projs_vec)
    def plot_pole_figure_detectors(self):
        self.pole_figure_detector_artists = self.View.plot_pole_figure_position(self.View.pole_proj_ax,
                                                                                self.Model.pole_figure_points,
                                                                                self.Model.detector_colors,
                                                                                equator())
    def plot_pole_figure_intensities(self):
        self.pole_figure_artists += self.View.plot_pole_figure_intensities(self.View.calc_pf_ax,
                                                                                self.Model.pole_figure_intensities,
                                                                                equator())
    def plot_goniometer(self):
        self.z_artist = [self.View.plot_line(self.View.lab_ax, self.Model.goniometer.ax1_eq, self.Model.r_dict['gonio_r']),
                         self.View.plot_line(self.View.lab_ax, self.Model.goniometer.ax1_eq[:,:self.Model.goniometer.phi_frac],
                                             self.Model.r_dict['gonio_r'], 'red'),
                         self.View.plot_goniometer_axis(self.View.lab_ax, 'Z',
                                                        self.Model.sample.position,
                                                        self.Model.goniometer.ax1,
                                                        self.Model.r_dict['gonio_v'])]
        self.xp_artist = [self.View.plot_line(self.View.lab_ax, self.Model.goniometer.ax2_eq, self.Model.r_dict['gonio_r']),
                          self.View.plot_line(self.View.lab_ax,
                                              self.Model.goniometer.ax2_eq[:, :self.Model.goniometer.theta_frac],
                                              self.Model.r_dict['gonio_r'], 'blue'),
                         self.View.plot_goniometer_axis(self.View.lab_ax, "X'",
                                                        self.Model.sample.position,
                                                        self.Model.goniometer.ax2,
                                                        self.Model.r_dict['gonio_v'])]
        self.zp_artist = [self.View.plot_line(self.View.lab_ax, self.Model.goniometer.ax3_eq, self.Model.r_dict['gonio_r']),
                         self.View.plot_line(self.View.lab_ax, self.Model.goniometer.ax3_eq[:, :self.Model.goniometer.psi_frac],
                                self.Model.r_dict['gonio_r'], 'green'),
                         self.View.plot_goniometer_axis(self.View.lab_ax, "Z'",
                                                        self.Model.sample.position,
                                                        self.Model.goniometer.ax3,
                                                        self.Model.r_dict['gonio_v'])]

    def plot_readouts(self):
        self.det_readout_plot_artist = self.View.show_detector_readout(self.View.det_ax,
                                                                       self.Model.detector_readout,
                                                                       self.Model.q_range,
                                                                       self.Model.detector_colors)

    def plot_q_probe(self):
        self.det_probe_artist = self.View.show_detector_probe(self.View.det_ax,
                                                              self.Model.q_probe,
                                                              self.Model.probe_window)

    def plot_all(self):
        # Lab Frame
        self.plot_lab_sample()
        self.plot_lab_beam_paths()
        self.plot_lab_Ks()
        self.plot_lab_components()
        self.plot_goniometer()
        self.plot_readouts()

        # Recip Frame
        self.plot_recip_sample()
        self.plot_recip_lattices()
        self.plot_ewald()
        self.plot_ewald_detector_Ks()
        self.plot_pole_sphere()

        # Pole Detector Frame
        self.plot_pole_figure_detectors()
        self.plot_q_probe()

        # Pole Figure Frame
        self.plot_pole_figure_intensities()

        self.View.fix_aspect()

    def run_experiment(self):
        for gonio_pos in self.Model.goniometer.exp_runs:
            self.upon_goniometer_change(*gonio_pos)
            self.Model.exp_data.state_data.append(self.Model.goniometer.to_state_data())
            self.Model.exp_data.detector_readouts.append([det.spectrum for det in self.Model.detectors])
        self.Model.exp_data.save_exp()

    def assign_file_name(self, val):
        self.Model.create_experiment()
        self.Model.exp_data.name = val

    def click_run_exp(self, val):
        self.remove_artist_set(self.pole_figure_artists)
        self.pole_figure_artists = []
        self.Model.reset_pole_figure()
        self.run_experiment()
        self.plot_pole_figure_intensities()
        self.View.fix_aspect()
        self.View.fig.canvas.draw()

    def init_view(self):
        self.View.setup_view()
        self.View.update_view_axes()
        self.View.add_goniometer_widgets(self.goniometer_update)
        self.View.add_lab_k_widgets(self.lab_k_update)
        self.View.add_sample_scale_widget(self.sample_scale_update)
        self.View.add_gonio_ring_scale_widget(self.gr_scale_update)
        self.View.add_gonio_axis_scale_widget(self.ga_scale_update)
        self.View.add_probe_pos_widget(self.update_probe_pos,
                                       self.Model.ewald_radii[0],
                                       self.Model.ewald_radii[1],
                                       self.Model.q_probe)
        self.View.add_run_experiment_button(self.click_run_exp)
        self.View.add_exp_name_textbox(self.assign_file_name)

    def upon_goniometer_change(self, phi, theta, psi):
        self.Model.goniometer.update_angles(phi, theta, psi)
        self.Model.sample.update()
        self.Model.update()

    def goniometer_update(self, val):
        self.upon_goniometer_change(self.View.slider_phi.val, self.View.slider_theta.val, self.View.slider_psi.val)
        self.View.update_view_axes()

        # Lab Frame
        self.update_lab_sample()
        self.remove_artist_set(self.z_artist)
        self.remove_artist_set(self.xp_artist)
        self.remove_artist_set(self.zp_artist)
        self.Model.goniometer.update_norms()
        self.Model.goniometer.update_equators()
        self.Model.goniometer.update_fracs()
        self.plot_goniometer()

        # Recip Frame
        self.plot_recip_sample()
        self.plot_recip_lattices()
        self.plot_ewald()
        self.plot_pole_sphere()
        self.plot_ewald_detector_Ks()

        # Pole Detector Frame
        self.plot_pole_figure_detectors()

        # Detector Readout
        self.View.det_ax.clear()
        self.plot_readouts()
        self.plot_q_probe()

        # Pole Figure Frame
        self.plot_pole_figure_intensities()


        self.View.fix_aspect()

    def lab_k_update(self, val):
        self.remove_artist_set(self.lab_K_vec_artists)
        self.Model.r_dict['k_vecs'] = val
        self.plot_lab_Ks()

    def sample_scale_update(self, val):
        self.remove_artist_set(self.lab_sample_artist)
        self.remove_artist_set(self.recip_sample_artist)
        self.Model.r_dict['recip_sample'] = val
        self.Model.r_dict['lab_sample'] = val
        if not self.Model.sample.primitive:
            self.Model.sample.scale = val
        self.plot_lab_sample()
        self.plot_recip_sample()
        self.View.fix_aspect()

    def gr_scale_update(self, val):
        self.remove_artist_set(self.z_artist)
        self.remove_artist_set(self.xp_artist)
        self.remove_artist_set(self.zp_artist)
        self.Model.r_dict['gonio_r'] = val
        self.plot_goniometer()
        self.View.fix_aspect()

    def ga_scale_update(self, val):
        self.remove_artist_set(self.z_artist)
        self.remove_artist_set(self.xp_artist)
        self.remove_artist_set(self.zp_artist)
        self.Model.r_dict['gonio_v'] = val
        self.plot_goniometer()
        self.View.fix_aspect()

    def update_probe_pos(self, val):
        self.remove_artist_set(self.det_probe_artist)
        self.remove_artist_set(self.pole_figure_artists)
        self.remove_artist_set(self.r_latt_artists)
        self.remove_artist_set(self.ewald_detector_vector_artists)
        self.p_ewald_sphere_artist.remove()
        self.pole_figure_artists = []
        self.Model.update_q_probe(val)
        self.Model.reset_pole_figure()
        self.Model.get_inplane_pole_figure()
        self.plot_pole_figure_intensities()
        self.plot_q_probe()
        self.plot_probe_ewald()
        self.plot_recip_lattices()
        self.plot_ewald_detector_Ks()

    def remove_artist_set(self, artists):
        if type(artists) != list:
            artists.remove()
        else:
            for a in artists:
                if type(a) != list:
                    a.remove()
                else:
                    self.remove_artist_set(a)


