from PIL.ImageOps import scale


class Presenter():
    def __init__(self, Model, View):
        self.Model = Model
        self.View = View
        self.init_view()
    def setup_view(self):
        self.View.setup_view(self.Model.pole_radius)
    def toggle_visible(self, artist):
        artist.set_visible(not artist._visible)
    def plot_ewald(self):
        sphere_data = self.Model.calc_ewald(self.Model.ewald_radius)
        self.ewald_sphere_artist = self.View.plot_sphere(self.View.recip_ax, sphere_data)
    def plot_recip_lattices(self):
        self.r_latt_artists = self.View.show_reciprocal_lattices(self.View.recip_ax,
                                                                 self.Model.sample.lab_space_rlatts,
                                                                 self.Model.sample.cell_colors,
                                                                 self.Model.sample.rl_alphas)
    def plot_recip_sample(self):
        self.recip_sample_artist = self.View.plot_cube(self.View.recip_ax, self.Model.get_sample(), ('red', 'blue', 'green'))
    def plot_lab_sample(self):
        self.lab_sample_artist = self.View.plot_cube(self.View.lab_ax, self.Model.get_sample(ratio=1), ('red', 'blue', 'green'))
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
                                                           scale = self.Model.ki_raw_scale)
    def plot_lab_components(self):
        source_repr, source_col, det_reprs, det_cols = self.Model.get_lab_frame()
        self.source_artist = self.View.plot_cube(self.View.lab_ax, source_repr, source_col)
        self.det_artists = [self.View.plot_cube(self.View.lab_ax,
                                                det_repr,
                                                det_cols[i]) for i, det_repr in enumerate(det_reprs)]
    def plot_pole_sphere(self):
        self.pole_sphere_artist = self.View.plot_sphere(self.View.recip_ax, self.Model.pole_cart_array)
        self.pole_sphere_eq_artist = self.View.plot_line(self.View.recip_ax, self.Model.eq_cart)
        self.pole_sphere_vector_artists = self.View.plot_pole_figure_vectors(self.View.recip_ax,
                                                                             self.Model.pole_Ks,
                                                                             self.Model.detector_colors,
                                                                             self.Model.pole_K_projs_vec)
    def plot_pole_figure_detectors(self):
        self.pole_figure_detector_artists = self.View.plot_pole_figure_position(self.View.pole_proj_ax,
                                                                                self.Model.sample_view_axis,
                                                                                self.Model.pole_figure_points,
                                                                                self.Model.detector_colors,
                                                                                self.Model.equator())
    def plot_pole_figure_intensities(self):
        self.pole_figure_artists = self.View.plot_pole_figure_intensities(self.View.calc_pf_ax,
                                                                                self.Model.sample_view_axis,
                                                                                self.Model.pole_figure_intensities,
                                                                                self.Model.equator())
    def plot_all(self):
        # Lab Frame
        self.plot_lab_sample()
        self.plot_lab_beam_paths()
        self.plot_lab_Ks()
        self.plot_lab_components()

        # Recip Frame
        self.plot_recip_sample()
        self.plot_recip_lattices()
        self.plot_ewald()
        self.plot_pole_sphere()

        # Pole Detector Frame
        self.plot_pole_figure_detectors()

        # Pole Figure Frame
        self.plot_pole_figure_intensities()

        self.View.fix_aspect()
    def init_view(self):
        self.View.setup_view(self.Model.pole_radius)
        self.View.update_view_axes()
