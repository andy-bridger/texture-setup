Some python code for dealing with Pole Figures for texture analysis.

Two main functionalities: Simulating an experiment `mantex_sim.py`; and investigating an experiment `NyRTex/run_NyRTex.py`

Both of these use the same underlying code base, scattering around the place but the primary logic for pole figures is in `model.py`, with a wrapper for both the simulation `mantex_vis.py` and the experiment visualisation `mantex_run.py`.

Simulating an experiment:

![Screenshot 2025-01-20 121502](https://github.com/user-attachments/assets/42880144-f623-4e9b-babf-81211100a175)

This allows you to set sample geometry with an STL file, set up detector, source and sample positions. You can assign one or more underlying crystal lattice to the sample and visualise how the ewald sphere and the pole figure projection interact to help understand what is happening within the experiment. You can either move the sample with the goniometer widget to see how the system changes or pre-load some angles and just run the experiment. Support for how changing the probe position in Q space alters the experiment practically and conceptually 

Visualising an experiment:

![Screenshot 2025-01-30 155712](https://github.com/user-attachments/assets/ba52ec48-96a3-45bf-89bc-6ec265887838)

This allows you to load in Sample Geometry (.STL), Sample Positions (as either a series of Rotation Matrix + Translation Vector or as goniometer angles that get converted under the hood), Apparatus Geometry (Source/ detector positions), Detector Binnings (list of boolean masks of which detectors should get binned together) and the detector signals (as an array of intensities and an array of q-space sample points) and it will let you visualise the setup of each run of the experiment and show how that has contributed to the final pole figure.

Next step is to implement this in `mantid`
