import os
import sys

import plotly.graph_objs as go


path_i = os.path.join(
    # Change this to the path of the file on your computer
    # "/mnt/f/Dropbox/__temp__/00_bands_data/band_disp_2.0.pickle",
    # "/mnt/f/Dropbox/__temp__/00_bands_data/band_disp_2.5.pickle",
    # "/mnt/f/Dropbox/__temp__/00_bands_data/band_disp_3.0.pickle",
    # "/mnt/f/Dropbox/__temp__/00_bands_data/band_disp_3.5.pickle",
    # "/mnt/f/Dropbox/__temp__/00_bands_data/band_disp_4.0.pickle",
    # "/mnt/f/Dropbox/__temp__/00_bands_data/band_disp_4.5.pickle",
    )

layout = go.Layout(
    # xaxis=dict(),

    yaxis=dict(

	# Y-axis range defined here
	range=[-4, 4],
	),
    )
