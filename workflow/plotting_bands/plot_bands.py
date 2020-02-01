# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.2
#   kernelspec:
#     display_name: Python [conda env:py27]
#     language: python
#     name: conda-env-py27-py
# ---

import plotly

# #############################################################################
import pickle; import os
path_i = os.path.join(
    # os.environ[""],
    "/home/raulf2012/Dropbox/band_disp.pickle",
    )
    # "/home/raulf2012/__temp__/band_disp.pickle")
with open(path_i, "rb") as fle:
    bands_data = pickle.load(fle)
# #############################################################################

# +
from dft_post_analysis.bands import plot_bands

traces, tmp = plot_bands(bands_data, plot_title="dijis")

# +
import chart_studio.plotly as py
import plotly.graph_objs as go


fig = go.Figure(data=traces)
fig.show()
# -

from plotly import io as pyio

# +
# from plotting
# import plotting.my_plotly
# from plotting.my_plotly import my_plotly_plot

pyio.write_html(
    fig,
    "outplot.html",
    # os.path.join(plot_dir, plot_name + ".html"),
    # config=None,
    # auto_play=True,
    # include_plotlyjs=True,
    # include_mathjax=False,
    # post_script=None,
    # full_html=True,
    # animation_opts=None,
    # validate=True,
    # default_width='100%',
    # default_height='100%',
    # auto_open=False,
    )

# -

from ase_modules

traces
