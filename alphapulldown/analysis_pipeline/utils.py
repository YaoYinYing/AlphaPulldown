import pathlib
from configparser import Interpolation
import IPython.display as display
import ipywidgets as widgets
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import os
#from analysis_pipeline.af2hyde_mod import plot_predicted_alignment_error
from af2plots.plotter import plotter


def display_pae_plots(subdir,figsize=(50, 50),top=5):
    """A function to display all the pae plots in the subdir"""
    subdir=pathlib.Path(subdir).resolve()
    job=subdir.stem
    images = sorted([i for i in os.listdir(subdir) if ".png" in i])
    # expected image file: <job>_PAE_plot_ranked_<rank_id>.png
    assert top <= len(images)
    if len(images) > 0:
        fig, axs = plt.subplots(1, top, figsize=figsize)
        for rank in range(top):
            img = plt.imread(os.path.join(subdir, f'{job}_PAE_plot_ranked_{rank}.png'))
            axs[rank].imshow(img,interpolation="nearest")
            axs[rank].axis("off")
        #plt.show()
    else:
        #plot_predicted_alignment_error(subdir)
        af2o = plotter()
        dd = af2o.parse_model_pickles(subdir)
        ff = af2o.plot_predicted_alignment_error(dd)

    plt.show()
