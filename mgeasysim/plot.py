import warnings
from mgeasysim import config as cf
import pandas as pd
from rapidfuzz import process, fuzz
import numpy as np
import subprocess
import os
import glob
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns


def mash_plot(matches):
    fig, ax = plt.subplots(figsize = (6,4))
    sns.histplot(matches['alt_ani'])
    sns.despine()
    ax.set_title('Mash of main matches to alts')
    plt.savefig(os.path.join(cf.OUTPUT, 'mash.png'), dpi=400, bbox_inches='tight')

def abundance_plot(abundances):
    fig = plt.figure(figsize=(8,4))
    ax = fig.add_subplot(1,2,1)
    fig.suptitle(f'{num_species} species example distribution')
    sns.histplot(abundances, bins=50)
    sns.despine()

    ax = fig.add_subplot(1,2,2)
    sns.histplot(abundances, bins=np.logspace(-6, 0, 50))
    ax.set_xscale('log')
    sns.despine()

    plt.savefig(os.path.join(cf.OUTPUT, 'abundances.png'), dpi=400, bbox_inches='tight')
