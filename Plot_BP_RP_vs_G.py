import glob
import os
from laspec import mrs

from astropy.table import Table
import astropy.units as u
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from laspec.mrs import MrsSpec
import random

import seaborn as sns

import dustmaps
from dustmaps.config import config
config['data_dir'] = "/Users/liujunhui/Desktop/2021workMac/202111012totallynew/dustmaps"

import dustmaps.bayestar

dustmaps.bayestar.fetch()

from dustmaps.bayestar import BayestarQuery

bayestar = BayestarQuery(version='bayestar2019')

from astropy.coordinates import SkyCoord, Distance
from astropy import units


dir_name = './20211205_GaiaEDR3_3arcsec.csv'
m8pd = pd.read_csv(dir_name)

c = SkyCoord(m8pd['combined_ra_obs'], m8pd['combined_dec_obs'], unit = 'deg')
d = 1000/m8pd['parallax'].values * units.pc

d[d<0] =0

coords = SkyCoord(c.galactic.l, c.galactic.b, distance = d, frame='galactic')

ebv = bayestar(coords, mode = 'median')

"""
Filter     G          G_BP      G_RP      B_T     V_T     J            H         Ks
?eff(?) 6437.70     5309.57   7709.85  4265.42  5332.38  12329.79   16395.59   21522.05
A?/AV   0.85926     1.06794   0.65199  1.35552  1.05047  0.29434    0.18128    0.11838
"""

AG = ebv*0.85926
ABP = ebv*1.06794
ARP = ebv*0.65199
AJ = ebv*0.29434
AH = ebv*0.18128
AK = ebv*0.11838

#%% absolute magnitude
dm = 5*np.log10(1000/m8pd['parallax']) - 5
MG = m8pd['phot_g_mean_mag'] - dm -AG
MBP = m8pd['phot_bp_mean_mag'] - dm -ABP
MRP = m8pd['phot_rp_mean_mag'] - dm -ARP

MBP_RP = MBP - MRP


map = sns.kdeplot(MBP_RP, MG)
plt.plot(MBP_RP, MG, '.')
plt.show()



















