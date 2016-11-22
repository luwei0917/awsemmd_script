#!/usr/bin/env python3
from pandas import DataFrame
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
#
# x = ['A']*300 + ['B']*400 + ['C']*300
# y = np.random.randn(1000)
# df = DataFrame({'Letter': x, 'N': y})
# grouped = df.groupby('Letter')
#
# for group in grouped:
#     plt.figure()
#     plt.hist(group[1].N)
#     plt.show()
gb_seq = '/Users/weilu/Research/server/nov14/ga_200/simulation/rerun_0/energy.dat'
ga_seq = '/Users/weilu/Research/server/nov14/ga_200/simulation/0/energy.dat'

# ga_seq = '/Users/weilu/Research/server/nov14/gb_200/simulation/rerun_0/energy.dat'
# gb_seq = '/Users/weilu/Research/server/nov14/gb_200/simulation/0/energy.dat'
df = pd.read_csv(ga_seq, sep='\s+')
# print(df)
df = df.drop('Step', 1)
df.drop('Shake', axis=1, inplace=True)
df.drop('Excluded', axis=1, inplace=True)
df.drop('AMH-Go', axis=1, inplace=True)
df.drop('Vec_FM', axis=1, inplace=True)
df.drop('SSB', axis=1, inplace=True)
df.drop('Membrane', axis=1, inplace=True)
df.drop('Electro.', axis=1, inplace=True)
matplotlib.style.use('ggplot')
# plt.style.use('presentation')
# print(df.describe())
# color = dict(boxes='DarkGreen', whiskers='DarkOrange', medians='DarkBlue', caps='Gray')
# color = dict(medians='Red')
# medianprops = dict(linestyle='-.', linewidth=2.5, color='firebrick')
# ax = df.boxplot(return_type='axes')
ax = df.hist()
# , patch_artist=True, medianprops=medianprops
# for subax in ax:
#     for median in subax['medians']:
#         median.set(color='Red')
# plt.setp(ax['medians'], color='red')
# for median in ax['medians']:
#     median.set(color='Red')

# ax.set_ylim(-300, 200)
# ax.set_ylim(-200, 100)
# ax.set_ylim(-100, 100)
# ax.set_ylim(-280, -180)
# print(matplotlib.style.available)
matplotlib.style.use('classic')
# /scratch/wl45/nov14/gb_200/simulation/rerun_0

# df_gb = pd.read_csv(gb_seq, sep='\s+')
# df_gb = df_gb.drop('Step', 1)
# df_gb.drop('Shake', axis=1, inplace=True)
# df_gb.drop('Excluded', axis=1, inplace=True)
# df_gb.drop('AMH-Go', axis=1, inplace=True)
# df_gb.drop('Vec_FM', axis=1, inplace=True)
# df_gb.drop('SSB', axis=1, inplace=True)
# df_gb.drop('Membrane', axis=1, inplace=True)
# df_gb.drop('Electro.', axis=1, inplace=True)
# ax = df_gb.boxplot(ax=ax, return_type='axes')
# ax.set_title("red is ga seq, blue is gb seq")
fig = ax.get_figure()
fig.savefig('figure.pdf')
# import random
# import numpy
# from matplotlib import pyplot
#
# x = [random.gauss(3,1) for _ in range(400)]
# y = [random.gauss(4,2) for _ in range(400)]
#
# bins = numpy.linspace(-10, 10, 100)
#
# pyplot.hist(x, bins, alpha=0.5, label='x')
# pyplot.hist(y, bins, alpha=0.5, label='y')
# pyplot.legend(loc='upper right')
# pyplot.show()


# df =pd.DataFrame({'col1':np.random.randn(100),'col2':np.random.randn(100)})
# df.hist(layout=(1,2))
# df.plot()
