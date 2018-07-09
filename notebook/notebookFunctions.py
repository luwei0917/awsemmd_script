import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy.interpolate import interp2d
import sys
import argparse
import os
import numpy as np
from numpy.random import uniform
import pandas as pd
from itertools import product
import datetime
import glob
import re
from numpy import ma
import networkx as nx
import scipy
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
# plt.rcParams['figure.figsize'] = (10,5)
# plt.rcParams['figure.figsize'] = (10,5.625)   # 16:9
plt.rcParams['figure.figsize'] = (10,6.180)    # golden ratio
# plt.rcParams['figure.figsize'] = (10*2,6.180*2)    #golden ratio

def temperature_exchange_table(data):
    b = data.groupby(["Run", "Temp"])["Step"].count().reset_index()
    c = b.pivot(index="Run", columns="Temp", values="Step").reset_index(drop=True)
    return c

def summarise_temperature_exchange_table(data):
    tmp = temperature_exchange_table(data)
    return (1 - np.isnan(tmp)).sum(axis=1)

def select(t, i=100):
    return t.groupby(["BiasTo", "Run"])["DisReal"].describe().query(f"count > {i}")

def raw_2d_plot(location):
    data = np.loadtxt(location)
    x = data[:,1]
    y = data[:,2]
    z = data[:,3]
    plt.scatter(x,y, c=z, cmap='rainbow')
    # plt.gray()
    plt.colorbar()
    return data

def show_images_all(all_data, temp=450, zmax=20, xlabel="xlabel", ylabel="ylabel", mode="2d_z_qw", force=0.2):
    plt.close('all')
    nrows = 4
    fig, ax = plt.subplots(nrows=nrows,ncols=3, figsize=(10*3,6.180*nrows), dpi=200)
    x = 1
    y = 2
    z = 3
    zmin = 0
    zmax=20
    titlefontsize = 28
    test = all_data.query(f"mode == '{mode}'").query(f"temp == '{temp}'").query(f"force == '{force}'")
    origin = test.query("perturbation == 'original'")
    one_change_data = test.query("change == 'rg'")
    up = one_change_data.query("upOrDown == 'p'")
    down = one_change_data.query("upOrDown == 'm'")
    plot_data = [down, origin, up]
    dic_upOrDown = {0:"down", 1:"origin", 2:"up"}
    test = test.query("change != 'none'")
    for idx, (change, one_change_data) in enumerate(test.groupby("change")):
        up = one_change_data.query("upOrDown == 'p'")
        down = one_change_data.query("upOrDown == 'm'")
        plot_data = [down, origin, up]
        for image_idx, pddata in enumerate(plot_data):
            data = pddata[["index", "x","y","f"]].values
    #         print(data)
            data = data[~np.isnan(data).any(axis=1)]  # remove rows with nan
            data = data[~(data[:,z] > zmax)]  # remove rows of data for z not in [zmin zmax]
            data = data[~(data[:,z] < zmin)]

            xi = np.linspace(min(data[:,x]), max(data[:,x]), 20)
            yi = np.linspace(min(data[:,y]), max(data[:,y]), 20)
            zi = griddata((data[:,x], data[:,y]), data[:,z], (xi[None,:], yi[:,None]), method='linear')
            # plt.contour(xi, yi, zi, 50, linewidths=0.25,colors='k')
            jet = cm = plt.get_cmap('jet')
    #         print(jet)
            # plt.contourf(xi, yi, zi, 20, cmap='rainbow')
    #         plt.figure()
            cs = ax[idx,image_idx].contourf(xi, yi, zi, 30, cmap='jet')
            # plt.xlim(xmin, xmax)
            cs.set_clim(zmin, zmax)
    #         fig.clim()
            ax[idx,image_idx].set_title(f"{change}: {dic_upOrDown[image_idx]}", fontsize=titlefontsize)
        fig.colorbar(cs, ax=ax[idx,2], shrink=1)
#     fig.colorbar()
    fig.suptitle(f'temp = {temp}', y=1.02, fontsize=titlefontsize*1.5)
#     fig.subplots_adjust(top=1.02)
    fig.tight_layout()

def getBound(location, res=30, zmin=0, zmax=20, x=1, y=2, z=3):
    # data = np.where(np.isnan(data), zmax, data)
    data = np.loadtxt(location)
    data = data[~np.isnan(data).any(axis=1)]  # remove rows with nan
    if zmin == -1:
        zmin = data[:,3].min()
    if zmax == -1:
        zmax = data[:,3].max()
    data = data[~(data[:,z] > zmax)]  # remove rows of data for z not in [zmin zmax]
    data = data[~(data[:,z] < zmin)]
    xmin = min(data[:,x])
    xmax = max(data[:,x])
    ymin = min(data[:,y])
    ymax = max(data[:,y])
    return(xmin,xmax,ymin,ymax)

def getxyz(data, res=30, zmin=0, zmax=20, x=1, y=2, z=3, xmin=-1,xmax=-1,ymin=-1,ymax=-1):
    # data = np.where(np.isnan(data), zmax, data)
    data = data[~np.isnan(data).any(axis=1)]  # remove rows with nan
    if zmin == -1:
        zmin = data[:,3].min()
    if zmax == -1:
        zmax = data[:,3].max()
    data = data[~(data[:,z] > zmax)]  # remove rows of data for z not in [zmin zmax]
    data = data[~(data[:,z] < zmin)]
    if xmin == -1:
        xi = np.linspace(min(data[:,x]), max(data[:,x]), res)
    else:
        xi = np.linspace(xmin, xmax, res)
    if ymin == -1:
        yi = np.linspace(min(data[:,y]), max(data[:,y]), res)
    else:
        yi = np.linspace(ymin, ymax, res)
    zi = griddata((data[:,x], data[:,y]), data[:,z], (xi[None,:], yi[:,None]), method='linear')
    return (xi,yi,zi)

def getxyz_2(data, res=30, zmin=0, zmax=20, x=1, y=2, z=3):
    # data = np.where(np.isnan(data), zmax, data)
    data = data[~np.isnan(data).any(axis=1)]  # remove rows with nan
    if zmin == -1:
        zmin = data[:,3].min()
    if zmax == -1:
        zmax = data[:,3].max()
    data = data[~(data[:,z] > zmax)]  # remove rows of data for z not in [zmin zmax]
    data = data[~(data[:,z] < zmin)]

    xi = np.linspace(min(data[:,x]), max(data[:,x]), res)
    yi = np.linspace(min(data[:,y]), max(data[:,y]), res)

    # fill in those nan with zmax
    tmp = np.ones((res**2, 4))*zmax
    pos = 0
    count = 0
    for i in range(res):
        for j in range(res):
            tmp[pos,0] = pos
            tmp[pos,1] = xi[i]
            tmp[pos,2] = yi[j]
            if count < data.shape[0] and pos == int(data[count, 0]):
                tmp[pos,3] = data[count, z]
                count += 1
            pos += 1
    zi = griddata((tmp[:,x], tmp[:,y]), tmp[:,z], (xi[None,:], yi[:,None]), method='linear')
    return (xi,yi,zi)

def getxyz_3(data, res=30, zmin=0, zmax=20, x=1, y=2, z=3):
    # data = np.where(np.isnan(data), zmax, data)
    data = data[~np.isnan(data).any(axis=1)]  # remove rows with nan
    if zmin == -1:
        zmin = data[:,3].min()
    if zmax == -1:
        zmax = data[:,3].max()
    data = data[~(data[:,z] > zmax)]  # remove rows of data for z not in [zmin zmax]
    data = data[~(data[:,z] < zmin)]

    xi = np.linspace(min(data[:,x]), max(data[:,x]), res)
    yi = np.linspace(min(data[:,y]), max(data[:,y]), res)

    # construct complete res*res array
    tmp = np.ones((res**2, 4))*np.nan
    pos = 0
    count = 0
    for i in range(res):
        for j in range(res):
            tmp[pos,0] = pos
            tmp[pos,1] = xi[i]
            tmp[pos,2] = yi[j]
            if count < data.shape[0] and pos == int(data[count, 0]):
                tmp[pos,3] = data[count, z]
                count += 1
            pos += 1

    # Assign zmax around for each point if it is nan before.
    new_tmp = tmp.copy()
    pos = 0
    count = 0
    for i in range(res):
        for j in range(res):
            if not np.isnan(tmp[pos,3]):
                new_tmp[pos] = tmp[pos]
                right_neighbour = i*res + (j+1)
                if j < res-1 and np.isnan(tmp[right_neighbour,3]):
                    new_tmp[right_neighbour, 3] = zmax
                left_neighbour = i*res + (j-1)
                if j>0 and np.isnan(tmp[left_neighbour,3]):
                    new_tmp[left_neighbour, 3] = zmax
                up_neighbour = (i+1)*res + j
                if i < res-1 and np.isnan(tmp[up_neighbour,3]):
                    new_tmp[up_neighbour, 3] = zmax
                down_neighbour = (i-1)*res + j
                if i > 0 and np.isnan(tmp[down_neighbour,3]):
                    new_tmp[down_neighbour, 3] = zmax
            pos += 1
    tmp = new_tmp[~np.isnan(new_tmp).any(axis=1)]
    zi = griddata((tmp[:,x], tmp[:,y]), tmp[:,z], (xi[None,:], yi[:,None]), method='linear')
    return (xi,yi,zi)

def plot2d(location, path, temp="450", res=30, zmin=0, zmax=30, z=3, xlabel="xlabel", ylabel="ylabel", title="", outname=None, **kargs):
    titlefontsize = 28
    data = np.loadtxt(location)
    xi, yi, zi = getxyz(data, zmin=zmin, zmax=zmax,res=res, z=z, **kargs)
    # V = ma.masked_array(zi, zi>40)
    # zi = np.where(np.isnan(zi), 1e6, zi)
    f_on_path = [zi[tuple(p)] for p in reversed(path)]
    # plt.contour(xi, yi, zi, 50, linewidths=0.25,colors='k')
    jet = cm = plt.get_cmap('jet')
    print(jet)
    # plt.contourf(xi, yi, zi, 20, cmap='rainbow')
    plt.figure()
    plt.contourf(xi, yi, zi, res, cmap='jet')
    plt.plot(xi[path[:,1]], yi[path[:,0]], 'r.-')

    # plt.xlim(xmin, xmax)
    # plt.clim(zmin, zmax)
    plt.colorbar()

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title, y=1.02, fontsize=titlefontsize)
    # plt.tight_layout()
    # plt.axis('equal')
    # plt.axes().set_aspect('equal')
    # plt.axes().set_aspect('scaled')
    if outname:
        plt.savefig(outname, dpi=300, bbox_inches='tight')
    # plt.show()
    # return (xi,yi,zi)
    return f_on_path

def plot2d_side_by_side(location1, location2, zmin=0, zmax=30, xlabel="xlabel", ylabel="ylabel", title1="", title2="", outname=None):
    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(20, 6.18))
    titlefontsize = 28
    # plt.contour(xi, yi, zi, 50, linewidths=0.25,colors='k')
    jet = cm = plt.get_cmap('jet')
    # print(jet)
    # plt.contourf(xi, yi, zi, 20, cmap='rainbow')
    data = np.loadtxt(location1)
    xi, yi, zi = getxyz(data, zmin=zmin, zmax=zmax)
    g = ax1.contourf(xi, yi, zi, 30, cmap='jet')
    # plt.xlim(xmin, xmax)
    g.set_clim(zmin, zmax)
    # plt.colorbar()

    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    ax1.set_title(title1, y=1.02, fontsize=titlefontsize)
    data = np.loadtxt(location2)
    xi, yi, zi = getxyz(data, zmin=zmin, zmax=zmax)
    g = ax2.contourf(xi, yi, zi, 30, cmap='jet')
    ax2.set_title(title2, y=1.02, fontsize=titlefontsize)
    g.set_clim(zmin, zmax)
    cbaxes = fig.add_axes([0.95, 0.1, 0.03, 0.8])
    plt.colorbar(g, cax=cbaxes)
    # plt.colorbar(g, ax=ax2)
    # plt.tight_layout()


def get_localQ(location, path, start=0, span=10):
    data = pd.read_table(location, sep='\s+', skiprows=1, names=["x", "y"] + ["Q" +str(i) for i in range(start, start+span)])
    d = data.dropna().values
    res = 30
    xi = np.linspace(min(d[:,0]), max(d[:,0]), res)
    yi = np.linspace(min(d[:,1]), max(d[:,1]), res)
    xv, yv = np.meshgrid(xi, yi)
    zi = griddata((d[:,0], d[:,1]), d[:,2:], (xv, yv), method='linear')
    nested_lst_of_tuples = [tuple(l) for l in path]
    tt = np.array([zi[l] for l in nested_lst_of_tuples])
    return tt



def shortest_path_2(location, temp="450", start=(4,5), end=-1, block=-1, res=30, zmin=0, zmax=30, xlabel="xlabel", ylabel="ylabel", title="AverageZ_Dis", save=False, plot1d=1, plot2d=True):
    data = np.loadtxt(location)
    xi, yi, zi = getxyz(data, res=res, zmin=zmin, zmax=zmax)
    zi = np.where(np.isnan(zi), 50, zi)

    if block != -1:
        x_low = np.searchsorted(xi, block[0])
        x_high = np.searchsorted(xi, block[1])
        y_low = np.searchsorted(yi, block[2])
        y_high = np.searchsorted(yi, block[3])
        zi[y_low:y_high, x_low:x_high] = 50
    V = ma.masked_array(zi, zi>40)
    G = nx.Graph()

    def func(u, v, d):
        node_u_wt = G.nodes[u].get('node_weight', 1)
        node_v_wt = G.nodes[v].get('node_weight', 1)
    #     edge_wt = d.get('weight', 1)
        edge_wt = 0
        return node_u_wt/2 + node_v_wt/2 + edge_wt
    n = len(xi)
    if end == -1:
        end = (n-5, n-5)
    # add nodes
    for i in range(n+1):
        for j in range(n+1):
            G.add_node((i,j))
    #         G.nodes[(i, j)]['node_weight'] = zi[i][j]
            # print(G)
            # continue
            # print("Hi")
            try:
                G.nodes[(i, j)]['node_weight'] = np.exp(zi[i][j])
    #             if zi[i][j] < 17:
    #                 G.nodes[(i, j)]['node_weight'] = zi[i][j]/10
    #             else:
    #                 G.nodes[(i, j)]['node_weight'] = zi[i][j]
            except IndexError:
                pass
    # add edges
    # connectivity = [(1,0), (0,1)]
    # connectivity = [(-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0)]
    connectivity = [(i,j) for i in [-1, 0, 1] for j in [-1, 0, 1] if (not (i == j == 0))]
    for i in range(n+1):
        for j in range(n+1):
            for (x,y) in connectivity:
                if i+x >= 0 and j+y >=0:
                    try:
                        G.add_edge((i,j), (i+x,j+y), weight=zi[i][j] + zi[i+x][j+y])
                        # G.add_edge((i,j), (i+x,j+y), weight=0)
                    except IndexError:
                        pass
                        # G.add_edge((i,j), (i+x,j+y), weight=max(zi))
    source = np.unravel_index(V.argmin(), V.shape)
    P = nx.dijkstra_path(G,start, end, weight=func)
    # P = nx.single_source_dijkstra(G,(4,5), weight=func)
    # P = nx.single_source_dijkstra(G,source, (4,5), weight=func)
    path = np.asarray(P)
    if plot2d:
        plt.contourf(xi, yi, V, res, cmap='jet')
        plt.plot(xi[path[:,1]], yi[path[:,0]], 'r.-')

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.clim(zmin, zmax)
        plt.colorbar()
        # plt.ylim(-8, 0)
        # plt.xlim(50, 70)
        # plt.clim(0, 20)
        # plt.colorbar()
    if save:
        plt.savefig(f"/Users/weilu/Dropbox/GlpG_paper_2018/figures/2d_{title}.png", dpi=300)
    f_on_path = [zi[tuple(p)] for p in reversed(path)]
    x_on_path = [xi[tuple(p)[1]] for p in reversed(path)]
    if plot1d == 2:
        plt.figure()
        x_on_path = np.array(x_on_path)
        d = pd.DataFrame(data={"x":x_on_path, "y":f_on_path})
        # mean the dupliation
        d = d.groupby("x").mean().reset_index().values
        x_smooth = np.linspace(d[:,0].min(), d[:,0].max(), 200)
        spl1 = scipy.interpolate.interp1d(d[:,0], d[:,1], kind="cubic")
        plt.plot(x_smooth, spl1(x_smooth))
        plt.xlabel("End to end distance(Å)")
        plt.ylabel("Free energy(kT)")
        # plt.plot(f_on_path)
        # plt.ylim([0,zmax])
        if save:
            plt.savefig(f"/Users/weilu/Dropbox/GlpG_paper_2018/figures/1d_path_{title}.png", dpi=300)
    if plot1d == 1:
        plt.figure()
        # x = np.array(range(len(f_on_path)))
        x = np.arange(len(f_on_path))
        x_smooth = np.linspace(x.min(), x.max(), 200)
        spl = scipy.interpolate.interp1d(x, f_on_path, kind="cubic")
        plt.plot(x_smooth, spl(x_smooth))
        # plt.plot(f_on_path)
        plt.ylim([0,zmax])
        if save:
            plt.savefig(f"/Users/weilu/Dropbox/GlpG_paper_2018/figures/1d_path_{title}.png", dpi=300)
    return (path, f_on_path)

def plot_shortest_path(location, path, res=30, zmin=0, zmax=30, z=3, xlabel="xlabel", ylabel="ylabel", title="", save=False, plot1d=True, plot2d=True, **kargs):
    data = np.loadtxt(location)
    xi, yi, zi = getxyz(data, res=res, zmin=zmin, zmax=zmax, z=z, **kargs)
    zi = np.where(np.isnan(zi), 50, zi)
    V = ma.masked_array(zi, zi>40)
    f_on_path = [zi[tuple(p)] for p in reversed(path)]
    if plot2d:
        plt.contourf(xi, yi, V, res, cmap='jet')
        plt.plot(xi[path[:,1]], yi[path[:,0]], 'r.-')
        plt.clim(zmin, zmax)
        # plt.colorbar()
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
    if save:
        plt.savefig(f"/Users/weilu/Dropbox/GlpG_paper_2018/figures/{title}.png", dpi=300)
    return (path, f_on_path)

def shortest_path(location, temp="450", start=(4,5), end=-1, block=-1, res=30, zmin=0, zmax=30, xlabel="xlabel", ylabel="ylabel", title="", save=False, plot1d=True, plot2d=True):
    data = np.loadtxt(location)
    xi, yi, zi = getxyz_3(data, res=res, zmin=zmin, zmax=zmax)
    zi = np.where(np.isnan(zi), 50, zi)
    V = ma.masked_array(zi, zi>40)
    if block != -1:
        x_low = np.searchsorted(xi, block[0])
        x_high = np.searchsorted(xi, block[1])
        y_low = np.searchsorted(yi, block[2])
        y_high = np.searchsorted(yi, block[3])
        zi[y_low:y_high, x_low:x_high] = 50

    G = nx.Graph()


    def func(u, v, d):
        node_u_wt = G.nodes[u].get('node_weight', 1)
        node_v_wt = G.nodes[v].get('node_weight', 1)
    #     edge_wt = d.get('weight', 1)
        edge_wt = 0
        return node_u_wt/2 + node_v_wt/2 + edge_wt
    n = len(xi)
    if end == -1:
        end = (n-5, n-5)
    # add nodes
    for i in range(n+1):
        for j in range(n+1):
            G.add_node((i,j))
    #         G.nodes[(i, j)]['node_weight'] = zi[i][j]
            try:
                G.nodes[(i, j)]['node_weight'] = np.exp(zi[i][j])
    #             if zi[i][j] < 17:
    #                 G.nodes[(i, j)]['node_weight'] = zi[i][j]/10
    #             else:
    #                 G.nodes[(i, j)]['node_weight'] = zi[i][j]
            except IndexError:
                pass
    # add edges
    # connectivity = [(1,0), (0,1)]
    # connectivity = [(-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0)]
    connectivity = [(i,j) for i in [-1, 0, 1] for j in [-1, 0, 1] if (not (i == j == 0))]
    for i in range(n+1):
        for j in range(n+1):
            for (x,y) in connectivity:
                if i+x >= 0 and j+y >=0:
                    try:
                        G.add_edge((i,j), (i+x,j+y), weight=zi[i][j] + zi[i+x][j+y])
                        # G.add_edge((i,j), (i+x,j+y), weight=0)
                    except IndexError:
                        pass
                        # G.add_edge((i,j), (i+x,j+y), weight=max(zi))
    source = np.unravel_index(V.argmin(), V.shape)
    P = nx.dijkstra_path(G,start, end, weight=func)
    # P = nx.single_source_dijkstra(G,(4,5), weight=func)
    # P = nx.single_source_dijkstra(G,source, (4,5), weight=func)
    path = np.asarray(P)
    if plot2d:
        plt.contourf(xi, yi, V, res, cmap='jet')
        plt.plot(xi[path[:,1]], yi[path[:,0]], 'r.-')
        plt.clim(zmin, zmax)
        plt.colorbar()
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
    if save:
        plt.savefig("/Users/weilu/Dropbox/GlpG_paper_2018/figures/2d_z6_qw.png", dpi=300)
    f_on_path = [zi[tuple(p)] for p in reversed(path)]
    if plot1d:
        plt.figure()
        plt.plot(f_on_path)
        plt.ylim([0,zmax])
        # if save:
        #     plt.savefig("/Users/weilu/papers/figures/shortest_path.png", dpi=300)
    return (path, f_on_path)

def plotPath(location, zmin=0, zmax=20, xlabel="xlabel", ylabel="ylabel", title="", outname=None):
    titlefontsize = 28
    data = np.loadtxt(location)
    xi, yi, zi = getxyz(data, zmin=zmin, zmax=zmax)
    # plt.contour(xi, yi, zi, 50, linewidths=0.25,colors='k')
    jet = cm = plt.get_cmap('jet')
    print(jet)
    # plt.contourf(xi, yi, zi, 20, cmap='rainbow')
    zi = np.where(np.isnan(zi), 50, zi)
    V = ma.masked_array(zi, zi>40)
    D, P = dijkstra(V)
    source = np.unravel_index(V.argmin(), V.shape)
    # source = (17,35)
    # print(P)
    # print(P)
    path = shortestPath(source, (3,3), P)
    # path = shortestPath((30,33), (3,3), P)
    plt.contourf(xi, yi, V, 30, cmap='jet')
    # plt.plot(path[:,1], path[:,0], 'r.-')
    plt.plot(xi[path[:,1]], yi[path[:,0]], 'r.-')
    # plt.xlim(xmin, xmax)
    plt.clim(zmin, zmax)
    plt.colorbar()

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title, y=1.02, fontsize=titlefontsize)
    # plt.tight_layout()
    # plt.axis('equal')
    # plt.axes().set_aspect('equal')
    # plt.axes().set_aspect('scaled')
    if outname:
        plt.savefig(outname, dpi=300, bbox_inches='tight')
    # plt.show()

def plotPath1d(location, zmin=0, zmax=20, xlabel="xlabel", ylabel="ylabel", title="", outname=None):
    titlefontsize = 28
    data = np.loadtxt(location)
    xi, yi, zi = getxyz(data, zmin=zmin, zmax=zmax)
    # plt.contour(xi, yi, zi, 50, linewidths=0.25,colors='k')
    jet = cm = plt.get_cmap('jet')
    print(jet)
    # plt.contourf(xi, yi, zi, 20, cmap='rainbow')
    zi = np.where(np.isnan(zi), 50, zi)
    V = ma.masked_array(zi, zi>40)
    D, P = dijkstra(V)
    source = np.unravel_index(V.argmin(), V.shape)
    # source = (17,35)
    # print(P)
    path = shortestPath(source, (3,3), P)
    f_on_path = [zi[tuple(p)] for p in path]
    # f_on_path = f_on_path[::-1]
    plt.figure()
    plt.plot(f_on_path)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title, y=1.02, fontsize=titlefontsize)
    # plt.tight_layout()
    # plt.axis('equal')
    # plt.axes().set_aspect('equal')
    # plt.axes().set_aspect('scaled')
    if outname:
        plt.savefig(outname, dpi=300, bbox_inches='tight')
    # plt.show()


def two_scales(ax1, time, data1, data2, c1, c2):
    """

    Parameters
    ----------
    ax : axis
        Axis to put two scales on

    time : array-like
        x-axis values for both datasets

    data1: array-like
        Data for left hand scale

    data2 : array-like
        Data for right hand scale

    c1 : color
        Color for line 1

    c2 : color
        Color for line 2

    Returns
    -------
    ax : axis
        Original axis
    ax2 : axis
        New twin axis
    """
    ax2 = ax1.twinx()

    ax1.plot(time, data1, color=c1)
    ax1.set_xlabel('')
    ax1.set_ylabel('FreeEnergy (kT)')

    ax2.plot(time, data2, color=c2)
    ax2.set_ylabel('Expected Value')
    # ax2.set_ylabel('Expected End-to-End distance (Å)')
    # print("hi")
    return ax1, ax2

# Change color of each axis
def color_y_axis(ax, color):
    """Color your axes."""
    for t in ax.get_yticklabels():
        t.set_color(color)
    return None

# # https://bougui505.github.io/2016/08/31/compute_the_shortest_path_on_a_grid_using_python.html
# def dijkstra(V):
#     mask = V.mask
#     visit_mask = mask.copy() # mask visited cells
#     m = np.ones_like(V) * np.inf
#     connectivity = [(i,j) for i in [-1, 0, 1] for j in [-1, 0, 1] if (not (i == j == 0))]
#     cc = np.unravel_index(V.argmin(), m.shape) # current_cell
#     m[cc] = 0
#     P = {}  # dictionary of predecessors
#     # while (~visit_mask).sum() > 0:
#     for _ in range(V.size):
#         # print cc
#         neighbors = [tuple(e) for e in np.asarray(cc) - connectivity
#                      if e[0] > 0 and e[1] > 0 and e[0] < V.shape[0] and e[1] < V.shape[1]]
#         neighbors = [ e for e in neighbors if not visit_mask[e] ]
#         tentative_distance = np.asarray([V[e]-V[cc] for e in neighbors])
#         for i,e in enumerate(neighbors):
#             d = tentative_distance[i] + m[cc]
#             if d < m[e]:
#                 m[e] = d
#                 P[e] = cc
#         visit_mask[cc] = True
#         m_mask = ma.masked_array(m, visit_mask)
#         cc = np.unravel_index(m_mask.argmin(), m.shape)
#     return m, P

# def shortestPath(start, end, P):
#     Path = []
#     step = end
#     while 1:
#         Path.append(step)
#         if step == start:
#             break
#         step = P[step]
#     Path.reverse()
#     return np.asarray(Path)


# def plot2d(location, temp="450", zmin=0, zmax=30, xlabel="xlabel", ylabel="ylabel", title="", outname=None):
#     titlefontsize = 28
#     filename = location + f"pmf-{temp}.dat"
#     data = np.loadtxt(filename)
#     xi, yi, zi = getxyz(data, zmin=zmin, zmax=zmax)
#     # plt.contour(xi, yi, zi, 50, linewidths=0.25,colors='k')
#     jet = cm = plt.get_cmap('jet')
#     print(jet)
#     # plt.contourf(xi, yi, zi, 20, cmap='rainbow')
#     plt.figure()
#     plt.contourf(xi, yi, zi, 30, cmap='jet')
#     # plt.xlim(xmin, xmax)
#     plt.clim(zmin, zmax)
#     plt.colorbar()

#     plt.xlabel(xlabel)
#     plt.ylabel(ylabel)
#     plt.title(title, y=1.02, fontsize=titlefontsize)
#     # plt.tight_layout()
#     # plt.axis('equal')
#     # plt.axes().set_aspect('equal')
#     # plt.axes().set_aspect('scaled')
#     if outname:
#         plt.savefig(outname, dpi=300, bbox_inches='tight')
#     # plt.show()