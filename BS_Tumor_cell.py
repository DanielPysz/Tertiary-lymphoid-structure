import pandas as pd
import numpy as np
from sklearn.neighbors import radius_neighbors_graph
import matplotlib.pyplot as plt
from sklearn.neighbors import KDTree
import networkx as nx
import matplotlib.pyplot as plt

df = pd.read_csv("C:/Users/danpy/OneDrive/Pulpit/if_data/0802_IF1.csv", index_col='cell.ID')
dic = pd.read_csv("C:/Users/danpy/OneDrive/Pulpit/IF1_phen_to_cell_mapping.csv")
dict = {}
for i,j in dic.iterrows():
    dict[j['phenotype']] = j['celltype']
x = list(df['nucleus.x'])
y = list(df['nucleus.y'])
CD15 = list(df["CD15.score.normalized"])
CK = list(df["CK.score.normalized"])
CD3 = list(df["CD3.score.normalized"])
CD11c = list(df["CD11c.score.normalized"])
CD20 = list(df["CD20.score.normalized"])
CD163 = list(df["CD163.score.normalized"])
phenotype = []
for i in np.arange(len(df)):
    ph = []
    if CD11c[i] > 1:
        ph.append("CD11c+")
    else:
        ph.append("CD11C-")
    if CD15[i] > 1:
        ph.append("CD15+")
    else:
        ph.append("CD15-")
    if CD163[i] > 1:
        ph.append("CD163+") 
    else:
        ph.append("CD163-")
    if CD20[i] > 1:
        ph.append("CD20+")
    else:
        ph.append("CD20-")
    if CD3[i] > 1:
        ph.append("CD3+")
    else:
        ph.append("CD3-")
    if CK[i] > 1:
        ph.append("CK+")
    else:
        ph.append("CK-")
    ph = "".join(ph)
    phenotype.append(ph)
df["phenotype"] = phenotype

type = []
for i in list(df['phenotype']):
    if i in list(dict.keys()):
        for j in list(dict.keys()):
            if i == j:
                type.append(dict[j]) 
    else:
        type.append('other')
df['celltype'] = type

ind = []
for i in list(df['celltype']):
    if str(i) == 'Bcell':
        ind.append(1)
    elif str(i) =='Tcell':
        ind.append(2)
    elif str(i) =='Macrophage':
        ind.append(3)
    else:
        ind.append(0)
df['index'] = ind  

def coord(x, y, cellt):
    global df
    coordinate = []
    for i, j in enumerate(list(df['celltype'])):
        if j == cellt:
            coordinate.append([x[i], y[i]])
    return coordinate

def scatter(x, y, cellt):
    global df
    scatx = []
    scaty = []
    for i,j in enumerate(list(df['celltype'])):
        if j == cellt:
            scatx.append(x[i])
            scaty.append(y[i])
    return scatx, scaty

coordinate = coord(x, y, "Bcell")
r = 30
graf = radius_neighbors_graph(coordinate, radius=r, mode='distance', include_self=False)

def visualize_graph(graph):
    nx_graph = nx.from_scipy_sparse_array(graph)
    pos = nx.spring_layout(nx_graph)
    nx.draw_networkx_nodes(nx_graph, pos, node_size=20, node_color='lightblue')
    nx.draw_networkx_edges(nx_graph, pos, edge_color='gray')
    nx.draw_networkx_labels(nx_graph, pos, font_color='black', font_size = 8)
    plt.axis('off')
    plt.show()

def extract_strongly_connected_components(graph):
    nx_graph = nx.from_scipy_sparse_array(graph)
    if nx.is_directed(nx_graph):
        components = list(nx.strongly_connected_components(nx_graph))
    else:
        components = list(nx.connected_components(nx_graph))
    return components

def component_plot(grafB, grafT):
    labx = []
    laby = []
    nx_graph = nx.from_scipy_sparse_array(grafB)
    pos = nx.spring_layout(nx_graph)
    components = extract_strongly_connected_components(grafB)
    for component in components:
        if len(component) > 1:
            for i in component:
                cellsk = pos.keys()
                for n, k in enumerate(cellsk):
                     if n == i:
                        labx.append(pos[k][0])
                        laby.append(pos[k][1])
    plt.scatter(labx,laby, color='red', label = 'BCell')
    labx = []
    laby = []
    nx_graph = nx.from_scipy_sparse_array(grafB)
    pos = nx.spring_layout(nx_graph)
    components = extract_strongly_connected_components(grafB)
    for component in components:
        if len(component) > 1:
            for i in component:
                cellsk = pos.keys()
                for n, k in enumerate(cellsk):
                     if n == i:
                        labx.append(pos[k][0])
                        laby.append(pos[k][1])
    plt.scatter(labx, laby, color = 'orange', label='Tcell')
    plt.legend()
    plt.show()

def scatter():
    Bx, By = scatter(x, y, "Bcell")
    Tx, Ty = scatter(x, y, "Tcell")
    Mx, My = scatter(x, y, 'Macrophage')
    plt.scatter(x, y, color='grey', label = 'other')
    plt.scatter(Mx, My, color = 'orange', label = 'Macrophage')
    plt.scatter(Bx, By, color='blue', label = 'B Cell')
    plt.scatter(Tx, Ty, color='red', label = 'T Cell')
    plt.legend(loc = 'upper left')
    plt.show()

grafB = radius_neighbors_graph(coord(x, y, "Bcell"), radius=r, mode='distance', include_self=True)
grafT = radius_neighbors_graph(coord(x, y, "Tcell"), radius=r, mode='distance', include_self=True)
component_plot(grafB, grafT)


