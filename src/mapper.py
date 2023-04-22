import numpy as np
import matplotlib.pyplot as plt

from src.cover_complex import *

from gudhi import SimplexTree
from gudhi import CoverComplex

import networkx as nx

# test the mapper on simple example :

## compute the mapper graph associated to this point cloud : 

t = np.linspace(0,10,300)
x_circle = np.sin(t)
y_circle = np.cos(t) + np.random.randn(len(t))/10

x_circle_2 = np.sin(t)
y_circle_2 = np.cos(t) -3 + np.random.randn(len(t))/10

x_center = np.random.randn(400)/10
y_center = np.random.randn(400)/3 - 1.5

x = np.concatenate((x_circle, x_circle_2, x_center))
y = np.concatenate((y_circle, y_circle_2, y_center))
plt.scatter(x,y)
plt.show()

X = np.concatenate((x.reshape(-1,1),y.reshape(-1,1)), axis=1)
X.shape

# filter = projection on y axis
X_filter = X[:,1:]

mapper = MapperComplex(input_type="point cloud", )
mapper.fit(X, filters=X_filter)
G = mapper.get_networkx()

# we obtain the expected graph : 2 linked circles !
plt.figure()
nx.draw_networkx(G)
plt.show()

## Compute the mapper on hi-c dataset

# Compute pairwise SCC matrix