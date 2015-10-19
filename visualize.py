import matplotlib.pyplot as plt
import numpy as np


def read_data(filename):
    data =[]
    weights = []
    with open(filename) as f:
        all_data = f.read()
        for line in all_data.splitlines():
            columns = line.split()
            data.append(
                    (float(columns[0]), float(columns[1])))
            if len(columns) > 2:
                weights.append(float(columns[2]))
    return data, weights

fig = plt.figure()
ax = fig.add_subplot(111)
x_data_label = "x"
y_data_label = "y"
title = "quadrature points"

ax.set_xlabel(x_data_label)
ax.set_ylabel(y_data_label)
ax.set_aspect("equal")

data, weights = read_data("data")
if len(weights) == len(data):
    print("Using weights")
else:
    print(len(weights))
    print(len(data))
for i, (x, y) in enumerate(data):
    if len(weights) == len(data):
        ax.plot(x, y, 'o', color = "black", markersize = (weights[i] * 2 * len(weights)))
    else:
        ax.plot(x, y, 'o', color = "black")


plt.title(title)
plt.show()
