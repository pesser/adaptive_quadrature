import matplotlib.pyplot as plt
import numpy as np


def read_data(filename):
    data =[]
    with open(filename) as f:
        all_data = f.read()
        for line in all_data.splitlines():
            columns = line.split()
            data.append(
                    (float(columns[0]), float(columns[1])))
    return data

fig = plt.figure()
ax = fig.add_subplot(111)
x_data_label = "x"
y_data_label = "y"
title = "quadrature points"

ax.set_xlabel(x_data_label)
ax.set_ylabel(y_data_label)
ax.set_aspect("equal")

data = read_data("data")
for x, y in data:
    print("{} {}".format(x, y))
    ax.plot(x, y, 'o', color = "black")


plt.title(title)
plt.show()
