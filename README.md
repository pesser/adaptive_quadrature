### Adaptive Quadrature Rules
Create a quadrature rule on a given cell, based on the refinement of this
cell in another triangulation.

### Examples
#### Diagonal
The grid to resolve the diagonal interface
![Diagonal grid](diagonal_grid.svg)
The generated quadrature rule
![Diagonal quadrature](diagonal_quadrature.png)

#### Parabola
The grid to resolve the parabola interface
![Diagonal grid](parabola_grid.svg)
The generated quadrature rule
![Diagonal quadrature](parabola_quadrature.png)

### Running
Make sure the environment variable `DEAL_II_DIR` contains the location of
`deal.II`, then

    mkdir build
    cd build
    cmake ..
    make qtest && ./multitria/qtest
    ./plot.sh deallog
