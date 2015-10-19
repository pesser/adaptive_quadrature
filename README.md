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
![Parabola grid](parabola_grid.svg)
The generated quadrature rule
![Parabola quadrature](parabola_quadrature.png)

#### Diagonal and Parabola
Combined quadrature rule for Diagonal and Parabola interface (`test5()`)
![Diagonal and parabola quadrature](diagonal_and_parabola_quadrature.png)

#### Cube and circle
Grid for cube interface
![Cube grid](cube.svg)
Grid for circle interface
![Circle grid](circle.svg)
Combined quadrature rule for Cube and Circle interface (`test_values.cc`)
![Cube and Circle quadrature](cube_and_circle.png)

### Running
Make sure the environment variable `DEAL_II_DIR` contains the location of
`deal.II`, then

    mkdir build
    cd build
    cmake ..
    make qtest && ./multitria/qtest
    ./plot.sh deallog
