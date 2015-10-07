#include <multitria/find_cells.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/base/logstream.h>
#include <fstream>

using namespace dealii;

template <int dim>
void output_grid(const Triangulation<dim>& tria,
                 std::string name,
                 const unsigned int nr)
{
  GridOut grid_out;
  std::stringstream filename;
  filename << name << "-" << nr << ".svg";
  std::ofstream out(filename.str());
  grid_out.write_svg(tria, out);
}

void test1()
{
  const int dim = 2;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);

  QGauss<dim> quadrature(2);
  Triangulation<dim>::cell_iterator root = tria.begin();
  Quadrature<dim> collected_quadrature = collect_quadratures(tria.begin(),
                                                             &quadrature);

  unsigned int n_qpoints = collected_quadrature.size();
  deallog << "Collected quadrature size: " << n_qpoints << std::endl;
  // quadrature weights should sum up to 1
  double sum = 0;
  for(unsigned int q = 0; q < collected_quadrature.size(); ++q)
  {
    sum += collected_quadrature.weight(q);
    const Point<dim>& q_point = collected_quadrature.point(q);
    for(unsigned int d = 0; d < dim; ++d)
    {
      deallog << q_point(d) << " ";
    }
    deallog << std::endl;
  }
  deallog << "Collected quadrature weight: " << sum << std::endl;
}


/* determine whether a given cell contains the interface described by the
 * graph of f
 */
  template <int dim>
bool contains_interface(typename Triangulation<dim>::cell_iterator cell,
                        std::function<double(double)> f,
                        unsigned int n_samples)
{
  // only 2d - see GeometryInfo for vertex numbering
  const Point<dim>& lower_left = cell->vertex(0);
  const Point<dim>& lower_right = cell->vertex(1);

  double min_x = lower_left(0);
  double max_x = lower_right(0);
  double delta_x = (max_x - min_x)/(double(n_samples) + 1.0);
  for(unsigned int s = 1; s <= n_samples; ++s)
  {
    double x = min_x + (double)s * delta_x;
    double y = f(x);
    Point<dim> p(x, y);
    if(cell->point_inside(p))
    {
      return true;
    }
  }
  return false;
}

  template <int dim>
void adapt_tria(Triangulation<dim>& tria,
                std::function<double(double)> interface)
{
  const unsigned int n_samples = 10;
  for(typename Triangulation<dim>::active_cell_iterator cell = tria.begin_active();
      cell != tria.end();
      ++cell)
  {
    if(contains_interface<dim>(cell, interface, n_samples))
    {
      cell->set_refine_flag();
    }
  }
  tria.execute_coarsening_and_refinement();
}

void test_interface_quadrature(std::function<double(double)> interface)
{
  const int dim = 2;
  const unsigned int n_adaptive_cycles = 4;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);
  for(unsigned int s = 0; s < n_adaptive_cycles; ++s)
  {
    adapt_tria(tria, interface);
    output_grid<dim>(tria, "test", s);
  }

  QGauss<dim> quadrature(2);
  Triangulation<dim>::cell_iterator root = tria.begin();
  Quadrature<dim> collected_quadrature = collect_quadratures(tria.begin(),
                                                             &quadrature);

  unsigned int n_qpoints = collected_quadrature.size();
  deallog << "Collected quadrature size: " << n_qpoints << std::endl;
  // quadrature weights should sum up to 1
  double sum = 0;
  for(unsigned int q = 0; q < collected_quadrature.size(); ++q)
  {
    sum += collected_quadrature.weight(q);
    const Point<dim>& q_point = collected_quadrature.point(q);
    for(unsigned int d = 0; d < dim; ++d)
    {
      deallog << q_point(d) << " ";
    }
    deallog << std::endl;
  }
  deallog << "Collected quadrature weight: " << sum << std::endl;
}

double diagonal(double x)
{
  return 1.0 - x;
}

void test2()
{
  test_interface_quadrature(diagonal);
}

double parabola(double x)
{
  return 0.7 - x*x;
}

void test3()
{
  test_interface_quadrature(parabola);
}

int main()
{
  std::ofstream logfile("deallog");
  deallog.attach(logfile);

  test3();

  return 0;
}
