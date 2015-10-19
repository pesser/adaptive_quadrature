#include <multitria/find_cells.h>
#include <multitria/find_values.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_q.h>
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

double diagonal2(double x)
{
  return x;
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


void test_multi_interface_quadrature(std::function<double(double)> interface1,
                                     std::function<double(double)> interface2)
{
  const int dim = 2;
  const unsigned int n_adaptive_cycles = 4;

  Triangulation<dim> tria1;
  Triangulation<dim> tria2;
  GridGenerator::hyper_cube(tria1);
  GridGenerator::hyper_cube(tria2);
  tria1.refine_global(1);
  tria2.refine_global(1);
  for(unsigned int s = 0; s < n_adaptive_cycles; ++s)
  {
    adapt_tria(tria1, interface1);
    adapt_tria(tria2, interface2);
    output_grid<dim>(tria1, "interface1_test", s);
    output_grid<dim>(tria2, "interface2_test", s);
  }

  std::vector<Point<dim> > qpoints(0);
  std::vector<double> qweights(0);
  std::vector<FEFunctionData<dim, Vector<double> >*> fef_datas;

  QGauss<dim> quadrature(2);

  // fe function dummies
  FE_Q<dim> fe1(1);
  DoFHandler<dim> dofh1(tria1);
  dofh1.distribute_dofs(fe1);
  FEValues<dim> fev1(fe1, quadrature, update_values);
  std::vector<double> fef1_values(0);
  Vector<double> input1(dofh1.n_dofs());
  FEFunctionData<dim, Vector<double> > fef1;
  fef1.values = &fef1_values;
  fef1.global_dofs = &input1;
  fef1.fev = &fev1;
  fef1.cell = dofh1.begin();

  FE_Q<dim> fe2(1);
  DoFHandler<dim> dofh2(tria2);
  dofh2.distribute_dofs(fe2);
  FEValues<dim> fev2(fe2, quadrature, update_values);
  std::vector<double> fef2_values(0);
  Vector<double> input2(dofh2.n_dofs());
  FEFunctionData<dim, Vector<double> > fef2;
  fef2.values = &fef2_values;
  fef2.global_dofs = &input2;
  fef2.fev = &fev2;
  fef2.cell = dofh2.begin();

  fef_datas.push_back(&fef1);
  fef_datas.push_back(&fef2);

  collect_quadratures_and_values(qpoints, qweights, fef_datas);

  unsigned int n_qpoints = qpoints.size();
  assert(qweights.size() == n_qpoints);
  deallog << "Collected quadrature size: " << n_qpoints << std::endl;
  // quadrature weights should sum up to 1
  double sum = 0;
  for(unsigned int q = 0; q < n_qpoints; ++q)
  {
    sum += qweights[q];
    const Point<dim>& q_point = qpoints[q];
    for(unsigned int d = 0; d < dim; ++d)
    {
      deallog << q_point(d) << " ";
    }
    deallog << qweights[q];
    deallog << std::endl;
  }
  deallog << "Collected quadrature weight: " << sum << std::endl;

  assert(fef1.values->size() == n_qpoints);
  assert(fef2.values->size() == n_qpoints);
}

void test4()
{
  test_multi_interface_quadrature(diagonal, diagonal2);
}

void test5()
{
  test_multi_interface_quadrature(diagonal, parabola);
}

int main()
{
  std::ofstream logfile("deallog");
  deallog.attach(logfile);

  test5();

  return 0;
}
