#include <multitria/find_cells.h>
#include <multitria/find_values.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/base/logstream.h>
#include <deal.II/numerics/vector_tools.h>
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


/* determine whether a given cell contains the interface described by the
 * boundary of the set represented by the given characteristic function
 */
  template <int dim>
bool contains_interface(typename Triangulation<dim>::cell_iterator cell,
                        const Function<dim, double>& cf,
                        unsigned int n_samples)
{
  unsigned int n_inside = 0;
  unsigned int n_outside = 0;
  QGauss<dim> quadrature(n_samples);

  for(unsigned int q = 0; q < quadrature.size(); ++q)
  {
    Point<dim> qpoint = StaticMappingQ1<dim>::mapping.transform_unit_to_real_cell(
        cell, quadrature.point(q));
    double value = cf.value(qpoint);
    if(value > 0.5)
    {
      ++n_inside;
    } else {
      ++n_outside;
    }
    if(n_inside > 0 && n_outside > 0)
    {
      return true;
    }
  }
  return false;
}

/* adapt tria based on characteristic function.
 */
  template <int dim>
void adapt_tria(Triangulation<dim>& tria,
                const Function<dim, double>& cf)
{
  const unsigned int n_samples = 60;
  for(typename Triangulation<dim>::active_cell_iterator cell = tria.begin_active();
      cell != tria.end();
      ++cell)
  {
    if(contains_interface<dim>(cell, cf, n_samples))
    {
      cell->set_refine_flag();
    }
  }
  tria.execute_coarsening_and_refinement();
}


template <int dim>
class Ball : public Function<dim>
{
  public:
    Ball(Point<dim> center, double r) : center(center), r(r)
  {}

    virtual double value(const dealii::Point<dim>& p,
                         const unsigned int) const
    {
      Tensor<1, dim> center_to_p = p - this->center;
      if(center_to_p.norm() < r)
      {
        return 1.0;
      }
      return 0.0;
    }

  protected:
    const Point<dim> center;
    const double r;
};


template <int dim>
class Cube : public Function<dim>
{
  public:
    Cube(Point<dim> center, double l) : center(center), l(l)
  {}

    virtual double value(const dealii::Point<dim>& p,
                         const unsigned int) const
    {
      Tensor<1, dim> center_to_p = p - this->center;
      for(unsigned int d = 0; d < dim; ++d)
      {
        if(std::abs(center_to_p[d]) >= 0.5 * this->l)
        {
          return 0.0;
        }
      }
      return 1.0;
    }

  protected:
    const Point<dim> center;
    const double l;
};


double test_quadrature_values(const Function<2, double>& cf1,
                              const Function<2, double>& cf2,
                              const unsigned int n_adaptive_cycles = 6)
{
  const int dim = 2;

  Triangulation<dim> tria1;
  Triangulation<dim> tria2;
  GridGenerator::hyper_cube(tria1);
  GridGenerator::hyper_cube(tria2);
  for(unsigned int s = 0; s < n_adaptive_cycles; ++s)
  {
    adapt_tria(tria1, cf1);
    adapt_tria(tria2, cf2);
    output_grid<dim>(tria1, "characteristic1_test", s);
    output_grid<dim>(tria2, "characteristic2_test", s);
  }

  std::vector<Point<dim> > qpoints(0);
  std::vector<double> qweights(0);
  std::vector<FEFunctionData<dim, Vector<double> >*> fef_datas;

  QGauss<dim> quadrature(2);

  // fe function dummies
  FE_DGQ<dim> fe1(1);
  DoFHandler<dim> dofh1(tria1);
  dofh1.distribute_dofs(fe1);
  FEValues<dim> fev1(fe1, quadrature, update_values);
  std::vector<double> fef1_values(0);
  Vector<double> input1(dofh1.n_dofs());
  VectorTools::interpolate(dofh1, cf1, input1);
  FEFunctionData<dim, Vector<double> > fef1;
  fef1.values = &fef1_values;
  fef1.global_dofs = &input1;
  fef1.fev = &fev1;
  fef1.cell = dofh1.begin();

  FE_DGQ<dim> fe2(1);
  DoFHandler<dim> dofh2(tria2);
  dofh2.distribute_dofs(fe2);
  FEValues<dim> fev2(fe2, quadrature, update_values);
  std::vector<double> fef2_values(0);
  Vector<double> input2(dofh2.n_dofs());
  VectorTools::interpolate(dofh2, cf2, input2);
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
  //deallog << "Collected quadrature size: " << n_qpoints << std::endl;
  // quadrature weights should sum up to 1
  double sum = 0;
  for(unsigned int q = 0; q < n_qpoints; ++q)
  {
    sum += qweights[q];
    const Point<dim>& q_point = qpoints[q];
    for(unsigned int d = 0; d < dim; ++d)
    {
      //deallog << q_point(d) << " ";
    }
    //deallog << qweights[q];
    //deallog << std::endl;
  }
  //deallog << "Collected quadrature weight: " << sum << std::endl;

  assert(fef1.values->size() == n_qpoints);
  assert(fef2.values->size() == n_qpoints);

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  FE_DGQ<dim> fe(0);
  DoFHandler<dim> dofh(tria);
  dofh.distribute_dofs(fe);
  Quadrature<dim> collected_q(qpoints, qweights);
  FEValues<dim> fev(fe, collected_q, update_values | update_JxW_values);
  Vector<double> input(dofh.n_dofs());
  assert(input.size() == 1);
  input[0] = 1.0; //constant one function for illustration
  fev.reinit(dofh.begin_active()); // reinit on active cell - just have one
  std::vector<double> values(collected_q.size());
  fev.get_function_values(input, values);

  double integral = 0.0;
  const std::vector<double> values_fe1 = *(fef_datas[0]->values);
  const std::vector<double> values_fe2 = *(fef_datas[1]->values);
  for(unsigned int q = 0; q < collected_q.size(); ++q)
  {
    integral += values[q] * values_fe1[q] * values_fe2[q] * fev.JxW(q);
  }
  deallog << collected_q.size() << " " << integral << std::endl;
  return integral;
}

double test_quadrature_values_global(const Function<2, double>& cf1,
                                     const Function<2, double>& cf2,
                                     const unsigned int n_global_cycles = 6)
{
  const int dim = 2;

  unsigned int n_qpoints = 2;
  for(unsigned int s = 0; s < n_global_cycles; ++s)
  {
    n_qpoints *= 2;
  }

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  FE_DGQ<dim> fe(0);
  DoFHandler<dim> dofh(tria);
  dofh.distribute_dofs(fe);
  QGauss<dim> quadrature(n_qpoints);
  FEValues<dim> fev(fe, quadrature,
                    update_values | update_JxW_values | update_quadrature_points);
  Vector<double> input(dofh.n_dofs());
  assert(input.size() == 1);
  input[0] = 1.0; //constant one function for illustration
  fev.reinit(dofh.begin_active()); // reinit on active cell - just have one
  std::vector<double> values(quadrature.size());
  fev.get_function_values(input, values);

  double integral = 0.0;
  std::vector<double> values_fe1(quadrature.size());
  std::vector<double> values_fe2(quadrature.size());
  cf1.value_list(fev.get_quadrature_points(), values_fe1);
  cf2.value_list(fev.get_quadrature_points(), values_fe2);
  for(unsigned int q = 0; q < quadrature.size(); ++q)
  {
    integral += values[q] * values_fe1[q] * values_fe2[q] * fev.JxW(q);
  }
  deallog << quadrature.size() << " " << integral << std::endl;
  return integral;
}

void test1()
{
  const unsigned int steps = 10;
  // perturbation to avoid grid alignment
  double delta = 0.01;
  // geometry
  double c = 0.6 + delta;
  double r = 0.2;
  Ball<2> ball(Point<2>(c, c), r);
  Cube<2> cube(Point<2>(c-r, c-r), 2*r);

  const double pi = 3.14159265359;
  const double expected_integral = pi * r * r / 4.0;
  deallog << "Expected: " << expected_integral << std::endl;
  deallog << "Adaptive:" << std::endl;
  for(unsigned int s = 0; s < steps; ++s)
  {
    double integral = test_quadrature_values(ball, cube, s);
    //deallog << "Error " << s << " " << std::abs(integral - expected_integral) << std::endl;
  }

  deallog << "Global:" << std::endl;
  for(unsigned int s = 0; s < steps; ++s)
  {
    double integral = test_quadrature_values_global(ball, cube, s);
    //deallog << "Error " << s << " " << std::abs(integral - expected_integral) << std::endl;
  }
}

int main()
{
  std::ofstream logfile("deallog");
  deallog.attach(logfile);

  test1();

  return 0;
}
