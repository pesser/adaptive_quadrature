#ifndef __multitria_find_values_h
#define __multitria_find_values_h

#include "find_cells.h"
#include <deal.II/fe/fe_values.h>

/**
 * Only for one scalar valued fe functions
 */
  template <int dim, class InputVector>
void collect_quadratures_and_values(
    std::vector<dealii::Point<dim> >& qpoints,
    std::vector<double>& qweights,
    std::vector<double>& values,
    dealii::FEValues<dim>& fev,
    const InputVector& input,
    typename dealii::DoFHandler<dim>::cell_iterator cell
    )
{
  if(cell->active())
  {
    // not refined, return copy of base quadrature
    qpoints.insert(qpoints.end(),
                   fev.get_quadrature().get_points().begin(),
                   fev.get_quadrature().get_points().end());
    qweights.insert(qweights.end(),
                   fev.get_quadrature().get_weights().begin(),
                   fev.get_quadrature().get_weights().end());
    fev.reinit(cell);
    fev.get_function_values(input, values);
    return;
  }
  for(unsigned int child = 0;
      child < dealii::GeometryInfo<dim>::max_children_per_cell;
      ++child)
  {
    // get child
    typename dealii::DoFHandler<dim>::cell_iterator child_cell =
      cell->child(child);
    // collect child's stuff
    collect_quadratures_and_values(qpoints, qweights,
                                   values,
                                   fev,
                                   input,
                                   child_cell);
  }
}

/*
 * Just a simple datastructure to collect different parts we need to
 * evaluate fe functions.
 */
template <int dim, class InputVector>
struct FEFunctionData
{
  std::vector<double>* values; // gets appended

  dealii::FEValues<dim>* fev;
  InputVector* global_dofs;
  typename dealii::DoFHandler<dim>::cell_iterator cell;
};

// Notes about this process:
// - We do not know a priori which quadrature rule we need to use,
//   thus it must be fast to calculate stuff for dynamic quadratures
// - Common to all calculations: Need to know local dof indices
//
// Two approaches:
// - Like we do here: Walk down cells, if a fe function is not further
//   refined, wait until quadrature rule on this cell is calculated from
//   refined fe functions. Afterwards, use this quadrature rule to evaluate
//   the non-refined fe function.
// Notice however that evaluation requires inverse Jacobians, which were
// actually already computed on the fine cells. Thus another approach would
// consists in decoupling fe evaluation on the unit cell (which would only
// require the global information about the dof indices and the quadrature
// points on the reference cell) and mapping evaluation on the physical cell
// or even all three (as in MatrixFree?)
// - coefficient index calculation
// - reference cell evaluation of basis
// - physical cell evaluation of mapping

/**
 * Only for scalar valued fe functions
 */
  template <int dim, class InputVector>
void collect_quadratures_and_values(
    std::vector<dealii::Point<dim> >& qpoints, // append
    std::vector<double>& qweights, // append
    std::vector<FEFunctionData<dim, InputVector>*> fef_datas) //append
{
  // split into refined and non refined fe functions
  std::vector<FEFunctionData<dim, InputVector>* > active_fef(0);
  std::vector<FEFunctionData<dim, InputVector>* > refined_fef(0);
  std::vector<typename dealii::DoFHandler<dim>::cell_iterator> refined_cells;
  for(auto fef : fef_datas)
  {
    if(fef->cell->active())
    {
      active_fef.push_back(fef);
    } else {
      refined_fef.push_back(fef);
      refined_cells.push_back(fef->cell);
    }
  }

  if(refined_fef.size() == 0)
  {
    // current cell not refined, use given quadrature and evaluate all functions
    Assert(active_fef.size() > 0, dealii::ExcInternalError());
    auto first_fef = *(active_fef.begin());
    qpoints.insert(qpoints.end(),
                   first_fef->fev->get_quadrature().get_points().begin(),
                   first_fef->fev->get_quadrature().get_points().end());
    qweights.insert(qweights.end(),
                    first_fef->fev->get_quadrature().get_weights().begin(),
                    first_fef->fev->get_quadrature().get_weights().end());
    for(auto fef : active_fef)
    {
      // TODO fevalues is horrible here - we would at the _very_ least need
      // a get_values that takes iterator ranges
      fef->fev->reinit(fef->cell);
      std::vector<double> tmp_values(fef->fev->n_quadrature_points);
      fef->fev->get_function_values(*(fef->global_dofs), tmp_values);
      fef->values->insert(fef->values->end(),
                          tmp_values.begin(),
                          tmp_values.end());
    }
    return;
  } else {
    // collect sub-quadrature points and values of refined functions
    std::vector<dealii::Point<dim> > sub_qpoints(0);
    std::vector<double> sub_qweights(0);
    // for each child, set correct cell in active_fef,
    // collect quadrature and values on child
    // merge quadratures and values
    for(unsigned int child = 0;
        child < dealii::GeometryInfo<dim>::max_children_per_cell;
        ++child)
    {
      // collect child quadrature
      std::vector<dealii::Point<dim> > sub_child_qpoints(0);
      std::vector<double> sub_child_qweights(0);
      for(unsigned int i = 0; i < refined_fef.size(); ++i)
      {
        refined_fef[i]->cell = refined_cells[i]->child(child);
      }
      collect_quadratures_and_values(
          sub_child_qpoints, sub_child_qweights,
          refined_fef);
      // now sub_child_qpoints contains quadrature on child and
      // refined fe functions had their values appended
      for(unsigned int q = 0; q < sub_child_qpoints.size(); ++q)
      {
        // append projected child quadrature to sub quadrature
        sub_qpoints.push_back(dealii::GeometryInfo<dim>::child_to_cell_coordinates(
                sub_child_qpoints[q], child));
        // append rescaled weights to sub weights
        sub_qweights.push_back(
            sub_child_qweights[q] /
            double(dealii::GeometryInfo<dim>::max_children_per_cell));
      }
    }
    // now sub_qpoints contains the quadrature rule for the current cell
    // use this to evaluate (and append) the active fe functions
    for(auto fef : active_fef)
    {
      // cannot use the given fevalues object as we need to use the
      // quadrature rule just constructed
      dealii::Quadrature<dim> sub_quadrature(sub_qpoints, sub_qweights);
      dealii::FEValues<dim> q_adjusted_fev (fef->fev->get_mapping(),
                                            fef->fev->get_fe(),
                                            sub_quadrature,
                                            fef->fev->get_update_flags());
      q_adjusted_fev.reinit(fef->cell);
      std::vector<double> tmp_values(sub_quadrature.size());
      q_adjusted_fev.get_function_values(*(fef->global_dofs), tmp_values);
      fef->values->insert(fef->values->end(),
                          tmp_values.begin(),
                          tmp_values.end());
    }
    // finally append sub_qpoints to qpoints
    qpoints.insert(qpoints.end(),
                   sub_qpoints.begin(),
                   sub_qpoints.end());
    qweights.insert(qweights.end(),
                    sub_qweights.begin(),
                    sub_qweights.end());
  }
}

#endif
