

#pragma once

#include "INSFVFluxKernel.h"
#include "INSFVMomentumResidualObject.h"

// Forward declare variable class
class INSFVVelocityVariable;

class INSFVEddyViscCSVReynoldsStress : public INSFVFluxKernel
{
public:
  static InputParameters validParams();

  INSFVEddyViscCSVReynoldsStress(const InputParameters & params);

  using INSFVFluxKernel::gatherRCData;
  void gatherRCData(const FaceInfo &) override final;

protected:
  /**
   * Routine to compute this object's strong residual (e.g. not multipled by area). This routine
   * should also populate the _ae and _an coefficients
   */
  ADReal computeStrongResidual();

  /// The dimension of the simulation
  const unsigned int _dim;

  /// index x|y|z
  const unsigned int _axis_index;

  /// x-velocity
  const INSFVVelocityVariable * const _u_var;
  /// y-velocity
  const INSFVVelocityVariable * const _v_var;
  /// z-velocity
  const INSFVVelocityVariable * const _w_var;

  /// Density
  const Moose::Functor<ADReal> & _rho;

  /// Turbulent eddy mixing length
  //const Moose::Functor<ADReal> & _mixing_len;
  const Moose::Functor<ADReal> & _eddy_visc_csv;

  /// Rhie-Chow element coefficient
  ADReal _ae = 0;

  /// Rhie-Chow neighbor coefficient
  ADReal _an = 0;
};
