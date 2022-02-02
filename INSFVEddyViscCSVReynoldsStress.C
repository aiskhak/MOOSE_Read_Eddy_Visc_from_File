#include "INSFVEddyViscCSVReynoldsStress.h"
#include "INSFVVelocityVariable.h"
#include "NS.h"

registerMooseObject("MooseApp", INSFVEddyViscCSVReynoldsStress);

InputParameters
INSFVEddyViscCSVReynoldsStress::validParams()
{
  InputParameters params = FVFluxKernel::validParams();
  params.addClassDescription("Computes force due to Reynolds stress term in incompressible RANS equations.");
  params.addRequiredCoupledVar("u", "Velocity in x direction.");
  params.addCoupledVar("v", "Velocity in y direction.");
  params.addCoupledVar("w", "Velocity in z direction.");
  params.addRequiredParam<MooseFunctorName>(NS::density, "fluid density");
  //params.addRequiredParam<MooseFunctorName>("mixing_length", "Turbulent eddy mixing length.");
  params.addRequiredParam<MooseFunctorName>("eddy_viscosity_csv", "Eddy viscosity from CSV file.");
  MooseEnum momentum_component("x=0 y=1 z=2");
  params.addRequiredParam<MooseEnum>("momentum_component", momentum_component, "Component of momentum equation that this kernel applies to.");
  params.set<unsigned short>("ghost_layers") = 2;
  return params;
}

INSFVEddyViscCSVReynoldsStress::INSFVEddyViscCSVReynoldsStress(const InputParameters & params)
  : FVFluxKernel(params),
    _dim(_subproblem.mesh().dimension()),
    _axis_index(getParam<MooseEnum>("momentum_component")),
    _u_var(dynamic_cast<const INSFVVelocityVariable *>(getFieldVar("u", 0))),
    _v_var(params.isParamValid("v") ? dynamic_cast<const INSFVVelocityVariable *>(getFieldVar("v", 0)) : nullptr),
    _w_var(params.isParamValid("w") ? dynamic_cast<const INSFVVelocityVariable *>(getFieldVar("w", 0)) : nullptr),
    _rho(getFunctor<ADReal>(NS::density)),
    //_mixing_len(getFunctor<ADReal>("mixing_length"))
    _eddy_visc_csv(getFunctor<ADReal>("eddy_viscosity_csv"))
{
#ifndef MOOSE_GLOBAL_AD_INDEXING
  mooseError("INSFV is not supported by local AD indexing. In order to use INSFV, please run the "
             "configure script in the root MOOSE directory with the configure option "
             "'--with-ad-indexing-type=global'");
#endif

  if (!_u_var)
    paramError("u", "u velocity must be an INSFVVelocityVariable.");

  if (_dim >= 2 && !_v_var)
    paramError("v", "In 2 or more dimensions, v velocity must be supplied and it must be an INSFVVelocityVariable.");

  if (_dim >= 3 && !_w_var)
    paramError("w", "In 3D, w velocity must be supplied and it must be an INSFVVelocityVariable.");
}

ADReal
INSFVEddyViscCSVReynoldsStress::computeQpResidual()
{
#ifdef MOOSE_GLOBAL_AD_INDEXING
  /*
  constexpr Real offset = 1e-15; // prevents explosion of sqrt(x) derivative to infinity

  const auto & grad_u = _u_var->adGradSln(*_face_info);
  ADReal symmetric_strain_tensor_norm = 2.0 * Utility::pow<2>(grad_u(0));
  if (_dim >= 2)
  {
    const auto & grad_v = _v_var->adGradSln(*_face_info);
    symmetric_strain_tensor_norm +=
        2.0 * Utility::pow<2>(grad_v(1)) + Utility::pow<2>(grad_v(0) + grad_u(1));
    if (_dim >= 3)
    {
      const auto & grad_w = _w_var->adGradSln(*_face_info);
      symmetric_strain_tensor_norm += 2.0 * Utility::pow<2>(grad_w(2)) +
                                      Utility::pow<2>(grad_u(2) + grad_w(0)) +
                                      Utility::pow<2>(grad_v(2) + grad_w(1));
    }
  }

  symmetric_strain_tensor_norm = std::sqrt(symmetric_strain_tensor_norm + offset);
  */
  // Interpolate the mixing length to the face
  //const ADReal mixing_len = _mixing_len(std::make_tuple(
  //    _face_info, Moose::FV::LimiterType::CentralDifference, true, faceArgSubdomains(_face_info)));

  // Compute the eddy diffusivity
  //ADReal eddy_diff = symmetric_strain_tensor_norm * mixing_len * mixing_len;

  // Interpolate the eddy viscosity to the face
  const ADReal eddy_diff = _eddy_visc_csv(std::make_tuple(
      _face_info, Moose::FV::LimiterType::CentralDifference, true, faceArgSubdomains(_face_info)));
  // Compute the dot product of the strain rate tensor and the normal vector
  // aka (grad_v + grad_v^T) * n_hat
  ADReal norm_strain_rate = gradUDotNormal();
  norm_strain_rate += _u_var->adGradSln(*_face_info)(_axis_index) * _normal(0);
  norm_strain_rate += _dim >= 2 ? _v_var->adGradSln(*_face_info)(_axis_index) * _normal(1) : 0;
  norm_strain_rate += _dim >= 3 ? _w_var->adGradSln(*_face_info)(_axis_index) * _normal(2) : 0;

  const ADReal rho = _rho(std::make_tuple(
      _face_info, Moose::FV::LimiterType::CentralDifference, true, faceArgSubdomains(_face_info)));

  // Return the turbulent stress contribution to the momentum equation
  return -1 * rho * eddy_diff * norm_strain_rate;

#else
  return 0;

#endif
}
