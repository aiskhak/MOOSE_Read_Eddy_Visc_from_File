#include "INSFVEddyViscCSVReynoldsStress.h"
#include "INSFVVelocityVariable.h"
#include "NS.h"
#include "SystemBase.h"

registerMooseObject("NavierStokesApp", INSFVEddyViscCSVReynoldsStress);

InputParameters
INSFVEddyViscCSVReynoldsStress::validParams()
{
  InputParameters params = INSFVFluxKernel::validParams();
  params.addClassDescription("Computes force due to Reynolds stress term in incompressible RANS equations.");
  params.addRequiredCoupledVar("u", "Velocity in x direction.");
  params.addCoupledVar("v", "Velocity in y direction.");
  params.addCoupledVar("w", "Velocity in z direction.");
  params.addRequiredParam<MooseFunctorName>(NS::density, "fluid density");
  //params.addRequiredParam<MooseFunctorName>("mixing_length", "Turbulent eddy mixing length.");
  params.addRequiredParam<MooseFunctorName>("eddy_viscosity_csv", "Eddy viscosity from CSV file.");
  MooseEnum momentum_component("x=0 y=1 z=2");
  params.addRequiredParam<MooseEnum>("momentum_component", momentum_component, "Component of momentum equation that this kernel applies to.");
  params.set<unsigned short>("ghost_layers") = 3;
  return params;
}

INSFVEddyViscCSVReynoldsStress::INSFVEddyViscCSVReynoldsStress(const InputParameters & params)
  : INSFVFluxKernel(params),
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
INSFVEddyViscCSVReynoldsStress::computeStrongResidual()
{
#ifdef MOOSE_GLOBAL_AD_INDEXING
  
  // compute dot product of strain rate tensor and normal vector (grad_v+grad_v^T)*n_hat
  const auto & grad_u = _u_var->adGradSln(*_face_info);
  const ADRealVectorValue * grad_v = nullptr;
  const ADRealVectorValue * grad_w = nullptr;
  ADReal norm_strain_rate = grad_u(_axis_index) * _normal(0);
  if (_dim >= 2)
  {
    grad_v = &_v_var->adGradSln(*_face_info);
    norm_strain_rate += (*grad_v)(_axis_index)*_normal(1);
    if (_dim >= 3)
    {
      grad_w = &_w_var->adGradSln(*_face_info);
      norm_strain_rate += (*grad_w)(_axis_index)*_normal(2);
    }
  }
  const ADRealVectorValue & var_grad = _index == 0 ? grad_u : (_index == 1 ? *grad_v : *grad_w);
  norm_strain_rate += var_grad * _normal;
  
  // compute eddy diffusivity
  /*
  ADReal symmetric_strain_tensor_norm = 2.0 * Utility::pow<2>(grad_u(0));
  if (_dim >= 2)
  {
    symmetric_strain_tensor_norm += 2.0 * Utility::pow<2>((*grad_v)(1)) + 
	                                Utility::pow<2>((*grad_v)(0) + grad_u(1));
    if (_dim >= 3)
	{
      symmetric_strain_tensor_norm += 2.0 * Utility::pow<2>((*grad_w)(2)) +
                                      Utility::pow<2>(grad_u(2) + (*grad_w)(0)) +
                                      Utility::pow<2>((*grad_v)(2) + (*grad_w)(1));
	}
  }
  constexpr Real offset = 1e-15; // prevents explosion of sqrt(x) derivative to infinity
  symmetric_strain_tensor_norm = std::sqrt(symmetric_strain_tensor_norm + offset);
  // interpolate mixing length to face
  const auto face = Moose::FV::makeCDFace(*_face_info, faceArgSubdomains());
  const ADReal mixing_len = _mixing_len(face);
  // compute the eddy diffusivity
  ADReal eddy_diff = symmetric_strain_tensor_norm * mixing_len * mixing_len;
  */

  // interpolate eddy diffusivity to face
  const auto face = Moose::FV::makeCDFace(*_face_info, faceArgSubdomains());
  const ADReal eddy_diff = _eddy_visc_csv(face);

  const ADReal rho = _rho(face);

  if (_face_type == FaceInfo::VarFaceNeighbors::ELEM ||
      _face_type == FaceInfo::VarFaceNeighbors::BOTH)
  {
    const auto dof_number = _face_info->elem().dof_number(_sys.number(), _var.number(), 0);
    // norm_strain_rate is a linear combination of degrees of freedom so it's safe to straight-up
    // index into the derivatives vector at the dof we care about
    _ae = norm_strain_rate.derivatives()[dof_number];
    _ae *= -rho * eddy_diff;
  }
  if (_face_type == FaceInfo::VarFaceNeighbors::NEIGHBOR ||
      _face_type == FaceInfo::VarFaceNeighbors::BOTH)
  {
    const auto dof_number = _face_info->neighbor().dof_number(_sys.number(), _var.number(), 0);
    _an = norm_strain_rate.derivatives()[dof_number];
    _an *= rho * eddy_diff;
  }

  // Return the turbulent stress contribution to the momentum equation
  return -1 * rho * eddy_diff * norm_strain_rate;

#else
  return 0;

#endif
}

void
INSFVEddyViscCSVReynoldsStress::gatherRCData(const FaceInfo & fi)
{
  if (skipForBoundary(fi))
    return;

  _face_info = &fi;
  _normal = fi.normal();
  _face_type = fi.faceType(_var.name());

  processResidual(computeStrongResidual() * (fi.faceArea() * fi.faceCoord()));

  if (_face_type == FaceInfo::VarFaceNeighbors::ELEM ||
      _face_type == FaceInfo::VarFaceNeighbors::BOTH)
    _rc_uo.addToA(&fi.elem(), _index, _ae * (fi.faceArea() * fi.faceCoord()));
  if (_face_type == FaceInfo::VarFaceNeighbors::NEIGHBOR ||
      _face_type == FaceInfo::VarFaceNeighbors::BOTH)
    _rc_uo.addToA(fi.neighborPtr(), _index, _an * (fi.faceArea() * fi.faceCoord()));
}
