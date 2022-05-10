# global parameters
mu  = 1.e-3
rho = 1.e3
advected_interp_method = 'upwind'
velocity_interp_method = 'rc'
mesh = 'tamu_2d_7.msh'
v_in_av = 0.673

[GlobalParams]
  rhie_chow_user_object = 'rc'
[]

[UserObjects]
  [rc]
    type = INSFVRhieChowInterpolator
    u = u
    v = v
    pressure = pressure
  []
[]

[Mesh]
  type = FileMesh
  file = ${mesh}
  dim  = 2
[]

[Problem]
  fv_bcs_integrity_check = true
  coord_type = 'RZ'
  #restart_file_base = tamu_2d_fv_out_cp/LATEST
[]

[Variables]
  [u]
    type = INSFVVelocityVariable
    initial_condition = 0
  []
  [v]
    type = INSFVVelocityVariable
    initial_condition = 0
  []
  [pressure]
    type = INSFVPressureVariable
  []
[]

[AuxVariables]
  #[mixing_length_aux_var]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  fv = true
  #[]
  #[wall_shear_stress_aux_var]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  fv = true
  #[]
  #[wall_yplus_aux_var]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  fv = true
  #[]
  #[eddy_viscosity_aux_var]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  fv = true
  #[]
  [eddy_viscosity_csv_aux_var]
    order = CONSTANT
    family = MONOMIAL
    fv = true
  []
  [vel_aux_var]
    order = CONSTANT
    family = MONOMIAL
    fv = true
  []
  [dudx_aux_var]
    type = MooseVariableFVReal
  []
  [dudy_aux_var]
    type = MooseVariableFVReal
  []
  [dvdx_aux_var]
    type = MooseVariableFVReal
  []
  [dvdy_aux_var]
    type = MooseVariableFVReal
  []
  [elvol_aux_var]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Functions]
  [./v_in]
    type = ParsedFunction
    vars = v_av
    vals = ${v_in_av}
    value = 0.5*v_av*(1/7+1)*(1/7+2)*(1-x/9.525e-3)^(1/7) # turbulent profile
  [../]
[]

[FVKernels]
  [mass]
    type = INSFVMassAdvection
    variable = pressure
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    vel = 'velocity'
    pressure = pressure
    u = u
    v = v
    mu = ${mu}
    rho = ${rho}
  []

  [u_time]
    type = INSFVMomentumTimeDerivative
	momentum_component = 'x'
    variable = 'u'
    rho = ${rho}
  []
  [u_advection]
    type = INSFVMomentumAdvection
	momentum_component = 'x'
    variable = u
    advected_quantity = 'rhou'
    vel = 'velocity'
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    pressure = pressure
    u = u
    v = v
    mu = ${mu}
    rho = ${rho}
  []
  [u_viscosity]
    type = INSFVMomentumDiffusion
    variable = u
    mu = ${mu}
    momentum_component = 'x'
  []
  [u_viscosity_rans]
    type = INSFVEddyViscCSVReynoldsStress #INSFVMixingLengthReynoldsStress
    variable = u
    rho = ${rho}
    #mixing_length = mixing_len
    eddy_viscosity_csv = eddy_viscosity_csv_aux_var
    momentum_component = 'x'
    u = u
    v = v
  []
  [u_pressure]
    type = INSFVMomentumPressure
    variable = u
    momentum_component = 'x'
    pressure = pressure
  []

  [v_time]
    type = INSFVMomentumTimeDerivative
	momentum_component = 'y'
    variable = v
    rho = ${rho}
  []
  [v_advection]
    type = INSFVMomentumAdvection
	momentum_component = 'y'
    variable = v
    advected_quantity = 'rhov'
    vel = 'velocity'
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    pressure = pressure
    u = u
    v = v
    mu = ${mu}
    rho = ${rho}
  []
  [v_viscosity]
    type = INSFVMomentumDiffusion
    variable = v
    mu = ${mu}
    momentum_component = 'y'
  []
  [v_viscosity_rans]
    type = INSFVEddyViscCSVReynoldsStress #INSFVMixingLengthReynoldsStress
    variable = v
    rho = ${rho}
    #mixing_length = mixing_len
    eddy_viscosity_csv = eddy_viscosity_csv_aux_var
    momentum_component = 'y'
    u = u
    v = v
  []
  [v_pressure]
    type = INSFVMomentumPressure
    variable = v
    momentum_component = 'y'
    pressure = pressure
  []
[]

[AuxKernels]
  #[mixing_len_aux_ker]
  #  type = WallDistanceMixingLengthAux
  #  walls = 'wall'
  #  variable = mixing_length_aux_var
  #  execute_on = 'initial'
  #  von_karman_const = 0.41
  #  delta = 1.e12 #9.525e-3
  #[]
  #[eddy_viscosity_aux_ker]
  #  type = INSFVMixingLengthTurbulentViscosityAux
  #  variable = eddy_viscosity_aux_var
  #  mixing_length = mixing_length_aux_var
  #  u = u
  #  v = v
  #[]
  [vel_mag_aux_ker]
    type = VectorMagnitudeAux
    variable = vel_aux_var
    x = u
    y = v
  []
  [eddy_viscosity_aux_ker]
    type = AuxVarFromCSVFile
    variable = eddy_viscosity_csv_aux_var
    file_name = 'turb_visc_nek.csv'
    header = true
  []
  [dudx_aux_ker]
    type = ADFunctorVectorElementalAux
    variable = 'dudx_aux_var'
    functor = 'grad_u'
    component = 0
  []
  [dudy_aux_ker]
    type = ADFunctorVectorElementalAux
    variable = 'dudy_aux_var'
    functor = 'grad_u'
    component = 1
  []
  [dvdx_aux_ker]
    type = ADFunctorVectorElementalAux
    variable = 'dvdx_aux_var'
    functor = 'grad_v'
    component = 0
  []
  [dvdy_aux_ker]
    type = ADFunctorVectorElementalAux
    variable = 'dvdy_aux_var'
    functor = 'grad_v'
    component = 1
  []
  [elvol_aux_ker]
    type = VolumeAux
    variable = elvol_aux_var
  []
[]

[FVBCs]
  [inlet-u]
    type = INSFVInletVelocityBC
    boundary = 'inlet'
    variable = u
    function = 0
  []
  [inlet-v]
    type = INSFVInletVelocityBC
    boundary = 'inlet'
    variable = v
    function = 'v_in'
  []
  [no-slip-wall-u]
    type = INSFVNoSlipWallBC
    boundary = 'wall'
    variable = u
    function = 0
  []
  [no-slip-wall-v]
    type = INSFVNoSlipWallBC
    boundary = 'wall'
    variable = v
    function = 0
  []
  [outlet-p]
    type = INSFVOutletPressureBC
    boundary = 'outlet'
    variable = pressure
    function = 0
  []
  [axis-u]
    type = INSFVSymmetryVelocityBC
    boundary = 'SYM'
    variable = u
    u = u
    v = v
    mu = ${mu}
    momentum_component = x
  []
  [axis-v]
    type = INSFVSymmetryVelocityBC
    boundary = 'SYM'
    variable = v
    u = u
    v = v
    mu = ${mu}
    momentum_component = y
  []
  [axis-p]
    type = INSFVSymmetryPressureBC
    boundary = 'SYM'
    variable = pressure
  []
[]

[Materials]
  #[ins_fv]
  #  type = INSFVMaterial
  #  u = 'u'
  #  v = 'v'
  #  pressure = 'pressure'
  #  rho = ${rho}
  #[]
  [grad_u_mat]
    type = ADGenericFunctorGradientMaterial
    prop_names = 'grad_u'
    prop_values = 'u'
  []
  [grad_v_mat]
    type = ADGenericFunctorGradientMaterial
    prop_names = 'grad_v'
    prop_values = 'v'
  []
[]

[Postprocessors]
  [in]
    type = SideIntegralVariablePostprocessor
    variable = v
    boundary = 'inlet'
  []
  [out]
    type = SideIntegralVariablePostprocessor
    variable = v
    boundary = 'outlet'
  []
[]

[VectorPostprocessors]
  [vpp]
    type = ElementValueSampler
    variable = 'u v pressure eddy_viscosity_csv_aux_var dudx_aux_var dudy_aux_var dvdx_aux_var dvdy_aux_var'
    sort_by = x
  []
  [elv]
    type = ElementValueSampler
    variable = 'elvol_aux_var'
    sort_by = x
    execute_on = 'initial'
  []
[]

[Preconditioning]
  [./SMP_PJFNK]
    type = SMP
    full = true
    solve_type = 'PJFNK'
    petsc_options_iname = '-pc_type -ksp_gmres_restart'
    petsc_options_value = 'lu 100'
  [../]
[]

[Executioner]
  type = Transient
  [./TimeStepper]
    type = IterationAdaptiveDT
    growth_factor = 1.5
    optimal_iterations = 8
    linear_iteration_ratio = 150
    dt = 5e-1
    cutback_factor = 0.8
    cutback_factor_at_failure = 0.8
  [../]
 #dt = 0.00001
 dtmin = 1e-6
 dtmax = 2
 nl_rel_tol = 1e-6
 nl_abs_tol = 1e-6
 nl_max_its = 20
 l_tol = 1e-5
 l_max_its = 100
 start_time = 0
 end_time = 4000
 num_steps = 200
 steady_state_detection = true
 steady_state_tolerance = 1.e-8
[]

[Outputs]
  exodus = true
  csv = true
  checkpoint = true
[]
