module mpas_microphysics
  use mpas_kind_types, only : RKIND
  use mpas_derived_types, only : MPAS_pool_type, MPAS_clock_type
  use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_get_array, mpas_pool_get_dimension
  use mpas_log, only : mpas_log_write

  implicit none

  real(kind=RKIND), parameter :: gravity = 9.80616_RKIND

  type(MPAS_pool_type), pointer :: mesh => null()
  type(MPAS_pool_type), pointer :: state => null()
  type(MPAS_pool_type), pointer :: diag => null()
  type(MPAS_pool_type), pointer :: diag_physics => null()

  integer, pointer :: nCells => null()
  integer, pointer :: nVertLevels => null()

  ! input: prw_vcd is the total condensation/ evaporation 
  real(kind=RKIND), dimension(:,:), pointer :: pressure => null()
  real(kind=RKIND), dimension(:,:), pointer :: prw_vcd => null()
      
  ! deposition process rates
  real(kind=RKIND), dimension(:,:), pointer :: prs_ide => null(), prs_sde => null()
  real(kind=RKIND), dimension(:,:), pointer :: prg_gde => null(), pri_inu => null()
  real(kind=RKIND), dimension(:,:), pointer :: pri_iha => null(), pri_ide => null()

! diagnostics
  real(kind=RKIND), dimension(:),   pointer :: condensation_vint => null()
  real(kind=RKIND), dimension(:),   pointer :: condensation_vint_acc => null()
  real(kind=RKIND), dimension(:), pointer :: evaporation_vint => null()
  real(kind=RKIND), dimension(:), pointer :: evaporation_vint_acc => null()
  real(kind=RKIND), dimension(:), pointer :: deposition_vint => null()
  real(kind=RKIND), dimension(:), pointer :: deposition_vint_acc => null()
  real(kind=RKIND), dimension(:), pointer :: sublimation_vint => null()
  real(kind=RKIND), dimension(:), pointer :: sublimation_vint_acc => null()

  logical :: arrays_associated = .false.
  logical :: first_call = .true.

  public :: mpas_microphysics_setup, mpas_microphysics_compute, mpas_microphysics_update, &
            mpas_microphysics_reset, mpas_microphysics_cleanup

contains

  subroutine mpas_microphysics_setup(configs, structs, clock)
    use mpas_log, only : mpas_log_write
    implicit none
    type(MPAS_pool_type), pointer :: configs, structs
    type(MPAS_clock_type), pointer :: clock
    character(len=200) :: msg
    logical :: assoc

    call mpas_pool_get_subpool(structs, 'mesh', mesh)
    call mpas_pool_get_subpool(structs, 'state', state)
    call mpas_pool_get_subpool(structs, 'diag', diag)
    call mpas_pool_get_subpool(structs, 'diag_physics', diag_physics)

    if (.not. associated(mesh)) then
      call mpas_log_write('FATAL ERROR: mesh pool not found!'); stop 1
    end if
    if (.not. associated(state)) then
      call mpas_log_write('FATAL ERROR: state pool not found!'); stop 1
    end if
    if (.not. associated(diag_physics)) then
      call mpas_log_write('FATAL ERROR: diag_physics pool not found!'); stop 1
    end if
    if (.not. associated(diag)) then
      call mpas_log_write('WARNING: diag pool not found; some output arrays will be missing')
    end if

    call mpas_pool_get_dimension(mesh, 'nCells', nCells)
    call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)

    call mpas_pool_get_array(diag_physics, 'prw_vcd', prw_vcd)
    call mpas_pool_get_array(diag_physics, 'prs_ide', prs_ide)
    call mpas_pool_get_array(diag_physics, 'prs_sde', prs_sde)
    call mpas_pool_get_array(diag_physics, 'prg_gde', prg_gde)
    call mpas_pool_get_array(diag_physics, 'pri_inu', pri_inu)
    call mpas_pool_get_array(diag_physics, 'pri_iha', pri_iha)
    call mpas_pool_get_array(diag_physics, 'pri_ide', pri_ide)

    if (associated(diag)) then
      call mpas_pool_get_array(diag, 'condensation_vint', condensation_vint)
      call mpas_pool_get_array(diag, 'condensation_vint_acc', condensation_vint_acc)
      call mpas_pool_get_array(diag, 'evaporation_vint', evaporation_vint)
      call mpas_pool_get_array(diag, 'evaporation_vint_acc', evaporation_vint_acc)
      call mpas_pool_get_array(diag, 'deposition_vint', deposition_vint)
      call mpas_pool_get_array(diag, 'deposition_vint_acc', deposition_vint_acc)
      call mpas_pool_get_array(diag, 'sublimation_vint', sublimation_vint)
      call mpas_pool_get_array(diag, 'sublimation_vint_acc', sublimation_vint_acc)
      call mpas_pool_get_array(diag, 'pressure', pressure)
    end if

      if (first_call) then
      if (associated(condensation_vint)) condensation_vint(:) = 0.0_RKIND
      if (associated(sublimation_vint)) sublimation_vint(:) = 0.0_RKIND
      if (associated(deposition_vint)) deposition_vint(:) = 0.0_RKIND
      if (associated(evaporation_vint)) evaporation_vint(:) = 0.0_RKIND
      if (associated(condensation_vint_acc)) condensation_vint_acc(:) = 0.0_RKIND
      if (associated(evaporation_vint_acc)) evaporation_vint_acc(:) = 0.0_RKIND
      if (associated(deposition_vint_acc)) deposition_vint_acc(:) = 0.0_RKIND
      if (associated(sublimation_vint_acc)) sublimation_vint_acc(:) = 0.0_RKIND
      call mpas_log_write('INFO: Initialized diagnostics for microphysics to zero.')
      first_call = .false.
    end if

    ! check if we actually associate these arrays 
    arrays_associated = associated(prw_vcd) .and. associated(pressure) .and. associated(evaporation_vint_acc) &
      .and. associated(deposition_vint_acc) .and. associated(sublimation_vint_acc) &
      .and. associated(condensation_vint_acc)

  end subroutine mpas_microphysics_setup


subroutine mpas_microphysics_update()
  integer :: iCell, k
  real(kind=RKIND) :: p_top, p_bot, dp, layer_mass
  real(kind=RKIND) :: condensation, evaporation
  real(kind=RKIND) :: sum_dep, sum_sub
  real(kind=RKIND) :: vert_inst_cr, vert_inst_evap, vert_inst_dep, vert_inst_sub

  if (.not. arrays_associated) return

  do iCell = 1, nCells
    vert_inst_cr    = 0.0_RKIND
    vert_inst_evap  = 0.0_RKIND
    vert_inst_dep   = 0.0_RKIND
    vert_inst_sub   = 0.0_RKIND

    do k = nVertLevels, 2, -1
      p_top = pressure(k, iCell)
      p_bot = pressure(k-1, iCell)
      dp = p_bot - p_top
      layer_mass = dp / gravity

      condensation = max(0.0_RKIND, prw_vcd(k, iCell))
      evaporation = -min(0.0_RKIND, prw_vcd(k, iCell))

      ! --- Condensation ---
      vert_inst_cr = vert_inst_cr + condensation * layer_mass
      
      ! --- Evaporation ---
      vert_inst_evap = vert_inst_evap + evaporation * layer_mass
      
      ! --- Deposition and Sublimation ---
      sum_dep = 0.0_RKIND
      sum_sub = 0.0_RKIND

      if (prs_ide(k,iCell) > 0.0_RKIND) sum_dep = sum_dep + prs_ide(k,iCell)
      if (prs_ide(k,iCell) < 0.0_RKIND) sum_sub = sum_sub - prs_ide(k,iCell)

      if (prs_sde(k,iCell) > 0.0_RKIND) sum_dep = sum_dep + prs_sde(k,iCell)
      if (prs_sde(k,iCell) < 0.0_RKIND) sum_sub = sum_sub - prs_sde(k,iCell)

      if (prg_gde(k,iCell) > 0.0_RKIND) sum_dep = sum_dep + prg_gde(k,iCell)
      if (prg_gde(k,iCell) < 0.0_RKIND) sum_sub = sum_sub - prg_gde(k,iCell)

      if (pri_inu(k,iCell) > 0.0_RKIND) sum_dep = sum_dep + pri_inu(k,iCell)
      if (pri_inu(k,iCell) < 0.0_RKIND) sum_sub = sum_sub - pri_inu(k,iCell)

      if (pri_iha(k,iCell) > 0.0_RKIND) sum_dep = sum_dep + pri_iha(k,iCell)
      if (pri_iha(k,iCell) < 0.0_RKIND) sum_sub = sum_sub - pri_iha(k,iCell)

      if (pri_ide(k,iCell) > 0.0_RKIND) sum_dep = sum_dep + pri_ide(k,iCell)
      if (pri_ide(k,iCell) < 0.0_RKIND) sum_sub = sum_sub - pri_ide(k,iCell)

      vert_inst_dep = vert_inst_dep + sum_dep * layer_mass
      vert_inst_sub = vert_inst_sub + sum_sub * layer_mass
    end do

    if (associated(condensation_vint_acc)) condensation_vint_acc(iCell) = &
         condensation_vint_acc(iCell) + vert_inst_cr
    if (associated(evaporation_vint_acc)) evaporation_vint_acc(iCell) = &
         evaporation_vint_acc(iCell) + vert_inst_evap
    if (associated(deposition_vint_acc)) deposition_vint_acc(iCell) = &
         deposition_vint_acc(iCell) + vert_inst_dep
    if (associated(sublimation_vint_acc)) sublimation_vint_acc(iCell) = &
         sublimation_vint_acc(iCell) + vert_inst_sub
  end do
end subroutine mpas_microphysics_update


subroutine mpas_microphysics_compute()
  integer :: iCell, k
  real(kind=RKIND) :: p_top, p_bot, dp, layer_mass
  real(kind=RKIND) :: vert_inst_cr, vert_inst_evap, vert_inst_dep, vert_inst_sub
  real(kind=RKIND) :: condensation, evaporation, sum_dep, sum_sub

  if (.not. arrays_associated) return

  do iCell = 1, nCells
    vert_inst_cr    = 0.0_RKIND
    vert_inst_evap  = 0.0_RKIND
    vert_inst_dep   = 0.0_RKIND
    vert_inst_sub   = 0.0_RKIND

    do k = nVertLevels, 2, -1
      p_top = pressure(k, iCell)
      p_bot = pressure(k-1, iCell)
      dp = p_bot - p_top
      layer_mass = dp / gravity

      condensation = max(0.0_RKIND, prw_vcd(k, iCell))
      evaporation = -min(0.0_RKIND, prw_vcd(k, iCell))
      vert_inst_cr = vert_inst_cr + condensation * layer_mass
      vert_inst_evap = vert_inst_evap + evaporation * layer_mass
     
      ! --- Deposition and Sublimation from six process rates ---
      sum_dep = 0.0_RKIND
      sum_sub = 0.0_RKIND

      if (prs_ide(k,iCell) > 0.0_RKIND) sum_dep = sum_dep + prs_ide(k,iCell)
      if (prs_ide(k,iCell) < 0.0_RKIND) sum_sub = sum_sub - prs_ide(k,iCell)

      if (prs_sde(k,iCell) > 0.0_RKIND) sum_dep = sum_dep + prs_sde(k,iCell)
      if (prs_sde(k,iCell) < 0.0_RKIND) sum_sub = sum_sub - prs_sde(k,iCell)

      if (prg_gde(k,iCell) > 0.0_RKIND) sum_dep = sum_dep + prg_gde(k,iCell)
      if (prg_gde(k,iCell) < 0.0_RKIND) sum_sub = sum_sub - prg_gde(k,iCell)

      if (pri_inu(k,iCell) > 0.0_RKIND) sum_dep = sum_dep + pri_inu(k,iCell)
      if (pri_inu(k,iCell) < 0.0_RKIND) sum_sub = sum_sub - pri_inu(k,iCell)

      if (pri_iha(k,iCell) > 0.0_RKIND) sum_dep = sum_dep + pri_iha(k,iCell)
      if (pri_iha(k,iCell) < 0.0_RKIND) sum_sub = sum_sub - pri_iha(k,iCell)

      if (pri_ide(k,iCell) > 0.0_RKIND) sum_dep = sum_dep + pri_ide(k,iCell)
      if (pri_ide(k,iCell) < 0.0_RKIND) sum_sub = sum_sub - pri_ide(k,iCell)

      vert_inst_dep = vert_inst_dep + sum_dep * layer_mass
      vert_inst_sub = vert_inst_sub + sum_sub * layer_mass
    end do

    ! Store vertical integrals for inst value 
    if (associated(condensation_vint)) condensation_vint(iCell) = vert_inst_cr
    if (associated(evaporation_vint))  evaporation_vint(iCell)  = vert_inst_evap
    if (associated(deposition_vint))        deposition_vint(iCell)        = vert_inst_dep
    if (associated(sublimation_vint))       sublimation_vint(iCell)       = vert_inst_sub
  end do
end subroutine mpas_microphysics_compute

      
  subroutine mpas_microphysics_reset()
    use mpas_atm_diagnostics_utils, only : MPAS_field_will_be_written
    implicit none

    call mpas_log_write('INFO: Resetting accumulators for extended diagnostics')
    if (associated(condensation_vint_acc)) then
      if (MPAS_field_will_be_written('condensation_vint_acc')) condensation_vint_acc(:) = 0.0_RKIND
    end if
    if (associated(evaporation_vint_acc)) then
      if (MPAS_field_will_be_written('evaporation_vint_acc')) evaporation_vint_acc(:) = 0.0_RKIND
    end if
    if (associated(deposition_vint_acc)) then
      if (MPAS_field_will_be_written('deposition_vint_acc')) deposition_vint_acc(:) = 0.0_RKIND
    end if
    if (associated(sublimation_vint_acc)) then
      if (MPAS_field_will_be_written('sublimation_vint_acc')) sublimation_vint_acc(:) = 0.0_RKIND
    end if
  end subroutine mpas_microphysics_reset

      
 subroutine mpas_microphysics_cleanup()
    nullify(mesh, state, diag, diag_physics)
    nullify(pressure, prw_vcd, prs_ide, prs_sde, prg_gde, pri_inu, pri_iha, pri_ide)
    nullify(condensation_vint, condensation_vint_acc)
    nullify(evaporation_vint, evaporation_vint_acc, deposition_vint, deposition_vint_acc, sublimation_vint, sublimation_vint_acc)
    nullify(nCells, nVertLevels)
    arrays_associated = .false.
  end subroutine mpas_microphysics_cleanup

end module mpas_microphysics
