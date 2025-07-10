module mpas_microphysics

  use mpas_kind_types, only : RKIND
  use mpas_derived_types, only : MPAS_pool_type, MPAS_clock_type
  use mpas_pool_routines, only : mpas_pool_get_subpool, &
                                 mpas_pool_get_array, &
                                 mpas_pool_get_dimension
  use mpas_log, only : mpas_log_write

  implicit none

  ! Physical constants
  real(kind=RKIND), parameter :: gravity = 9.80616_RKIND

  ! Pools
  type(MPAS_pool_type), pointer :: mesh => null()
  type(MPAS_pool_type), pointer :: state => null()
  type(MPAS_pool_type), pointer :: diag => null()
  type(MPAS_pool_type), pointer :: diag_physics => null()

  ! Dimensions
  integer, pointer :: nCells => null()
  integer, pointer :: nVertLevels => null()

  ! Input arrays
  real(kind=RKIND), dimension(:,:), pointer :: pressure => null() 
  real(kind=RKIND), dimension(:,:), pointer :: prw_vcd => null() ! ! from microphysics diag_physics pool

  ! Accumulated 3D field for condensation_rate_acc
  real(kind=RKIND), dimension(:,:), pointer :: condensation_rate_acc => null() 

  ! Output arrays for diagnostics
  real(kind=RKIND), dimension(:,:), pointer :: condensation_rate => null()    ! instantaneous (equal to prw_vcd)
  real(kind=RKIND), dimension(:), pointer :: condensation_rate_vint => null()
  real(kind=RKIND), dimension(:), pointer :: condensation_rate_vint_acc => null()

  ! Internal state
  logical :: arrays_associated = .false.
  logical :: first_call = .true.

  public :: mpas_microphysics_setup, &
            mpas_microphysics_compute, &
            mpas_microphysics_reset, &
            mpas_microphysics_cleanup

contains


subroutine mpas_microphysics_setup(configs, structs, clock)
  use mpas_log, only : mpas_log_write
  implicit none

  type(MPAS_pool_type), pointer :: configs
  type(MPAS_pool_type), pointer :: structs
  type(MPAS_clock_type), pointer :: clock

  character(len=100) :: msg
  logical :: assoc

  ! Connect pools
  call mpas_pool_get_subpool(structs, 'mesh', mesh)
  call mpas_pool_get_subpool(structs, 'state', state)
  call mpas_pool_get_subpool(structs, 'diag', diag)
  call mpas_pool_get_subpool(structs, 'diag_physics', diag_physics)

  if (.not. associated(mesh)) then
    call mpas_log_write('FATAL ERROR: mesh pool not found!')
    stop 1
  end if

  if (.not. associated(state)) then
    call mpas_log_write('FATAL ERROR: state pool not found!')
    stop 1
  end if

  if (.not. associated(diag_physics)) then
    call mpas_log_write('FATAL ERROR: diag_physics pool not found!')
    stop 1
  end if

  if (.not. associated(diag)) then
    call mpas_log_write('WARNING: diag pool not found; output arrays may not be available')
  end if

  ! Get dimensions
  call mpas_pool_get_dimension(mesh, 'nCells', nCells)
  call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)

  ! Get input array
  call mpas_pool_get_array(diag_physics, 'prw_vcd', prw_vcd)
  if (.not. associated(prw_vcd)) then
    call mpas_log_write('FATAL ERROR: prw_vcd array not associated after retrieval.')
    stop 1
  else
    write(msg, '(A,I0,A,I0)') 'INFO: prw_vcd successfully associated. Size: ', &
         size(prw_vcd,1), ' x ', size(prw_vcd,2)
    call mpas_log_write(trim(msg))
  end if

  ! Get output arrays
  if (associated(diag)) then
    call mpas_pool_get_array(diag, 'condensation_rate', condensation_rate)
    call mpas_pool_get_array(diag, 'condensation_rate_vint', condensation_rate_vint)
    call mpas_pool_get_array(diag, 'condensation_rate_vint_acc', condensation_rate_vint_acc)
    call mpas_pool_get_array(diag, 'condensation_rate_acc', condensation_rate_acc)
    call mpas_pool_get_array(diag, 'pressure', pressure)

    assoc = associated(condensation_rate)
    write(msg, '(A,L1)') 'INFO: condensation_rate associated: ', assoc
    call mpas_log_write(trim(msg))

    assoc = associated(condensation_rate_acc)
    write(msg, '(A,L1)') 'INFO: condensation_rate_acc associated: ', assoc
    call mpas_log_write(trim(msg))

    assoc = associated(condensation_rate_vint)
    write(msg, '(A,L1)') 'INFO: condensation_rate_vint associated: ', assoc
    call mpas_log_write(trim(msg))

    assoc = associated(condensation_rate_vint_acc)
    write(msg, '(A,L1)') 'INFO: condensation_rate_vint_acc associated: ', assoc
    call mpas_log_write(trim(msg))

    assoc = associated(pressure)
    write(msg, '(A,L1)') 'INFO: pressure associated: ', assoc
    call mpas_log_write(trim(msg))
  end if

  ! Initialize arrays on first call
  if (first_call) then
    if (associated(condensation_rate))        condensation_rate        = 0.0_RKIND
    if (associated(condensation_rate_acc))    condensation_rate_acc    = 0.0_RKIND
    if (associated(condensation_rate_vint))   condensation_rate_vint   = 0.0_RKIND
    if (associated(condensation_rate_vint_acc)) condensation_rate_vint_acc = 0.0_RKIND

    first_call = .false.
    call mpas_log_write('INFO: Initialized condensation_rate arrays to zero.')
  end if

  arrays_associated = associated(prw_vcd) .and. associated(condensation_rate_acc) .and. &
                      associated(condensation_rate) .and. associated(condensation_rate_vint) .and. &
                      associated(condensation_rate_vint_acc) .and. associated(pressure)

  if (.not. arrays_associated) then
    call mpas_log_write('WARNING: Some condensation_rate arrays are not associated:')
    if (.not. associated(prw_vcd))               call mpas_log_write('  prw_vcd NOT associated')
    if (.not. associated(condensation_rate_acc)) call mpas_log_write('  condensation_rate_acc NOT associated')
    if (.not. associated(condensation_rate))     call mpas_log_write('  condensation_rate NOT associated')
    if (.not. associated(condensation_rate_vint))call mpas_log_write('  condensation_rate_vint NOT associated')
    if (.not. associated(condensation_rate_vint_acc)) call mpas_log_write('  condensation_rate_vint_acc NOT associated')
    if (.not. associated(pressure))              call mpas_log_write('  pressure NOT associated')
  end if

end subroutine mpas_microphysics_setup



subroutine mpas_microphysics_update()
  integer :: iCell, k
  real(kind=RKIND) :: p_top, p_bot, dp, layer_mass
  real(kind=RKIND) :: prw_val, vert_int_inst

  if (.not. arrays_associated) return

  do iCell = 1, nCells
    vert_int_inst = 0.0_RKIND

    do k = nVertLevels, 2, -1
      p_top = pressure(k, iCell)
      p_bot = pressure(k-1, iCell)
      dp = p_bot - p_top
      layer_mass = dp / gravity

      prw_val = max(0.0_RKIND, prw_vcd(k, iCell))

      ! Store instantaneous condensation_rate (instantaneous diagnostic)
      if (associated(condensation_rate)) condensation_rate(k, iCell) = prw_val

      ! Accumulate instantaneous condensation rate over time
      if (associated(condensation_rate_acc)) condensation_rate_acc(k, iCell) = &
           condensation_rate_acc(k, iCell) + prw_val

      ! Accumulate vertical integral of instantaneous condensation rate
      vert_int_inst = vert_int_inst + prw_val * layer_mass
    end do

    ! Accumulate vertical integral over time
    if (associated(condensation_rate_vint_acc)) condensation_rate_vint_acc(iCell) = &
         condensation_rate_vint_acc(iCell) + vert_int_inst
  end do
      
end subroutine mpas_microphysics_update


subroutine mpas_microphysics_compute()
  integer :: iCell, k
  real(kind=RKIND) :: p_top, p_bot, dp, layer_mass
  real(kind=RKIND) :: vert_int_inst
      
  character(len=256) :: message
  real(kind=RKIND) :: min_cr, max_cr
  real(kind=RKIND) :: min_vint, max_vint
  real(kind=RKIND) :: prw_val
     
  if (.not. arrays_associated) return

  do iCell = 1, nCells
    vert_int_inst = 0.0_RKIND

    ! Compute vertical integral of instantaneous condensation_rate
    do k = nVertLevels, 2, -1
      p_top = pressure(k, iCell)
      p_bot = pressure(k-1, iCell)
      dp = p_bot - p_top
      layer_mass = dp / gravity

      ! Clip negative values (since they are the evaporation rates) 
      prw_val = max(0.0_RKIND, condensation_rate(k, iCell))
      vert_int_inst = vert_int_inst + prw_val * layer_mass
    end do

    ! Store instantaneous vertical integral diagnostic
    if (associated(condensation_rate_vint)) condensation_rate_vint(iCell) = vert_int_inst
      end do
      
end subroutine mpas_microphysics_compute


subroutine mpas_microphysics_reset()
  use mpas_atm_diagnostics_utils, only : MPAS_field_will_be_written
  implicit none

  if (associated(condensation_rate_acc)) then
    if (MPAS_field_will_be_written('condensation_rate_acc')) then
      condensation_rate_acc(:,:) = 0.0_RKIND
    end if
  end if

  if (associated(condensation_rate_vint_acc)) then
    if (MPAS_field_will_be_written('condensation_rate_vint_acc')) then
      condensation_rate_vint_acc(:) = 0.0_RKIND
    end if
  end if
end subroutine mpas_microphysics_reset


subroutine mpas_microphysics_cleanup()
    nullify(mesh, state, diag, diag_physics)
    nullify(pressure, prw_vcd, condensation_rate_acc)
    nullify(condensation_rate, condensation_rate_vint, condensation_rate_vint_acc)
    nullify(nCells, nVertLevels)
    arrays_associated = .false.
end subroutine mpas_microphysics_cleanup

end module mpas_microphysics
