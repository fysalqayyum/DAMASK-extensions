! SPDX-License-Identifier: AGPL-3.0-or-later
!--------------------------------------------------------------------------------------------------
!> @brief Local (per-cell) damage integration for the grid solver.
!
!> Reproduces DAMASK 2.0.3's `damage local` homogenization, which DAMASK 3 dropped.
!> DAMASK2's `damage_local_updateState` reads, in full:
!>
!>     phi = max(residualStiffness, min(1.0_pReal, phi + subdt*phiDot))
!>
!> i.e. explicit forward Euler on dphi/dt = f(phi), evaluated independently in every
!> cell, with no gradient term and no mobility. This module is that update and nothing
!> else. It deliberately does NOT solve
!>
!>     mu*dphi/dt = div(K grad phi) + f(phi)
!>
!> which is what grid_damage_spectral does. There is therefore no length scale, no
!> PETSc solve, and no FFT: cells never exchange information, so a cell reaching the
!> phi_min floor cannot make the damage operator singular. That is the whole point --
!> the nonlocal solver terminates a few increments after damage nucleates, because a
!> failed cell drives the Green operator singular, whereas the published DAMASK2 runs
!> ran to completion.
!>
!> The source term f(phi) is taken from homogenization_f_phi, i.e. the same
!> phase_f_phi that the spectral solver uses. phase_damage.f90 documents it as
!> "identical to DAMASK2 source_damage_isoDuctile_getRateAndItsTangent", so the
!> constitutive response here is DAMASK2's, unchanged.
!>
!> CAVEAT, to be stated in any publication using this module: a local damage model
!> with softening is mesh dependent. Refining the grid narrows the damage band toward
!> a single cell and drives the dissipated fracture energy toward zero. This module
!> exists to reproduce a published DAMASK2 study at its original resolution, not
!> because local damage is the better model.
!--------------------------------------------------------------------------------------------------
module grid_damage_local
#include <petsc/finclude/petscsys.h>
  use PETScSys
#ifndef PETSC_EXPOSES_MPI
  use MPI_f08
#endif

  use prec
  use parallelization
  use IO
  use misc
  use CLI
  use HDF5_utilities
  use HDF5
  use spectral_utilities
  use discretization_grid
  use homogenization
  use types
  use config

#ifndef PETSC_EXPOSES_MPIF90
  implicit none(type,external)
#else
  implicit none
#endif
  private

  type :: tNumerics
    real(pREAL) :: &
      phi_min, &                                                                                    !< non-zero residual damage (DAMASK2: residualStiffness)
      eps_damage_atol, &                                                                            !< absolute tolerance for staggered damage convergence
      eps_damage_rtol                                                                               !< relative tolerance for staggered damage convergence
    logical :: &
      irreversible                                                                                  !< forbid healing (phi <= phi_lastInc)
  end type tNumerics

  type(tNumerics) :: num

  real(pREAL), dimension(:,:,:), allocatable :: &
    phi, &                                                                                          !< current damage field
    phi_lastInc, &                                                                                  !< field of previous increment
    phi_stagInc                                                                                     !< field of previous staggered iteration

  public :: &
    grid_damage_local_init, &
    grid_damage_local_solution, &
    grid_damage_local_restartWrite, &
    grid_damage_local_forward

contains

!--------------------------------------------------------------------------------------------------
!> @brief Allocate all necessary fields and fill them with data, potentially from restart file.
!--------------------------------------------------------------------------------------------------
subroutine grid_damage_local_init(num_grid_damage)

  type(tDict), pointer, intent(in) :: num_grid_damage

  integer(HID_T) :: fileHandle, groupHandle
  real(pREAL), dimension(1,product(cells(1:2))*cells3) :: tempN
  character(len=:), allocatable :: extmsg


  print'(/,1x,a)', '<<<+-  grid_local_damage init  -+>>>'

  print'(/,1x,a)', 'local (per-cell) damage, reproducing DAMASK 2.0.3 `damage local`'
  print'(  1x,a)', 'no gradient term, no length scale -- see grid_damage_local.f90 header'

  if (.not. homogenization_damage_active()) call IO_error(501,ext_msg='damage')

!-------------------------------------------------------------------------------------------------
! read numerical parameters and do sanity checks
! Same dict path and key names as the spectral solver (solver -> grid -> damage) so a
! case can be switched between the two by editing only the load file. N_iter_max and
! the PETSc options are meaningless here and are ignored if present.
  num%eps_damage_atol = num_grid_damage%get_asReal('eps_abs_phi',defaultVal=1.0e-2_pREAL)
  num%eps_damage_rtol = num_grid_damage%get_asReal('eps_rel_phi',defaultVal=1.0e-6_pREAL)
  num%phi_min         = num_grid_damage%get_asReal('phi_min',    defaultVal=1.0e-6_pREAL)
  ! DAMASK2 clamped only to [residualStiffness, 1] and did NOT forbid healing, so the
  ! default here is .false. to match it. The spectral solver always enforces
  ! irreversibility (grid_damage_spectral.f90:370, min(...,phi_lastInc)).
  num%irreversible    = num_grid_damage%get_asBool('phi_irreversible',defaultVal=.false.)

  extmsg = ''
  if (num%eps_damage_atol <= 0.0_pREAL) extmsg = trim(extmsg)//' eps_abs_phi'
  if (num%eps_damage_rtol <= 0.0_pREAL) extmsg = trim(extmsg)//' eps_rel_phi'
  if (num%phi_min <= 0.0_pREAL .or. num%phi_min > 1.0_pREAL) &
                                           extmsg = trim(extmsg)//' phi_min'

  if (extmsg /= '') call IO_error(301,ext_msg=trim(extmsg))

!--------------------------------------------------------------------------------------------------
! init fields
  allocate(phi        (cells(1),cells(2),cells3), source=1.0_pREAL)
  allocate(phi_lastInc(cells(1),cells(2),cells3), source=1.0_pREAL)
  allocate(phi_stagInc(cells(1),cells(2),cells3), source=1.0_pREAL)

  restartRead: if (CLI_restartInc /= -1) then
    print'(/,1x,a,1x,i0)', 'loading restart data of increment', CLI_restartInc

    fileHandle  = HDF5_openFile(CLI_jobName//'_restart.hdf5','r')
    groupHandle = HDF5_openGroup(fileHandle,'solver')
    call HDF5_read(tempN,groupHandle,'phi',.false.)
    phi = reshape(tempN,[cells(1),cells(2),cells3])
    call HDF5_read(tempN,groupHandle,'phi_lastinc',.false.)
    phi_lastInc = reshape(tempN,[cells(1),cells(2),cells3])
    phi_stagInc = phi_lastInc
    call HDF5_closeGroup(groupHandle)
    call HDF5_closeFile(fileHandle)
  else
    phi = discretization_grid_getScalarInitialCondition('phi')
    phi_lastInc = phi
    phi_stagInc = phi_lastInc
  end if restartRead

  call homogenization_set_phi(reshape(phi,[product(cells(1:2))*cells3]))

end subroutine grid_damage_local_init


!--------------------------------------------------------------------------------------------------
!> @brief Integrate the local damage evolution over one time increment.
!
!> DAMASK2:  phi = max(residualStiffness, min(1.0, phi + subdt*phiDot))
!>
!> phiDot is evaluated at the current phi (not phi_lastInc) and the ODE is integrated
!> with adaptive sub-stepping: at each sub-step the largest dt that keeps every cell
!> above phi_min is used.  This prevents the catastrophic single-step crash observed
!> with high-N isoductile (phi^N blowup drives phi from ~0.7 to phi_min in one step,
!> degrading C66*phi**2 to ~1e-6 and killing the mechanical solver).  DAMASK 2's
!> damage_local had equivalent sub-stepping via its inner loop.
!--------------------------------------------------------------------------------------------------
function grid_damage_local_solution(Delta_t) result(solution)

  real(pREAL), intent(in) :: Delta_t                                                                !< increment in time for current solution

  type(tSolutionState) :: solution
  integer :: i, j, k, ce
  real(pREAL) :: phi_min_, phi_max_, stagNorm
  real(pREAL) :: t_remain, dt_sub, dt_sub_safe, floor_tol
  real(pREAL), dimension(cells(1),cells(2),cells3) :: phi_dot
  integer(MPI_INTEGER_KIND) :: err_MPI


!-----------------------------------------------------------------------------------------------
! Sub-step the forward-Euler update.  phi_dot = f(phi) = 1 - phi*D is evaluated at the
! current phi and integrated with the largest safe sub-step that keeps every cell above
! phi_min.  This is mathematically equivalent to DAMASK 2's damage_local inner loop.
!-----------------------------------------------------------------------------------------------
  phi = phi_lastInc                                                                               ! start from last-increment field
  t_remain = Delta_t
  floor_tol = 100.0_pREAL*epsilon(1.0_pREAL)
  do while (t_remain > 0.0_pREAL)

    ! evaluate phiDot = f(phi) at current phi
    ce = 0
    do k = 1, cells3;  do j = 1, cells(2);  do i = 1, cells(1)
      ce = ce + 1
      phi_dot(i,j,k) = homogenization_f_phi(phi(i,j,k),ce)
    end do; end do; end do

    ! A cell already at the residual-stiffness floor must not constrain the
    ! next safe step to zero.  Freeze only further degradation; retain a
    ! positive phi_dot so the historical reversible local model can heal.
    where (phi <= num%phi_min + floor_tol .and. phi_dot < 0.0_pREAL)
      phi_dot = 0.0_pREAL
    end where

    ! find largest sub-step that keeps every cell >= phi_min
    dt_sub = t_remain
    do k = 1, cells3;  do j = 1, cells(2);  do i = 1, cells(1)
      if (phi_dot(i,j,k) < 0.0_pREAL .and. phi(i,j,k) > num%phi_min + floor_tol) then
        dt_sub_safe = (num%phi_min - phi(i,j,k)) / phi_dot(i,j,k)
        if (dt_sub_safe < dt_sub) dt_sub = dt_sub_safe
      end if
    end do; end do; end do

    ! advance phi by dt_sub
    ce = 0
    do k = 1, cells3;  do j = 1, cells(2);  do i = 1, cells(1)
      ce = ce + 1
      phi(i,j,k) = phi(i,j,k) + dt_sub * phi_dot(i,j,k)
    end do; end do; end do

    phi = max(num%phi_min, min(1.0_pREAL, phi))
    where (phi <= num%phi_min + floor_tol) phi = num%phi_min
    if (dt_sub >= t_remain*(1.0_pREAL - floor_tol)) then
      t_remain = 0.0_pREAL
    else
      t_remain = t_remain - dt_sub
    end if
  end do

  if (num%irreversible) phi = min(phi,phi_lastInc)

  call homogenization_set_phi(reshape(phi,[product(cells(1:2))*cells3]))

!--------------------------------------------------------------------------------------------------
! There is no iterative solve, so the field solution itself never fails to converge --
! any non-convergence of the increment comes from the mechanical solver. The staggered
! check is still needed so the mechanical/damage staggered loop terminates.
  solution%converged        = .true.
  solution%iterationsNeeded = 1

  stagNorm = maxval(abs(phi - phi_stagInc))
  phi_min_ = minval(phi)
  phi_max_ = maxval(phi)
  call MPI_Allreduce(MPI_IN_PLACE,stagNorm,1_MPI_INTEGER_KIND,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)
  call MPI_Allreduce(MPI_IN_PLACE,phi_min_,1_MPI_INTEGER_KIND,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)
  call MPI_Allreduce(MPI_IN_PLACE,phi_max_,1_MPI_INTEGER_KIND,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)

  solution%stagConverged = stagNorm < max(num%eps_damage_atol, num%eps_damage_rtol*phi_max_)
  call MPI_Allreduce(MPI_IN_PLACE,solution%stagConverged,1_MPI_INTEGER_KIND,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)

  phi_stagInc = phi

  print'(/,1x,a)', '... local damage updated ..........................................'
  ! Identical format to grid_damage_spectral_solution: existing log post-processing
  ! greps "Delta Damage" and takes field 4 as the minimum phi.
  print'(/,1x,a,f8.6,2x,f8.6,2x,e11.4)', 'Minimum|Maximum|Delta Damage      = ', phi_min_, phi_max_, stagNorm
  print'(/,1x,a)', '==========================================================================='
  flush(IO_STDOUT)

end function grid_damage_local_solution


!--------------------------------------------------------------------------------------------------
!> @brief Set DAMASK data to current solver status.
!--------------------------------------------------------------------------------------------------
subroutine grid_damage_local_forward(cutBack)

  logical, intent(in) :: cutBack


  if (cutBack) then
    call homogenization_set_phi(reshape(phi_lastInc,[product(cells(1:2))*cells3]))
    phi         = phi_lastInc
    phi_stagInc = phi_lastInc
  else
    phi_lastInc = phi
  end if

end subroutine grid_damage_local_forward


!--------------------------------------------------------------------------------------------------
!> @brief Write current solver and constitutive data for restart to file.
!--------------------------------------------------------------------------------------------------
subroutine grid_damage_local_restartWrite()

  integer(HID_T) :: fileHandle, groupHandle


  print'(1x,a)', 'saving damage solver data required for restart'; flush(IO_STDOUT)

  fileHandle  = HDF5_openFile(CLI_jobName//'_restart.hdf5','a')
  groupHandle = HDF5_openGroup(fileHandle,'solver')
  call HDF5_write(reshape(phi,[1,product(cells(1:2))*cells3]),groupHandle,'phi')
  call HDF5_write(reshape(phi_lastInc,[1,product(shape(phi_lastInc))]),groupHandle,'phi_lastinc')
  call HDF5_closeGroup(groupHandle)
  call HDF5_closeFile(fileHandle)

end subroutine grid_damage_local_restartWrite

end module grid_damage_local
