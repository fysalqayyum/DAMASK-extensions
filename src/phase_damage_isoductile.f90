! SPDX-License-Identifier: AGPL-3.0-or-later
!--------------------------------------------------------------------------------------------------
!> @brief material subroutine incorporating isotropic ductile damage source mechanism
!> @details Ported to DAMASK 3.1.0 from the DAMASK v2.0.3 grid-solver source
!!          "source_damage_isoDuctile.f90" (P. Shanthraj, L. Sharma, Max-Planck-Institut fuer
!!          Eisenforschung GmbH). DAMASK 3.1.0 upstream only ships isobrittle/anisobrittle damage;
!!          this restores the plastic-strain-driven ductile nucleation criterion used by DAMASK 2
!!          material.config files (source damage_isoductile), unchanged in formulation.
!!
!!          Local damage driving force D accumulates as
!!            dot(D) = dot_gamma_sum / phi**N / gamma_crit
!!          where dot_gamma_sum is the sum of the absolute shear rates over all active slip
!!          (and, if present, twin) systems, N is a rate-sensitivity exponent, and gamma_crit is
!!          the critical accumulated plastic shear. The local part of the nonlocal damage source
!!          term is f = 1 - D*phi, identical to DAMASK 2's source_damage_isoDuctile_getRateAndItsTangent.
!--------------------------------------------------------------------------------------------------
submodule (phase:damage) isoductile

  type :: tParameters                                                                               !< container type for internal constitutive parameters
    real(pREAL) :: &
      gamma_crit, &                                                                                 !< critical accumulated plastic shear (DAMASK2: isoductile_criticalplasticstrain)
      N                                                                                             !< damage rate sensitivity exponent (DAMASK2: isoductile_ratesensitivity)
    character(len=pSTRLEN), allocatable, dimension(:) :: &
      output
  end type tParameters

  type(tParameters), dimension(:), allocatable :: param                                             !< containers of constitutive parameters (len Ninstances)

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
module function isoductile_init() result(mySources)

  logical, dimension(:), allocatable :: mySources

  type(tDict), pointer :: &
    phases, &
    phase, &
    src
  integer :: Nmembers,ph
  character(len=:), allocatable :: &
    refs, &
    extmsg


  mySources = source_active('isoductile')
  if (count(mySources) == 0) return

  print'(/,1x,a)', '<<<+-  phase:damage:isoductile init  -+>>>'
  print'(/,1x,a,1x,i0)', '# phases:',count(mySources); flush(IO_STDOUT)


  phases => config_material%get_dict('phase')
  allocate(param(size(phases)))
  extmsg = ''

  do ph = 1, size(phases)
    if (mySources(ph)) then
      phase => phases%get_dict(ph)
      src => phase%get_dict('damage')

      associate(prm => param(ph))

        print'(/,1x,a,1x,i0,a)', 'phase',ph,': '//phases%key(ph)
        refs = config_listReferences(src,indent=3)
        if (len(refs) > 0) print'(/,1x,a)', refs

        prm%N          = src%get_asReal('N')
        prm%gamma_crit = src%get_asReal('gamma_crit')

#if defined (__GFORTRAN__)
        prm%output = output_as1dStr(src)
#else
        prm%output = src%get_as1dStr('output',defaultVal=emptyStrArray)
#endif

        if (prm%N          <= 0.0_pREAL) extmsg = trim(extmsg)//' N'
        if (prm%gamma_crit <= 0.0_pREAL) extmsg = trim(extmsg)//' gamma_crit'

        Nmembers = count(material_ID_phase==ph)
        call phase_allocateState(damageState(ph),Nmembers,1,1,0)
        damageState(ph)%atol = src%get_asReal('atol_phi',defaultVal=1.0e-3_pREAL)
        if (any(damageState(ph)%atol < 0.0_pREAL)) extmsg = trim(extmsg)//' atol_phi'

      end associate

      if (extmsg /= '') call IO_error(211,ext_msg=trim(extmsg)//'(damage_isoDuctile)')
    end if

  end do

end function isoductile_init


!--------------------------------------------------------------------------------------------------
!> @brief calculates the rate of the local damage driving force state variable
!> @details identical formulation to DAMASK2 source_damage_isoDuctile_dotState:
!!          dot(D) = sum(|shear rates|) / phi**N / gamma_crit
!!          Guard: phi is clamped below at phi_clamp to prevent the phi**N singularity
!!          from driving dot(D) → ∞ when phi → 0.  phi_clamp is set to the solver's
!!          phi_min default (1e-3) which is the physical floor for the stiffness
!!          degradation.  Without this guard, phi drops from ~0.7 to phi_min in a single
!!          increment, degrading the elastic stiffness by phi**2 ≈ 1e-6 and crashing the
!!          mechanical solver (error 950).
!--------------------------------------------------------------------------------------------------
module subroutine isoductile_dotState(ph, en)

  integer, intent(in) :: &
    ph,en

  real(pREAL), parameter :: &
    phi_clamp = 1.0e-3_pREAL


  associate(prm => param(ph))
    damageState(ph)%dotState(1,en) = &
      plastic_dotGammaSum(ph,en) / max(damage_phi(ph,en), phi_clamp)**prm%N / prm%gamma_crit
  end associate

end subroutine isoductile_dotState


!--------------------------------------------------------------------------------------------------
!> @brief Write results to HDF5 output file.
!--------------------------------------------------------------------------------------------------
module subroutine isoductile_result(phase,group)

  integer,          intent(in) :: phase
  character(len=*), intent(in) :: group

  integer :: o


  associate(prm => param(phase), stt => damageState(phase)%state)
    outputsLoop: do o = 1,size(prm%output)
      select case(trim(prm%output(o)))
        case ('D')
          call result_writeDataset(stt,group,trim(prm%output(o)),'ductile damage driving force','-')
      end select
    end do outputsLoop
  end associate

end subroutine isoductile_result

end submodule isoductile
