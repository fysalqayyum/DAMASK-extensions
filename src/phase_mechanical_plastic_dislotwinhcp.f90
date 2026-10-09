! SPDX-License-Identifier: AGPL-3.0-or-later
!--------------------------------------------------------------------------------------------------
!> @author Faisal Qayyum, University of Tabuk
!> @brief material subroutine incorporating dislocation and deformation twinning physics for
!!        hexagonal close-packed crystals
!> @details Extends the dislotwin approach of Wong et al. (Acta Materialia 118:140-151, 2016)
!! to hP lattices. Twin nucleation feeds on dislocation reactions between basal and prismatic
!! slip systems; the nucleation pair for each twin system is derived at initialization from
!! Schmid tensor alignment instead of the hard-coded fcc pair table. The fcc-specific
!! cross-slip suppression probability (P_ncs) and the cos(60 deg) term in the stress for
!! infinite separation of partials do not apply to the ~90 deg hP partial geometry and are
!! omitted. No transformation (TRIP) or shear banding contributions.
!--------------------------------------------------------------------------------------------------
submodule(phase:plastic) dislotwinhcp

  type :: tParameters
    real(pREAL) :: &
      Q_cl       = 1.0_pREAL, &                                                                     !< activation energy for dislocation climb
      omega      = 1.0_pREAL, &                                                                     !< frequency factor for dislocation climb
      D          = 1.0_pREAL, &                                                                     !< grain size
      i_tw       = 1.0_pREAL, &                                                                     !< adjustment parameter to calculate MFP for twinning
      L_tw       = 1.0_pREAL, &                                                                     !< length of twin nuclei
      x_c        = 1.0_pREAL                                                                        !< critical distance for formation of twin nucleus
    type(tPolynomial) :: &
      Gamma_sf                                                                                      !< stacking fault energy
    real(pREAL),               allocatable, dimension(:) :: &
      b_sl, &                                                                                       !< magnitude of Burgers vector (m) for each slip system
      b_tw, &                                                                                       !< magnitude of Burgers vector (m) for each twin system
      Q_sl,&                                                                                        !< activation energy for glide (J) for each slip system
      v_0, &                                                                                        !< dislocation velocity prefactor (m/s) for each slip system
      t_tw, &                                                                                       !< twin thickness (m) for each twin system
      i_sl, &                                                                                       !< Adj. parameter for distance between 2 forest dislocations for each slip system
      p, &                                                                                          !< p-exponent in glide velocity
      q, &                                                                                          !< q-exponent in glide velocity
      r, &                                                                                          !< exponent in twin nucleation rate
      tau_0, &                                                                                      !< strength due to elements in solid solution
      gamma_char_tw, &                                                                              !< characteristic shear for twins
      B, &                                                                                          !< drag coefficient
      d_caron                                                                                       !< distance of spontaneous annhihilation
    real(pREAL),               allocatable, dimension(:,:) :: &
      h_sl_sl, &                                                                                    !< components of slip-slip interaction matrix
      h_sl_tw, &                                                                                    !< components of slip-twin interaction matrix
      h_tw_tw, &                                                                                    !< components of twin-twin interaction matrix
      n0_sl, &                                                                                      !< slip system normal
      forestProjection
    real(pREAL),               allocatable, dimension(:,:,:) :: &
      P_sl, &
      P_tw
    integer :: &
      sum_N_sl, &                                                                                   !< total number of active slip systems
      sum_N_tw                                                                                      !< total number of active twin systems
    integer,                   allocatable, dimension(:) :: &
      N_tw
    integer,                   allocatable, dimension(:,:) :: &
      hcp_twinNucleationSlipPair                                                                    !< (basal,prismatic) slip pair feeding each twin nucleus
    character(len=:),          allocatable :: &
      isotropic_bound
    character(len=pSTRLEN),    allocatable, dimension(:) :: &
      output
    logical :: &
      extendedDislocations = .false., &                                                             !< consider split into partials for climb calculation
      omitDipoles = .false.                                                                         !< flag controlling consideration of dipole formation
    character(len=:),          allocatable, dimension(:) :: &
      systems_sl, &
      systems_tw
  end type tParameters                                                                              !< container type for internal constitutive parameters

  type :: tIndexDotState
    integer, dimension(2) :: &
      rho_mob, &
      rho_dip, &
      gamma_sl, &
      f_tw
  end type tIndexDotState

  type :: tDislotwinhcpState
    real(pREAL),                  dimension(:,:),   pointer :: &
      rho_mob, &
      rho_dip, &
      gamma_sl, &
      f_tw
  end type tDislotwinhcpState

  type :: tDislotwinhcpDependentState
    real(pREAL),                  dimension(:,:),   allocatable :: &
      Lambda_sl, &                                                                                  !< mean free path between 2 obstacles seen by a moving dislocation
      Lambda_tw, &                                                                                  !< mean free path between 2 obstacles seen by a growing twin
      tau_pass                                                                                      !< threshold stress for slip
  end type tDislotwinhcpDependentState

!--------------------------------------------------------------------------------------------------
! containers for parameters and state
  type(tParameters),                 allocatable, dimension(:) :: param
  type(tIndexDotState),              allocatable, dimension(:) :: indexDotState
  type(tDislotwinhcpState),          allocatable, dimension(:) :: state
  type(tDislotwinhcpDependentState), allocatable, dimension(:) :: dependentState

contains


!--------------------------------------------------------------------------------------------------
!> @brief Perform module initialization.
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
module function plastic_dislotwinhcp_init() result(myPlasticity)

  logical, dimension(:), allocatable :: myPlasticity
  integer :: &
    ph, i, s, &
    Nmembers, &
    sizeState, sizeDotState, &
    startIndex, endIndex
  integer,     dimension(:), allocatable :: &
    N_sl
  real(pREAL), allocatable, dimension(:) :: &
    f_edge, &                                                                                       !< edge character fraction of total dislocation density
    rho_mob_0, &                                                                                    !< initial unipolar dislocation density per slip system
    rho_dip_0, &                                                                                    !< initial dipole dislocation density per slip system
    overlap                                                                                         !< Schmid tensor alignment between slip and twin systems
  character(len=:), allocatable :: &
    refs, &
    extmsg
  type(tDict), pointer :: &
    phases, &
    phase, &
    mech, &
    pl


  myPlasticity = plastic_active('dislotwinhcp')
  if (count(myPlasticity) == 0) return

  print'(/,1x,a)', '<<<+-  phase:mechanical:plastic:dislotwinhcp init  -+>>>'

  print'(/,1x,a)', 'A. Ma and F. Roters, Acta Materialia 52(12):3603–3612, 2004'
  print'(  1x,a)', 'https://doi.org/10.1016/j.actamat.2004.04.012'

  print'(/,1x,a)', 'S.L. Wong et al., Acta Materialia 118:140–151, 2016'
  print'(  1x,a)', 'https://doi.org/10.1016/j.actamat.2016.07.032'

  print'(/,1x,a)', 'extended to hP lattices with twin nucleation from basal-prismatic slip pairs'

  print'(/,1x,a,1x,i0)', '# phases:',count(myPlasticity); flush(IO_STDOUT)

  phases => config_material%get_dict('phase')
  allocate(param(size(phases)))
  allocate(indexDotState(size(phases)))
  allocate(state(size(phases)))
  allocate(dependentState(size(phases)))
  extmsg = ''

  do ph = 1, size(phases)
    if (.not. myPlasticity(ph)) cycle

    associate(prm => param(ph), &
              stt => state(ph), dst => dependentState(ph), &
              idx_dot => indexDotState(ph))

    phase => phases%get_dict(ph)
    mech  => phase%get_dict('mechanical')
    pl    => mech%get_dict('plastic')

    print'(/,1x,a,1x,i0,a)', 'phase',ph,': '//phases%key(ph)
    refs = config_listReferences(pl,indent=3)
    if (len(refs) > 0) print'(/,1x,a)', refs

    if (phase_lattice(ph) /= 'hP') extmsg = trim(extmsg)//' dislotwinhcp requires hP lattice'

#if defined (__GFORTRAN__)
    prm%output = output_as1dStr(pl)
#else
    prm%output = pl%get_as1dStr('output',defaultVal=emptyStrArray)
#endif

   prm%isotropic_bound = pl%get_asStr('isotropic_bound',defaultVal='isostrain')

!--------------------------------------------------------------------------------------------------
! slip related parameters
    N_sl         = pl%get_as1dInt('N_sl',defaultVal=emptyIntArray)
    prm%sum_N_sl = sum(abs(N_sl))
    slipActive: if (prm%sum_N_sl > 0) then
      prm%systems_sl = crystal_labels_slip(N_sl,phase_lattice(ph))
      prm%P_sl       = crystal_SchmidMatrix_slip(N_sl,phase_lattice(ph),phase_cOverA(ph))
      prm%n0_sl      = crystal_slip_normal(N_sl,phase_lattice(ph),phase_cOverA(ph))

      prm%extendedDislocations = pl%get_asBool('extend_dislocations',defaultVal=prm%extendedDislocations)
      prm%omitDipoles          = pl%get_asBool('omit_dipoles',       defaultVal=prm%omitDipoles)

      prm%Q_cl                 = pl%get_asReal('Q_cl')

      f_edge       = math_expand(pl%get_as1dReal('f_edge',    requiredSize=size(N_sl), &
                                                 defaultVal=[(0.5_pREAL,i=1,size(N_sl))]),N_sl)

#ifdef __GFORTRAN__
      rho_mob_0    = pl%get_as1dReal('rho_mob_0', requiredChunks=N_sl)
      rho_dip_0    = pl%get_as1dReal('rho_dip_0', requiredChunks=N_sl)
#else
      rho_mob_0    = math_expand(pl%get_as1dReal('rho_mob_0', requiredSize=size(N_sl)),N_sl)
      rho_dip_0    = math_expand(pl%get_as1dReal('rho_dip_0', requiredSize=size(N_sl)),N_sl)
#endif

      prm%v_0      = math_expand(pl%get_as1dReal('v_0',       requiredSize=size(N_sl)),N_sl)
      prm%b_sl     = math_expand(pl%get_as1dReal('b_sl',      requiredSize=size(N_sl)),N_sl)
      prm%Q_sl     = math_expand(pl%get_as1dReal('Q_sl',      requiredSize=size(N_sl)),N_sl)
      prm%i_sl     = math_expand(pl%get_as1dReal('i_sl',      requiredSize=size(N_sl)),N_sl)
      prm%p        = math_expand(pl%get_as1dReal('p_sl',      requiredSize=size(N_sl)),N_sl)
      prm%q        = math_expand(pl%get_as1dReal('q_sl',      requiredSize=size(N_sl)),N_sl)
      prm%tau_0    = math_expand(pl%get_as1dReal('tau_0',     requiredSize=size(N_sl)),N_sl)
      prm%B        = math_expand(pl%get_as1dReal('B',         requiredSize=size(N_sl), &
                                                 defaultVal=[(0.0_pREAL,i=1,size(N_sl))]),N_sl)
      prm%d_caron  = prm%b_sl *  pl%get_asReal('D_a')

      prm%h_sl_sl = crystal_interaction_SlipBySlip(N_sl,pl%get_as1dReal('h_sl-sl'),phase_lattice(ph))

      prm%forestProjection = spread(          f_edge,1,prm%sum_N_sl) &
                           * crystal_forestProjection_edge (N_sl,phase_lattice(ph),phase_cOverA(ph)) &
                           + spread(1.0_pREAL-f_edge,1,prm%sum_N_sl) &
                           * crystal_forestProjection_screw(N_sl,phase_lattice(ph),phase_cOverA(ph))

      ! multiplication factor according to crystal structure (nearest neighbors bcc vs fcc/hex)
      ! details: Argon & Moffat, Acta Metallurgica, Vol. 29, pg 293 to 299, 1981
      prm%omega = pl%get_asReal('omega', defaultVal=1000.0_pREAL) * 12.0_pREAL

      ! sanity checks
      if (    prm%Q_cl          <= 0.0_pREAL)          extmsg = trim(extmsg)//' Q_cl'
      if (any(rho_mob_0         <  0.0_pREAL))         extmsg = trim(extmsg)//' rho_mob_0'
      if (any(rho_dip_0         <  0.0_pREAL))         extmsg = trim(extmsg)//' rho_dip_0'
      if (any(prm%v_0           <  0.0_pREAL))         extmsg = trim(extmsg)//' v_0'
      if (any(prm%b_sl          <= 0.0_pREAL))         extmsg = trim(extmsg)//' b_sl'
      if (any(prm%Q_sl          <= 0.0_pREAL))         extmsg = trim(extmsg)//' Q_sl'
      if (any(prm%i_sl          <= 0.0_pREAL))         extmsg = trim(extmsg)//' i_sl'
      if (any(prm%B             <  0.0_pREAL))         extmsg = trim(extmsg)//' B'
      if (any(prm%d_caron       <  0.0_pREAL))         extmsg = trim(extmsg)//' d_caron(D_a,b_sl)'
      if (any(prm%p<=0.0_pREAL .or. prm%p>1.0_pREAL))  extmsg = trim(extmsg)//' p_sl'
      if (any(prm%q< 1.0_pREAL .or. prm%q>2.0_pREAL))  extmsg = trim(extmsg)//' q_sl'
    else slipActive
      rho_mob_0 = emptyRealArray
      rho_dip_0 = emptyRealArray
      allocate(prm%v_0, &
               prm%b_sl, &
               prm%Q_sl, &
               prm%i_sl, &
               prm%p, &
               prm%q, &
               prm%tau_0, &
               prm%B, &
               prm%d_caron, &
               source=emptyRealArray)
      allocate(prm%forestProjection(0,0), &
               prm%h_sl_sl(0,0))
    end if slipActive

!--------------------------------------------------------------------------------------------------
! twin related parameters
    prm%N_tw = pl%get_as1dInt('N_tw', defaultVal=emptyIntArray)
    prm%sum_N_tw = sum(abs(prm%N_tw))
    twinActive: if (prm%sum_N_tw > 0) then
      prm%systems_tw    = crystal_labels_twin(prm%N_tw,phase_lattice(ph))
      prm%P_tw          = crystal_SchmidMatrix_twin(prm%N_tw,phase_lattice(ph),phase_cOverA(ph))
      prm%gamma_char_tw = crystal_characteristicShear_Twin(prm%N_tw,phase_lattice(ph),phase_cOverA(ph))

      prm%L_tw             = pl%get_asReal('L_tw')
      prm%i_tw             = pl%get_asReal('i_tw')

      prm%b_tw = math_expand(pl%get_as1dReal('b_tw', requiredSize=size(prm%N_tw)),prm%N_tw)
      prm%t_tw = math_expand(pl%get_as1dReal('t_tw', requiredSize=size(prm%N_tw)),prm%N_tw)
      prm%r    = math_expand(pl%get_as1dReal('p_tw', requiredSize=size(prm%N_tw)),prm%N_tw)

      prm%h_tw_tw = crystal_interaction_TwinByTwin(prm%N_tw,pl%get_as1dReal('h_tw-tw'), &
                                                   phase_lattice(ph))

      ! sanity checks
      if (    prm%L_tw          < 0.0_pREAL)  extmsg = trim(extmsg)//' L_tw'
      if (    prm%i_tw          < 0.0_pREAL)  extmsg = trim(extmsg)//' i_tw'
      if (any(prm%b_tw          < 0.0_pREAL)) extmsg = trim(extmsg)//' b_tw'
      if (any(prm%t_tw          < 0.0_pREAL)) extmsg = trim(extmsg)//' t_tw'
      if (any(prm%r             < 0.0_pREAL)) extmsg = trim(extmsg)//' p_tw'
    else twinActive
      allocate(prm%gamma_char_tw, &
               prm%b_tw, &
               prm%t_tw, &
               prm%r, &
               source=emptyRealArray)
      allocate(prm%h_tw_tw(0,0))
    end if twinActive

!--------------------------------------------------------------------------------------------------
! parameters required for several mechanisms and their interactions
    if (prm%sum_N_sl + prm%sum_N_tw > 0) &
      prm%D    = pl%get_asReal('D')

    if (prm%sum_N_tw > 0) then
      prm%x_c  = pl%get_asReal('x_c')
      if (prm%x_c  < 0.0_pREAL)  extmsg = trim(extmsg)//' x_c'
    end if

    if (prm%sum_N_tw > 0 .or. prm%extendedDislocations) &
      prm%Gamma_sf = polynomial(pl,'Gamma_sf','T')

    slipAndTwinActive: if (prm%sum_N_sl * prm%sum_N_tw > 0) then
      prm%h_sl_tw = crystal_interaction_SlipByTwin(N_sl,prm%N_tw,pl%get_as1dReal('h_sl-tw'), &
                                                   phase_lattice(ph))

      ! Twin nucleation feeds on dislocation reactions between basal and prismatic slip;
      ! both families must be fully active.
      if (size(N_sl) < 2) then
        extmsg = trim(extmsg)//' N_sl: dislotwinhcp requires basal and prismatic slip families'
      elseif (N_sl(1) /= 3 .or. N_sl(2) /= 3) then
        extmsg = trim(extmsg)//' N_sl: dislotwinhcp twin nucleation requires N_sl(basal)=3 and N_sl(prismatic)=3'
      else
        ! For each twin system, select the (basal, prismatic) slip pair whose Schmid tensors
        ! align best with the twin Schmid tensor. This generalizes the static fcc pair table.
        allocate(prm%hcp_twinNucleationSlipPair(2,prm%sum_N_tw))
        allocate(overlap(prm%sum_N_sl))
        do i = 1, prm%sum_N_tw
          overlap = [(abs(math_tensordot(prm%P_sl(1:3,1:3,s),prm%P_tw(1:3,1:3,i))), &
                      s = 1, prm%sum_N_sl)]
          prm%hcp_twinNucleationSlipPair(1,i) =            maxloc(overlap(1:3),1)                   ! best-aligned basal system
          prm%hcp_twinNucleationSlipPair(2,i) = N_sl(1) +  maxloc(overlap(4:6),1)                   ! best-aligned prismatic system
        end do
        deallocate(overlap)
      end if
    elseif (prm%sum_N_tw > 0) then slipAndTwinActive
      extmsg = trim(extmsg)//' N_sl: dislotwinhcp twin nucleation requires active slip'
    end if slipAndTwinActive

!--------------------------------------------------------------------------------------------------
! allocate state arrays
    Nmembers  = count(material_ID_phase == ph)
    sizeDotState = size(['rho_mob ','rho_dip ','gamma_sl']) * prm%sum_N_sl &
                 + size(['f_tw'])                           * prm%sum_N_tw
    sizeState = sizeDotState

    call phase_allocateState(plasticState(ph),Nmembers,sizeState,sizeDotState,0)
    deallocate(plasticState(ph)%dotState)

!--------------------------------------------------------------------------------------------------
! state aliases and initialization
    startIndex = 1
    endIndex   = prm%sum_N_sl
    idx_dot%rho_mob = [startIndex,endIndex]
    stt%rho_mob => plasticState(ph)%state(startIndex:endIndex,:)
    stt%rho_mob = spread(rho_mob_0,2,Nmembers)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asReal('atol_rho',defaultVal=1.0_pREAL)
    if (any(plasticState(ph)%atol(startIndex:endIndex) < 0.0_pREAL)) extmsg = trim(extmsg)//' atol_rho'

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_sl
    idx_dot%rho_dip = [startIndex,endIndex]
    stt%rho_dip => plasticState(ph)%state(startIndex:endIndex,:)
    stt%rho_dip = spread(rho_dip_0,2,Nmembers)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asReal('atol_rho',defaultVal=1.0_pREAL)

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_sl
    idx_dot%gamma_sl = [startIndex,endIndex]
    stt%gamma_sl => plasticState(ph)%state(startIndex:endIndex,:)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asReal('atol_gamma',defaultVal=1.0e-6_pREAL)
    if (any(plasticState(ph)%atol(startIndex:endIndex) < 0.0_pREAL)) extmsg = trim(extmsg)//' atol_gamma'

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_tw
    idx_dot%f_tw = [startIndex,endIndex]
    stt%f_tw => plasticState(ph)%state(startIndex:endIndex,:)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asReal('atol_f_tw',defaultVal=1.0e-6_pREAL)
    if (any(plasticState(ph)%atol(startIndex:endIndex) < 0.0_pREAL)) extmsg = trim(extmsg)//' atol_f_tw'

    allocate(dst%tau_pass (prm%sum_N_sl,Nmembers),source=0.0_pREAL)
    allocate(dst%Lambda_sl(prm%sum_N_sl,Nmembers),source=0.0_pREAL)
    allocate(dst%Lambda_tw(prm%sum_N_tw,Nmembers),source=0.0_pREAL)

    end associate

!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
    if (extmsg /= '') call IO_error(211,ext_msg=trim(extmsg))

  end do

end function plastic_dislotwinhcp_init


!--------------------------------------------------------------------------------------------------
!> @brief Return the homogenized elasticity matrix.
!--------------------------------------------------------------------------------------------------
module function plastic_dislotwinhcp_homogenizedC(ph,en) result(homogenizedC)

  integer,     intent(in) :: &
    ph, en
  real(pREAL), dimension(6,6) :: &
    homogenizedC, &
    C
  real(pREAL), dimension(:,:,:), allocatable :: &
    C66_tw
  integer :: i
  real(pREAL) :: f_matrix


  C = elastic_C66(ph,en)

  associate(prm => param(ph), stt => state(ph))

    f_matrix = 1.0_pREAL &
             - sum(stt%f_tw(1:prm%sum_N_tw,en))

    homogenizedC = f_matrix * C

    twinActive: if (prm%sum_N_tw > 0) then
      C66_tw    = crystal_C66_twin(prm%N_tw,C,phase_lattice(ph),phase_cOverA(ph))
      do i = 1, prm%sum_N_tw
        homogenizedC = homogenizedC &
                     + stt%f_tw(i,en)*C66_tw(1:6,1:6,i)
      end do
     end if twinActive

  end associate

end function plastic_dislotwinhcp_homogenizedC


!--------------------------------------------------------------------------------------------------
!> @brief Calculate plastic velocity gradient and its tangent.
!--------------------------------------------------------------------------------------------------
module subroutine dislotwinhcp_LpAndItsTangent(Lp,dLp_dMp,Mp,ph,en)

  real(pREAL), dimension(3,3),     intent(out) :: Lp
  real(pREAL), dimension(3,3,3,3), intent(out) :: dLp_dMp
  real(pREAL), dimension(3,3),     intent(in)  :: Mp
  integer,                         intent(in)  :: ph,en

  integer :: i,k,l,m,n
  real(pREAL) :: &
    f_matrix, &
    T
  real(pREAL), dimension(param(ph)%sum_N_sl) :: &
    dot_gamma_sl,ddot_gamma_dtau_sl
  real(pREAL), dimension(param(ph)%sum_N_tw) :: &
    dot_gamma_tw,ddot_gamma_dtau_tw


  T = thermal_T(ph,en)
  Lp = 0.0_pREAL
  dLp_dMp = 0.0_pREAL

  associate(prm => param(ph), stt => state(ph))

    f_matrix = 1.0_pREAL &
             - sum(stt%f_tw(1:prm%sum_N_tw,en))

    call kinetics_sl(Mp,T,ph,en,dot_gamma_sl,ddot_gamma_dtau_sl)
    slipContribution: do i = 1, prm%sum_N_sl
      Lp = Lp + dot_gamma_sl(i)*prm%P_sl(1:3,1:3,i)
      forall (k=1:3,l=1:3,m=1:3,n=1:3) &
        dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                         + ddot_gamma_dtau_sl(i) * prm%P_sl(k,l,i) * prm%P_sl(m,n,i)
    end do slipContribution

    if (prm%sum_N_tw > 0) call kinetics_tw(Mp,T,dot_gamma_sl,ph,en,dot_gamma_tw,ddot_gamma_dtau_tw)
    twinContibution: do i = 1, prm%sum_N_tw
      Lp = Lp + dot_gamma_tw(i)*prm%P_tw(1:3,1:3,i)
      forall (k=1:3,l=1:3,m=1:3,n=1:3) &
        dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                         + ddot_gamma_dtau_tw(i)* prm%P_tw(k,l,i)*prm%P_tw(m,n,i)
    end do twinContibution

    Lp      = Lp      * f_matrix
    dLp_dMp = dLp_dMp * f_matrix

    end associate

end subroutine dislotwinhcp_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief Calculate the rate of change of microstructure.
!--------------------------------------------------------------------------------------------------
module function dislotwinhcp_dotState(Mp,ph,en) result(dotState)

  real(pREAL), dimension(3,3),  intent(in):: &
    Mp                                                                                              !< Mandel stress
  integer,                      intent(in) :: &
    ph, &
    en
  real(pREAL), dimension(plasticState(ph)%sizeDotState) :: &
    dotState

  integer :: i
  real(pREAL) :: &
    f_matrix, &
    d_hat, &
    v_cl, &                                                                                         !< climb velocity
    tau, &
    sigma_cl, &                                                                                     !< climb stress
    b_d                                                                                             !< ratio of Burgers vector to stacking fault width
  real(pREAL), dimension(param(ph)%sum_N_sl) :: &
    dot_rho_dip_formation, &
    dot_rho_dip_climb, &
    dot_gamma_sl
  real(pREAL), dimension(param(ph)%sum_N_tw) :: &
    dot_gamma_tw
  real(pREAL) :: &
    mu, nu, &
    T


  associate(prm => param(ph), stt => state(ph), dst => dependentState(ph), &
            dot_rho_mob => dotState(indexDotState(ph)%rho_mob(1):indexDotState(ph)%rho_mob(2)), &
            dot_rho_dip => dotState(indexDotState(ph)%rho_dip(1):indexDotState(ph)%rho_dip(2)), &
            abs_dot_gamma_sl => dotState(indexDotState(ph)%gamma_sl(1):indexDotState(ph)%gamma_sl(2)), &
            dot_f_tw => dotState(indexDotState(ph)%f_tw(1):indexDotState(ph)%f_tw(2)))

    mu = elastic_mu(ph,en,prm%isotropic_bound)
    nu = elastic_nu(ph,en,prm%isotropic_bound)
    T = thermal_T(ph,en)

    f_matrix = 1.0_pREAL &
             - sum(stt%f_tw(1:prm%sum_N_tw,en))

    call kinetics_sl(Mp,T,ph,en,dot_gamma_sl)
    abs_dot_gamma_sl = abs(dot_gamma_sl)

    slipState: do i = 1, prm%sum_N_sl
      tau = math_tensordot(Mp,prm%P_sl(1:3,1:3,i))

      significantSlipStress: if (dEq0(tau) .or. prm%omitDipoles) then
        d_hat = dst%Lambda_sl(i,en)
        dot_rho_dip_formation(i) = 0.0_pREAL
      else significantSlipStress
        d_hat = 3.0_pREAL*mu*prm%b_sl(i)/(16.0_pREAL*PI*abs(tau))
        d_hat = math_clip(d_hat, right = dst%Lambda_sl(i,en))
        d_hat = math_clip(d_hat, left  = prm%d_caron(i))

        dot_rho_dip_formation(i) = 2.0_pREAL*(d_hat-prm%d_caron(i))/prm%b_sl(i) &
                                 * stt%rho_mob(i,en)*abs_dot_gamma_sl(i)
      end if significantSlipStress

      if (dEq(d_hat,prm%d_caron(i))) then
        dot_rho_dip_climb(i) = 0.0_pREAL
      else
        ! Argon & Moffat, Acta Metallurgica, Vol. 29, pg 293 to 299, 1981
        sigma_cl = dot_product(prm%n0_sl(1:3,i),matmul(Mp,prm%n0_sl(1:3,i)))
        if (prm%extendedDislocations) then
          b_d = 24.0_pREAL*PI*(1.0_pREAL - nu)/(2.0_pREAL + nu) * prm%Gamma_sf%at(T) / (mu*prm%b_sl(i))
        else
          b_d = 1.0_pREAL
        end if
        v_cl = 2.0_pREAL*prm%omega*b_d**2*exp(-prm%Q_cl/(K_B*T)) &
             * (exp(abs(sigma_cl)*prm%b_sl(i)**3/(K_B*T)) - 1.0_pREAL)
        dot_rho_dip_climb(i) = 4.0_pREAL*v_cl*stt%rho_dip(i,en) &
                             / (d_hat-prm%d_caron(i))
      end if
    end do slipState

    dot_rho_mob = abs_dot_gamma_sl/(prm%b_sl*dst%Lambda_sl(:,en)) &
                - dot_rho_dip_formation &
                - 2.0_pREAL*prm%d_caron/prm%b_sl * stt%rho_mob(:,en)*abs_dot_gamma_sl

    dot_rho_dip = dot_rho_dip_formation &
                - 2.0_pREAL*prm%d_caron/prm%b_sl * stt%rho_dip(:,en)*abs_dot_gamma_sl &
                - dot_rho_dip_climb

    if (prm%sum_N_tw > 0) call kinetics_tw(Mp,T,abs_dot_gamma_sl,ph,en,dot_gamma_tw)
    dot_f_tw = f_matrix*dot_gamma_tw/prm%gamma_char_tw

  end associate

end function dislotwinhcp_dotState


!--------------------------------------------------------------------------------------------------
!> @brief Sum of the absolute shear rates over all active slip and twin systems.
!> @details Equivalent to DAMASK2's plasticState(phase)%slipRate consumed by ductile damage
!!          source models; recomputes the same kinetics as dislotwinhcp_dotState.
!--------------------------------------------------------------------------------------------------
module function dislotwinhcp_dotGammaSum(Mp,ph,en) result(dotGammaSum)

  real(pREAL), dimension(3,3),  intent(in):: &
    Mp                                                                                              !< Mandel stress
  integer,                      intent(in) :: &
    ph, &
    en
  real(pREAL) :: &
    dotGammaSum, &
    T
  real(pREAL), dimension(param(ph)%sum_N_sl) :: &
    dot_gamma_sl
  real(pREAL), dimension(param(ph)%sum_N_tw) :: &
    dot_gamma_tw


  associate(prm => param(ph))

    T = thermal_T(ph,en)

    call kinetics_sl(Mp,T,ph,en,dot_gamma_sl)
    dotGammaSum = sum(abs(dot_gamma_sl))

    if (prm%sum_N_tw > 0) then
      call kinetics_tw(Mp,T,abs(dot_gamma_sl),ph,en,dot_gamma_tw)
      dotGammaSum = dotGammaSum + sum(abs(dot_gamma_tw))
    end if

  end associate

end function dislotwinhcp_dotGammaSum


!--------------------------------------------------------------------------------------------------
!> @brief Calculate derived quantities from state.
!--------------------------------------------------------------------------------------------------
module subroutine dislotwinhcp_dependentState(ph,en)

  integer,       intent(in) :: &
    ph, &
    en

  real(pREAL) :: &
    sumf_tw
  real(pREAL), dimension(param(ph)%sum_N_sl) :: &
    inv_lambda_sl
  real(pREAL), dimension(param(ph)%sum_N_tw) :: &
    inv_lambda_tw_tw, &                                                                             !< 1/mean free distance between 2 twin stacks from different systems seen by a growing twin
    f_over_t_tw
  real(pREAL) :: &
    mu


  associate(prm => param(ph), stt => state(ph), dst => dependentState(ph))

    mu = elastic_mu(ph,en,prm%isotropic_bound)
    sumf_tw = sum(stt%f_tw(1:prm%sum_N_tw,en))

    !* rescaled volume fraction for topology
    f_over_t_tw = stt%f_tw(1:prm%sum_N_tw,en)/prm%t_tw                                              ! this is per system

    inv_lambda_sl = sqrt(matmul(prm%forestProjection,stt%rho_mob(:,en)+stt%rho_dip(:,en)))/prm%i_sl
    if (prm%sum_N_tw > 0 .and. prm%sum_N_sl > 0) &
      inv_lambda_sl = inv_lambda_sl + matmul(prm%h_sl_tw,f_over_t_tw)/(1.0_pREAL-sumf_tw)
    dst%Lambda_sl(:,en) = prm%D / (1.0_pREAL+prm%D*inv_lambda_sl)

    inv_lambda_tw_tw = matmul(prm%h_tw_tw,f_over_t_tw)/(1.0_pREAL-sumf_tw)
    dst%Lambda_tw(:,en) = prm%i_tw*prm%D/(1.0_pREAL+prm%D*inv_lambda_tw_tw)

    !* threshold stress for dislocation motion
    dst%tau_pass(:,en) = mu*prm%b_sl* sqrt(matmul(prm%h_sl_sl,stt%rho_mob(:,en)+stt%rho_dip(:,en)))

  end associate

end subroutine dislotwinhcp_dependentState


!--------------------------------------------------------------------------------------------------
!> @brief Write results to HDF5 output file.
!--------------------------------------------------------------------------------------------------
module subroutine plastic_dislotwinhcp_result(ph,group)

  integer,          intent(in) :: ph
  character(len=*), intent(in) :: group

  integer :: ou


  associate(prm => param(ph), stt => state(ph), dst => dependentState(ph))

    do ou = 1,size(prm%output)

      select case(trim(prm%output(ou)))

        case('rho_mob')
          call result_writeDataset(stt%rho_mob,group,trim(prm%output(ou)), &
                                   'mobile dislocation density','1/m²',prm%systems_sl)
        case('rho_dip')
          call result_writeDataset(stt%rho_dip,group,trim(prm%output(ou)), &
                                   'dislocation dipole density','1/m²',prm%systems_sl)
        case('gamma_sl')
          call result_writeDataset(stt%gamma_sl,group,trim(prm%output(ou)), &
                                   'plastic shear','1',prm%systems_sl)
        case('Lambda_sl')
          call result_writeDataset(dst%Lambda_sl,group,trim(prm%output(ou)), &
                                   'mean free path for slip','m',prm%systems_sl)
        case('tau_pass')
          call result_writeDataset(dst%tau_pass,group,trim(prm%output(ou)), &
                                   'passing stress for slip','Pa',prm%systems_sl)

        case('f_tw')
          call result_writeDataset(stt%f_tw,group,trim(prm%output(ou)), &
                                   'twinned volume fraction','m³/m³',prm%systems_tw)
        case('Lambda_tw')
          call result_writeDataset(dst%Lambda_tw,group,trim(prm%output(ou)), &
                                   'mean free path for twinning','m',prm%systems_tw)

      end select

    end do

  end associate

end subroutine plastic_dislotwinhcp_result


!--------------------------------------------------------------------------------------------------
!> @brief Calculate shear rates on slip systems, their derivatives with respect to resolved
!         stress, and the resolved stress.
!> @details Derivatives and resolved stress are calculated only optionally.
! NOTE: Contrary to common convention, here the result (i.e. intent(out)) variables have to be put
! at the end since some of them are optional.
!--------------------------------------------------------------------------------------------------
pure subroutine kinetics_sl(Mp,T,ph,en, &
                            dot_gamma_sl,ddot_gamma_dtau_sl,tau_sl)

  real(pREAL), dimension(3,3),  intent(in) :: &
    Mp                                                                                              !< Mandel stress
  real(pREAL),                  intent(in) :: &
    T                                                                                               !< temperature
  integer,                      intent(in) :: &
    ph, &
    en
  real(pREAL), dimension(param(ph)%sum_N_sl), intent(out) :: &
    dot_gamma_sl
  real(pREAL), dimension(param(ph)%sum_N_sl), optional, intent(out) :: &
    ddot_gamma_dtau_sl, &
    tau_sl

  real(pREAL), dimension(param(ph)%sum_N_sl) :: &
    ddot_gamma_dtau
  real(pREAL), dimension(param(ph)%sum_N_sl) :: &
    tau, &
    stressRatio, &
    StressRatio_p, &
    Q_kB_T, &
    v_wait_inverse, &                                                                               !< inverse of the effective velocity of a dislocation waiting at obstacles (unsigned)
    v_run_inverse, &                                                                                !< inverse of the velocity of a free moving dislocation (unsigned)
    dV_wait_inverse_dTau, &
    dV_run_inverse_dTau, &
    dV_dTau, &
    tau_eff                                                                                         !< effective resolved stress
  integer :: i


  associate(prm => param(ph), stt => state(ph), dst => dependentState(ph))

    tau = [(math_tensordot(Mp,prm%P_sl(1:3,1:3,i)),i = 1, prm%sum_N_sl)]

    tau_eff = abs(tau)-dst%tau_pass(:,en)

    significantStress: where(tau_eff > tol_math_check)
      ! Unlike cubic crystals, hP combines very soft (basal) with very hard (pyramidal) systems,
      ! so the resolved stress on a soft system routinely exceeds its barrier strength tau_0.
      ! Clip the ratio at 1: the barrier is fully overcome (athermal limit, v_wait = v_0) and
      ! the dislocation velocity becomes drag-controlled via v_run.
      stressRatio    = math_clip(tau_eff/prm%tau_0, 0.0_pREAL, 1.0_pREAL)
      StressRatio_p  = stressRatio** prm%p
      Q_kB_T = prm%Q_sl/(K_B*T)
      v_wait_inverse = exp(Q_kB_T*(1.0_pREAL-StressRatio_p)** prm%q) &
                     / prm%v_0
      v_run_inverse  = prm%B/(tau_eff*prm%b_sl)

      dot_gamma_sl = sign(stt%rho_mob(:,en)*prm%b_sl/(v_wait_inverse+v_run_inverse),tau)

      dV_wait_inverse_dTau = -1.0_pREAL * v_wait_inverse * prm%p * prm%q * Q_kB_T &
                           * (stressRatio**(prm%p-1.0_pREAL)) &
                           * (1.0_pREAL-StressRatio_p)**(prm%q-1.0_pREAL) &
                           / prm%tau_0
      dV_run_inverse_dTau  = -1.0_pREAL * v_run_inverse/tau_eff
      dV_dTau              = -1.0_pREAL * (dV_wait_inverse_dTau+dV_run_inverse_dTau) &
                           / (v_wait_inverse+v_run_inverse)**2
      ddot_gamma_dtau = dV_dTau*stt%rho_mob(:,en)*prm%b_sl
    else where significantStress
      dot_gamma_sl    = 0.0_pREAL
      ddot_gamma_dtau = 0.0_pREAL
    end where significantStress

  end associate

  if (present(ddot_gamma_dtau_sl)) ddot_gamma_dtau_sl = ddot_gamma_dtau
  if (present(tau_sl))             tau_sl             = tau

end subroutine kinetics_sl


!--------------------------------------------------------------------------------------------------
!> @brief Calculate shear rates on twin systems and their derivatives with respect to resolved
!         stress.
!> @details Derivatives are calculated only optionally.
!! Nuclei form at the intersection of basal and prismatic dislocations; the feeding pair for each
!! twin system was selected at initialization by Schmid tensor alignment. Compared to the fcc
!! formulation: (1) the partials bounding an hP twin nucleus are ~90 deg apart, so the
!! cos(60 deg)/x0 term in the stress for infinite separation (tau_r) vanishes; (2) there is no
!! cross-slip based nucleation suppression (P_ncs).
! NOTE: Contrary to common convention, here the result (i.e. intent(out)) variables have to be put
! at the end since some of them are optional.
!--------------------------------------------------------------------------------------------------
pure subroutine kinetics_tw(Mp,T,abs_dot_gamma_sl,ph,en,&
                            dot_gamma_tw,ddot_gamma_dtau_tw)

  real(pREAL), dimension(3,3),  intent(in) :: &
    Mp                                                                                              !< Mandel stress
  real(pREAL),                  intent(in) :: &
    T                                                                                               !< temperature
  integer,                      intent(in) :: &
    ph, &
    en
  real(pREAL), dimension(param(ph)%sum_N_sl), intent(in) :: &
    abs_dot_gamma_sl
  real(pREAL), dimension(param(ph)%sum_N_tw), intent(out) :: &
    dot_gamma_tw
  real(pREAL), dimension(param(ph)%sum_N_tw), optional, intent(out) :: &
    ddot_gamma_dtau_tw

  real(pREAL) :: &
    tau, tau_r, tau_hat, &
    dot_N_0, &
    x0, V, &
    Gamma_sf, &
    mu, nu, &
    P, dP_dtau
  integer, dimension(2) :: &
    s
  integer :: i


  associate(prm => param(ph), stt => state(ph), dst => dependentState(ph))

    mu = elastic_mu(ph,en,prm%isotropic_bound)
    nu = elastic_nu(ph,en,prm%isotropic_bound)
    Gamma_sf = prm%Gamma_sf%at(T)

    tau_hat = 3.0_pREAL*prm%b_tw(1)*mu/prm%L_tw &
            + Gamma_sf/(3.0_pREAL*prm%b_tw(1))
    x0 = mu*prm%b_sl(1)**2*(2.0_pREAL+nu)/(Gamma_sf*8.0_pREAL*PI*(1.0_pREAL-nu))
    tau_r = mu*prm%b_sl(1)/(2.0_pREAL*PI)/(x0+prm%x_c)                                              ! hP: cos(90 deg) = 0 removes the second term of the fcc expression

    do i = 1, prm%sum_N_tw
      tau = math_tensordot(Mp,prm%P_tw(1:3,1:3,i))

      if (tau > tol_math_check .and. tau < tau_r) then
        P = exp(-(tau_hat/tau)**prm%r(i))
        dP_dTau = prm%r(i) * (tau_hat/tau)**prm%r(i)/tau * P

        s = prm%hcp_twinNucleationSlipPair(1:2,i)
        dot_N_0 = sum(abs_dot_gamma_sl(s(2:1:-1))*(stt%rho_mob(s,en)+stt%rho_dip(s,en)))/(prm%L_tw*3.0_pREAL)

        V = PI/4.0_pREAL*dst%Lambda_tw(i,en)**2*prm%t_tw(i)
        dot_gamma_tw(i) = V*dot_N_0*P*prm%gamma_char_tw(i)
        if (present(ddot_gamma_dtau_tw)) &
          ddot_gamma_dtau_tw(i) = V*dot_N_0*dP_dtau*prm%gamma_char_tw(i)
      else
        dot_gamma_tw(i) = 0.0_pREAL
        if (present(ddot_gamma_dtau_tw)) ddot_gamma_dtau_tw(i) = 0.0_pREAL
      end if
    end do

  end associate

end subroutine kinetics_tw

end submodule dislotwinhcp
