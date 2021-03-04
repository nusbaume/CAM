module atm_import_export

  use shr_kind_mod  , only: r8 => shr_kind_r8, cl=>shr_kind_cl
  use time_manager  , only: get_nstep
  use cam_logfile   , only: iulog
  use spmd_utils    , only: masterproc

  implicit none

  integer     ,parameter :: debug = 0 ! internal debug level
  character(*),parameter :: F01 = "('(cam_import_export) ',a, i8,2x,i8,2x,d21.14)"

contains

  subroutine atm_import( x2a, cam_in, cam_out, restart_init )

    !-----------------------------------------------------------------------
    use cam_cpl_indices
    use camsrfexch        , only: cam_in_t, cam_out_t
    use phys_grid         , only: get_ncols_p, get_rlat_p, get_rlon_p
    use ppgrid            , only: begchunk, endchunk
    use shr_const_mod     , only: shr_const_stebol
    use shr_sys_mod       , only: shr_sys_abort
    use seq_drydep_mod    , only: n_drydep
    use shr_fire_emis_mod , only: shr_fire_emis_mechcomps_n
    use co2_cycle         , only: c_i, co2_readFlux_ocn, co2_readFlux_fuel
    use co2_cycle         , only: co2_transport, co2_time_interp_ocn, co2_time_interp_fuel
    use co2_cycle         , only: data_flux_ocn, data_flux_fuel
    use physconst         , only: mwco2
    use time_manager      , only: is_first_step

    !Water isotopes:
    use water_tracer_vars, only: wtrc_nsrfvap, wtrc_iasrfvap, wtrc_indices, wtrc_species, &
                                 trace_water
    use water_tracers    , only: wtrc_ratio
    use water_isotopes   , only: isph2o, isph216o, isphdo, isph218o

    !
    ! Arguments
    !
    real(r8)      , intent(in)    :: x2a(:,:)
    type(cam_in_t), intent(inout) :: cam_in(begchunk:endchunk)
    type(cam_out_t), intent(in)   :: cam_out(begchunk:endchunk)
    logical, optional, intent(in) :: restart_init
    !
    ! Local variables
    !
    integer            :: i,j,lat,n,c,ig  ! indices
    integer            :: ncols         ! number of columns
    logical, save      :: first_time = .true.
    integer, parameter :: ndst = 2
    integer, target    :: spc_ndx(ndst)
    integer, pointer   :: dst_a5_ndx, dst_a7_ndx
    integer, pointer   :: dst_a1_ndx, dst_a3_ndx
    integer            :: nstep
    logical            :: overwrite_flds

    !water isotopes:
    real(r8)           :: R  !water tracer ratio

    real(r8)           :: wtlat
    real(r8)           :: wtlon
    real(r8), parameter:: radtodeg = 180.0_r8/SHR_CONST_PI

    !-----------------------------------------------------------------------

    overwrite_flds = .true.
    ! don't overwrite fields if invoked during the initialization phase
    ! of a 'continue' or 'branch' run type with data from .rs file
    if (present(restart_init)) overwrite_flds = .not. restart_init

    ! ccsm sign convention is that fluxes are positive downward

    ig=1
    do c=begchunk,endchunk
       ncols = get_ncols_p(c)

       do i =1,ncols
          if (overwrite_flds) then
             cam_in(c)%wsx(i)    = -x2a(index_x2a_Faxx_taux,ig)
             cam_in(c)%wsy(i)    = -x2a(index_x2a_Faxx_tauy,ig)
             cam_in(c)%shf(i)    = -x2a(index_x2a_Faxx_sen, ig)
             cam_in(c)%cflx(i,1) = -x2a(index_x2a_Faxx_evap,ig)
          endif
          cam_in(c)%lhf(i)       = -x2a(index_x2a_Faxx_lat, ig)
          cam_in(c)%lwup(i)      = -x2a(index_x2a_Faxx_lwup,ig)

          ! Iterate over the isotopes that need to go to the surface, and try
          ! to match the ones specified with the surface field names.
          !
          ! NOTE: isph2o is total water, so is the same as Q
          !

          !Need to set this before doing water tracers:
          cam_in(c)%landfrac(i)  = x2a(index_x2a_Sf_lfrac, ig)

          !Need to define lat/lon for water tracers:
          wtlat = get_rlat_p(c,i)*radtodeg
          wtlon = get_rlon_p(c,i)*radtodeg

          if (trace_water) then
             do j = 1, wtrc_nsrfvap
               select case(wtrc_species(wtrc_iasrfvap(j)))
                 case (isph2o)
                   if(j .eq. 1) then !H2O tracer?
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                   else !water tag?
                     !H2O water tags:
                     !----------------
                     if(-x2a(index_x2a_Faxx_evap,ig) .lt. 0._r8) then !dew/frost?
                       !calculate surface vapor ratio:
                       R = wtrc_ratio(wtrc_species(wtrc_iasrfvap(j)),cam_out(c)%qbot(i,wtrc_indices(wtrc_iasrfvap(j))),&
                                     cam_out(c)%qbot(i,wtrc_indices(wtrc_iasrfvap(1))))
                       cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = R*-x2a(index_x2a_Faxx_evap,ig)
                     else
                       if(j .eq. 5) then
                         !All Land:
                         cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) =  cam_in(c)%landfrac(i)*-x2a(index_x2a_Faxx_evap,ig)
                       else if(j .eq. 6) then
                         !North Pacific:
                         if((wtlat > 30._r8) .and. (wtlat < 60._r8)) then
                           if((wtlon > 30._r8) .and. (wtlon <= 250._r8)) then
                             !NOTE:  Could also use ocnfrac here if you wanted no sea-ice fluxes.
                             cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = (1._r8-cam_in(c)%landfrac(i))*-x2a(index_x2a_Faxx_evap,ig)
                           else
                             cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                           end if
                         else
                           cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                         end if
                         !NOTE:  If you have more non-isotopic water tags, you would put them here!
                       end if !water tracers
                     end if   !dew/frost
                   !----------------
                   end if !H2O tracer
                 case (isph216o)
                   !cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap_16O,ig)
                   !Given that there are currently no isotopic land-surface fluxes, make sure to set H216O to bulk water:
                   cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                 case (isphdo)
                   cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap_HDO,ig)
                 case (isph218o)
                   cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap_18O,ig)
               end select
             end do
          end if

          cam_in(c)%asdir(i)     =  x2a(index_x2a_Sx_avsdr, ig)
          cam_in(c)%aldir(i)     =  x2a(index_x2a_Sx_anidr, ig)
          cam_in(c)%asdif(i)     =  x2a(index_x2a_Sx_avsdf, ig)
          cam_in(c)%aldif(i)     =  x2a(index_x2a_Sx_anidf, ig)
          cam_in(c)%ts(i)        =  x2a(index_x2a_Sx_t,     ig)
          cam_in(c)%sst(i)       =  x2a(index_x2a_So_t,     ig)
          cam_in(c)%snowhland(i) =  x2a(index_x2a_Sl_snowh, ig)
          cam_in(c)%snowhice(i)  =  x2a(index_x2a_Si_snowh, ig)
          cam_in(c)%tref(i)      =  x2a(index_x2a_Sx_tref,  ig)
          cam_in(c)%qref(i)      =  x2a(index_x2a_Sx_qref,  ig)
          cam_in(c)%u10(i)       =  x2a(index_x2a_Sx_u10,   ig)
          cam_in(c)%icefrac(i)   =  x2a(index_x2a_Sf_ifrac, ig)
          cam_in(c)%ocnfrac(i)   =  x2a(index_x2a_Sf_ofrac, ig)
!	  cam_in(c)%landfrac(i)  =  x2a(index_x2a_Sf_lfrac, ig)

          if ( associated(cam_in(c)%ram1) ) &
               cam_in(c)%ram1(i) =  x2a(index_x2a_Sl_ram1 , ig)
          if ( associated(cam_in(c)%fv) ) &
               cam_in(c)%fv(i)   =  x2a(index_x2a_Sl_fv   , ig)
          if ( associated(cam_in(c)%soilw) ) &
               cam_in(c)%soilw(i) =  x2a(index_x2a_Sl_soilw, ig)
          if ( associated(cam_in(c)%dstflx) ) then
             cam_in(c)%dstflx(i,1) = x2a(index_x2a_Fall_flxdst1, ig)
             cam_in(c)%dstflx(i,2) = x2a(index_x2a_Fall_flxdst2, ig)
             cam_in(c)%dstflx(i,3) = x2a(index_x2a_Fall_flxdst3, ig)
             cam_in(c)%dstflx(i,4) = x2a(index_x2a_Fall_flxdst4, ig)
          endif
          if ( associated(cam_in(c)%meganflx) ) then
             cam_in(c)%meganflx(i,1:shr_megan_mechcomps_n) = &
                  x2a(index_x2a_Fall_flxvoc:index_x2a_Fall_flxvoc+shr_megan_mechcomps_n-1, ig)
          endif

          ! Fire emission fluxes
          if ( associated(cam_in(c)%fireflx) .and. associated(cam_in(c)%fireztop) ) then
             cam_in(c)%fireflx(i,:shr_fire_emis_mechcomps_n) = &
                  x2a(index_x2a_Fall_flxfire:index_x2a_Fall_flxfire+shr_fire_emis_mechcomps_n-1, ig)
             cam_in(c)%fireztop(i) = x2a(index_x2a_Sl_ztopfire, ig)
          endif

          ! dry dep velocities
          if ( index_x2a_Sl_ddvel/=0 .and. n_drydep>0 ) then
             cam_in(c)%depvel(i,:n_drydep) = &
                  x2a(index_x2a_Sl_ddvel:index_x2a_Sl_ddvel+n_drydep-1, ig)
          endif
          !
          ! fields needed to calculate water isotopes to ocean evaporation processes
          !
          cam_in(c)%ustar(i) = x2a(index_x2a_So_ustar,ig)
          cam_in(c)%re(i)    = x2a(index_x2a_So_re   ,ig)
          cam_in(c)%ssq(i)   = x2a(index_x2a_So_ssq  ,ig)
          !
          ! bgc scenarios
          !
          if (index_x2a_Fall_fco2_lnd /= 0) then
             cam_in(c)%fco2_lnd(i) = -x2a(index_x2a_Fall_fco2_lnd,ig)
          end if
          if (index_x2a_Faoo_fco2_ocn /= 0) then
             cam_in(c)%fco2_ocn(i) = -x2a(index_x2a_Faoo_fco2_ocn,ig)
          end if
          if (index_x2a_Faoo_fdms_ocn /= 0) then
             cam_in(c)%fdms(i)     = -x2a(index_x2a_Faoo_fdms_ocn,ig)
          end if

          ig=ig+1

       end do
    end do

    ! Get total co2 flux from components,
    ! Note - co2_transport determines if cam_in(c)%cflx(i,c_i(1:4)) is allocated

    if (co2_transport().and.overwrite_flds) then

       ! Interpolate in time for flux data read in
       if (co2_readFlux_ocn) then
          call co2_time_interp_ocn
       end if
       if (co2_readFlux_fuel) then
          call co2_time_interp_fuel
       end if

       ! from ocn : data read in or from coupler or zero
       ! from fuel: data read in or zero
       ! from lnd : through coupler or zero
       do c=begchunk,endchunk
          ncols = get_ncols_p(c)
          do i=1,ncols

             ! all co2 fluxes in unit kgCO2/m2/s ! co2 flux from ocn
             if (index_x2a_Faoo_fco2_ocn /= 0) then
                cam_in(c)%cflx(i,c_i(1)) = cam_in(c)%fco2_ocn(i)
             else if (co2_readFlux_ocn) then
                ! convert from molesCO2/m2/s to kgCO2/m2/s
                cam_in(c)%cflx(i,c_i(1)) = &
                     -data_flux_ocn%co2flx(i,c)*(1._r8- cam_in(c)%landfrac(i)) &
                     *mwco2*1.0e-3_r8
             else
                cam_in(c)%cflx(i,c_i(1)) = 0._r8
             end if

             ! co2 flux from fossil fuel
             if (co2_readFlux_fuel) then
                cam_in(c)%cflx(i,c_i(2)) = data_flux_fuel%co2flx(i,c)
             else
                cam_in(c)%cflx(i,c_i(2)) = 0._r8
             end if

             ! co2 flux from land (cpl already multiplies flux by land fraction)
             if (index_x2a_Fall_fco2_lnd /= 0) then
                cam_in(c)%cflx(i,c_i(3)) = cam_in(c)%fco2_lnd(i)
             else
                cam_in(c)%cflx(i,c_i(3)) = 0._r8
             end if

             ! merged co2 flux
             cam_in(c)%cflx(i,c_i(4)) = cam_in(c)%cflx(i,c_i(1)) + &
                                        cam_in(c)%cflx(i,c_i(2)) + &
                                        cam_in(c)%cflx(i,c_i(3))
          end do
       end do
    end if
    !
    ! if first step, determine longwave up flux from the surface temperature
    !
    if (first_time) then
       if (is_first_step()) then
          do c=begchunk, endchunk
             ncols = get_ncols_p(c)
             do i=1,ncols
                cam_in(c)%lwup(i) = shr_const_stebol*(cam_in(c)%ts(i)**4)
             end do
          end do
       end if
       first_time = .false.
    end if

    !-----------------------------------------------------------------
    ! Debug output
    !-----------------------------------------------------------------

    if (debug > 0 .and. masterproc) then
       nstep = get_nstep()
       ig=1
       do c=begchunk, endchunk
          ncols = get_ncols_p(c)
          do i=1,ncols
             write(iulog,F01)'import: nstep, ig, Faxx_tauy = ',nstep,ig,x2a(index_x2a_Faxx_tauy ,ig)
             write(iulog,F01)'import: nstep, ig, Faxx_taux = ',nstep,ig,x2a(index_x2a_Faxx_taux ,ig)
             write(iulog,F01)'import: nstep, ig, Faxx_shf  = ',nstep,ig,x2a(index_x2a_Faxx_sen  ,ig)
             write(iulog,F01)'import: nstep, ig, Faxx_lhf  = ',nstep,ig,x2a(index_x2a_Faxx_lat  ,ig)
             write(iulog,F01)'import: nstep, ig, Sx_asdir  = ',nstep,ig,x2a(index_x2a_Sx_avsdr  ,ig)
             write(iulog,F01)'import: nstep, ig, Sx_aldir  = ',nstep,ig,x2a(index_x2a_Sx_anidr  ,ig)
             write(iulog,F01)'import: nstep, ig, Sx_asdif  = ',nstep,ig,x2a(index_x2a_Sx_avsdf  ,ig)
             write(iulog,F01)'import: nstep, ig, Sx_aldif  = ',nstep,ig,x2a(index_x2a_Sx_anidf  ,ig)
             write(iulog,F01)'import: nstep, ig, Sx_t      = ',nstep,ig,x2a(index_x2a_Sx_t      ,ig)
             write(iulog,F01)'import: nstep, ig, Sl_snowh  = ',nstep,ig,x2a(index_x2a_Sl_snowh  ,ig)
             write(iulog,F01)'import: nstep, ig, Si_snowh  = ',nstep,ig,x2a(index_x2a_Si_snowh  ,ig)
             write(iulog,F01)'import: nstep, ig, Sf_ifrac  = ',nstep,ig,x2a(index_x2a_Sf_ifrac  ,ig)
             write(iulog,F01)'import: nstep, ig, Sf_ofrac  = ',nstep,ig,x2a(index_x2a_Sf_ofrac  ,ig)
             write(iulog,F01)'import: nstep, ig, Sf_lfrac  = ',nstep,ig,x2a(index_x2a_Sf_lfrac  ,ig)
             if (.not. first_time .and. .not. is_first_step()) then
                write(iulog,F01)'import: nstep, ig, Faxa_lwup = ',nstep,ig,x2a(index_x2a_Faxx_lwup, ig)
             else
                write(iulog,F01)'import: nstep, ig, Faxa_lwup = ',nstep,ig,cam_in(c)%lwup(i)
             end if
             ig = ig + 1
          end do
       end do
    end if

  end subroutine atm_import

  !===============================================================================

  subroutine atm_export( cam_out, a2x )

    !-------------------------------------------------------------------
    use camsrfexch, only: cam_out_t
    use phys_grid , only: get_ncols_p
    use ppgrid    , only: begchunk, endchunk
    use cam_cpl_indices

    !Water isotopes:
    use water_tracer_vars, only: wtrc_nsrfvap, wtrc_iasrfvap, wtrc_indices, wtrc_species, &
                                 trace_water
    use water_tracers    , only: wtrc_ratio
    use water_isotopes   , only: isph2o, isph216o, isphdo, isph218o

    !
    ! Arguments
    !
    type(cam_out_t), intent(in)    :: cam_out(begchunk:endchunk)
    real(r8)       , intent(inout) :: a2x(:,:)
    !
    ! Local variables
    !
    integer :: avsize, avnat
    integer :: i,j,m,c,n,ig       ! indices
    integer :: ncols            ! Number of columns
    integer :: nstep
   !water tracers:
    logical :: pass16, passD, pass18 !logicals that prevent the passing of water tag infromation to iCLM4.
    !-----------------------------------------------------------------------

    ! Copy from component arrays into chunk array data structure
    ! Rearrange data from chunk structure into lat-lon buffer and subsequently
    ! create attribute vector

    ig=1
    do c=begchunk, endchunk
       ncols = get_ncols_p(c)
       do i=1,ncols
          a2x(index_a2x_Sa_pslv   ,ig) = cam_out(c)%psl(i)
          a2x(index_a2x_Sa_z      ,ig) = cam_out(c)%zbot(i)
          a2x(index_a2x_Sa_topo   ,ig) = cam_out(c)%topo(i)
          a2x(index_a2x_Sa_u      ,ig) = cam_out(c)%ubot(i)
          a2x(index_a2x_Sa_v      ,ig) = cam_out(c)%vbot(i)
          a2x(index_a2x_Sa_tbot   ,ig) = cam_out(c)%tbot(i)
          a2x(index_a2x_Sa_ptem   ,ig) = cam_out(c)%thbot(i)
          a2x(index_a2x_Sa_pbot   ,ig) = cam_out(c)%pbot(i)
          a2x(index_a2x_Sa_shum   ,ig) = cam_out(c)%qbot(i,1)
          !water tracers/isotopes:
          !----------------------
          !
          ! Iterate over the isotopes that need to go to the surface, and try
          ! to match the ones specified with the surface field names.
          !
          ! NOTE: isph2o is total water, so is the same as Q
          if(trace_water) then
            a2x(index_a2x_Sa_shum_16O   ,ig) = 0._r8
            a2x(index_a2x_Sa_shum_HDO   ,ig) = 0._r8
            a2x(index_a2x_Sa_shum_18O   ,ig) = 0._r8

           !logical to prevent surface vapor from tags being passed on. -JN
            pass16 = .true.
            passD  = .true.
            pass18 = .true.

            do j = 1, wtrc_nsrfvap
              select case(wtrc_species(wtrc_iasrfvap(j)))
                case (isph216o)
                  if(pass16) then !pass on H216O?
                    a2x(index_a2x_Sa_shum_16O   ,ig) = cam_out(c)%qbot(i,wtrc_indices(wtrc_iasrfvap(j)))
                    pass16 = .false.
                  end if
                case (isphdo)
                  if(passD) then !pass on HDO?
                    a2x(index_a2x_Sa_shum_HDO   ,ig) = cam_out(c)%qbot(i,wtrc_indices(wtrc_iasrfvap(j)))
                    passD = .false.
                  end if
                case (isph218o)
                  if(pass18) then !pass on H218O?
                    a2x(index_a2x_Sa_shum_18O   ,ig) = cam_out(c)%qbot(i,wtrc_indices(wtrc_iasrfvap(j)))
                    pass18 = .false.
                  end if
              end select
            end do
          end if
          !----------------------
	  a2x(index_a2x_Sa_dens   ,ig) = cam_out(c)%rho(i)
          a2x(index_a2x_Faxa_swnet,ig) = cam_out(c)%netsw(i)
          a2x(index_a2x_Faxa_lwdn ,ig) = cam_out(c)%flwds(i)
          a2x(index_a2x_Faxa_rainc,ig) = (cam_out(c)%precc(i)-cam_out(c)%precsc(i))*1000._r8
          a2x(index_a2x_Faxa_rainl,ig) = (cam_out(c)%precl(i)-cam_out(c)%precsl(i))*1000._r8
          a2x(index_a2x_Faxa_snowc,ig) = cam_out(c)%precsc(i)*1000._r8
          a2x(index_a2x_Faxa_snowl,ig) = cam_out(c)%precsl(i)*1000._r8
          a2x(index_a2x_Faxa_swndr,ig) = cam_out(c)%soll(i)
          a2x(index_a2x_Faxa_swvdr,ig) = cam_out(c)%sols(i)
          a2x(index_a2x_Faxa_swndf,ig) = cam_out(c)%solld(i)
          a2x(index_a2x_Faxa_swvdf,ig) = cam_out(c)%solsd(i)
          !water tracers/isotopes:
          !----------------------
          if(trace_water) then
          !NOTE:  converting m/s to kg/m2/s here too(may need to convert snow to equiv. water???):
            a2x(index_a2x_Faxa_rainl_16O,ig)=cam_out(c)%precrl_16O(i)*1000._r8
            a2x(index_a2x_Faxa_snowl_16O,ig)=cam_out(c)%precsl_16O(i)*1000._r8
            a2x(index_a2x_Faxa_rainc_16O,ig)=cam_out(c)%precrc_16O(i)*1000._r8
            a2x(index_a2x_Faxa_snowc_16O,ig)=cam_out(c)%precsc_16O(i)*1000._r8
            a2x(index_a2x_Faxa_rainl_HDO,ig)=cam_out(c)%precrl_HDO(i)*1000._r8
            a2x(index_a2x_Faxa_snowl_HDO,ig)=cam_out(c)%precsl_HDO(i)*1000._r8
            a2x(index_a2x_Faxa_rainc_HDO,ig)=cam_out(c)%precrc_HDO(i)*1000._r8
            a2x(index_a2x_Faxa_snowc_HDO,ig)=cam_out(c)%precsc_HDO(i)*1000._r8
            a2x(index_a2x_Faxa_rainl_18O,ig)=cam_out(c)%precrl_18O(i)*1000._r8
            a2x(index_a2x_Faxa_snowl_18O,ig)=cam_out(c)%precsl_18O(i)*1000._r8
            a2x(index_a2x_Faxa_rainc_18O,ig)=cam_out(c)%precrc_18O(i)*1000._r8
            a2x(index_a2x_Faxa_snowc_18O,ig)=cam_out(c)%precsc_18O(i)*1000._r8
          end if
          !----------------------

          ! aerosol deposition fluxes
          a2x(index_a2x_Faxa_bcphidry,ig) = cam_out(c)%bcphidry(i)
          a2x(index_a2x_Faxa_bcphodry,ig) = cam_out(c)%bcphodry(i)
          a2x(index_a2x_Faxa_bcphiwet,ig) = cam_out(c)%bcphiwet(i)
          a2x(index_a2x_Faxa_ocphidry,ig) = cam_out(c)%ocphidry(i)
          a2x(index_a2x_Faxa_ocphodry,ig) = cam_out(c)%ocphodry(i)
          a2x(index_a2x_Faxa_ocphiwet,ig) = cam_out(c)%ocphiwet(i)
          a2x(index_a2x_Faxa_dstwet1,ig)  = cam_out(c)%dstwet1(i)
          a2x(index_a2x_Faxa_dstdry1,ig)  = cam_out(c)%dstdry1(i)
          a2x(index_a2x_Faxa_dstwet2,ig)  = cam_out(c)%dstwet2(i)
          a2x(index_a2x_Faxa_dstdry2,ig)  = cam_out(c)%dstdry2(i)
          a2x(index_a2x_Faxa_dstwet3,ig)  = cam_out(c)%dstwet3(i)
          a2x(index_a2x_Faxa_dstdry3,ig)  = cam_out(c)%dstdry3(i)
          a2x(index_a2x_Faxa_dstwet4,ig)  = cam_out(c)%dstwet4(i)
          a2x(index_a2x_Faxa_dstdry4,ig)  = cam_out(c)%dstdry4(i)

          if (index_a2x_Sa_co2prog /= 0) then
             a2x(index_a2x_Sa_co2prog,ig) = cam_out(c)%co2prog(i) ! atm prognostic co2
          end if
          if (index_a2x_Sa_co2diag /= 0) then
             a2x(index_a2x_Sa_co2diag,ig) = cam_out(c)%co2diag(i) ! atm diagnostic co2
          end if
          if (index_a2x_Faxa_nhx > 0 ) then
             a2x(index_a2x_Faxa_nhx,ig) = cam_out(c)%nhx_nitrogen_flx(i)
          endif
          if (index_a2x_Faxa_noy > 0 ) then
             a2x(index_a2x_Faxa_noy,ig) = cam_out(c)%noy_nitrogen_flx(i)
          endif

          ig=ig+1
       end do
    end do

    !-----------------------------------------------------------------
    ! Debug output
    !-----------------------------------------------------------------

    if (debug > 0 .and. masterproc) then
       nstep = get_nstep()
       ig=1
       do c=begchunk, endchunk
          ncols = get_ncols_p(c)
          do i=1,ncols
             write(iulog,F01)'export: nstep, ig, Sa_z          = ',nstep,ig,a2x(index_a2x_Sa_z,ig)
             write(iulog,F01)'export: nstep, ig, Sa_topo       = ',nstep,ig,a2x(index_a2x_Sa_topo,ig)
             write(iulog,F01)'export: nstep, ig, Sa_u          = ',nstep,ig,a2x(index_a2x_Sa_u,ig)
             write(iulog,F01)'export: nstep, ig, Sa_v          = ',nstep,ig,a2x(index_a2x_Sa_v,ig)
             write(iulog,F01)'export: nstep, ig, Sa_tbot       = ',nstep,ig,a2x(index_a2x_Sa_tbot,ig)
             write(iulog,F01)'export: nstep, ig, Sa_ptem       = ',nstep,ig,a2x(index_a2x_Sa_ptem,ig)
             write(iulog,F01)'export: nstep, ig, Sa_pbot       = ',nstep,ig,a2x(index_a2x_Sa_pbot,ig)
             write(iulog,F01)'export: nstep, ig, Sa_shum       = ',nstep,ig,a2x(index_a2x_Sa_shum,ig)
             write(iulog,F01)'export: nstep, ig, Sa_dens       = ',nstep,ig,a2x(index_a2x_Sa_dens,ig)
             write(iulog,F01)'export: nstep, ig, Faxa_swnet    = ',nstep,ig,a2x(index_a2x_Faxa_swnet,ig)
             write(iulog,F01)'export: nstep, ig, Faxa_lwdn     = ',nstep,ig,a2x(index_a2x_Faxa_lwdn,ig)
             write(iulog,F01)'export: nstep, ig, Faxa_rainc    = ',nstep,ig,a2x(index_a2x_Faxa_rainc,ig)
             write(iulog,F01)'export: nstep, ig, Faxa_rainl    = ',nstep,ig,a2x(index_a2x_Faxa_rainl,ig)
             write(iulog,F01)'export: nstep, ig, Faxa_snowc    = ',nstep,ig,a2x(index_a2x_Faxa_snowc,ig)
             write(iulog,F01)'export: nstep, ig, Faxa_snowl    = ',nstep,ig,a2x(index_a2x_Faxa_snowl,ig)
             write(iulog,F01)'export: nstep, ig, Faxa_swndr    = ',nstep,ig,a2x(index_a2x_Faxa_swndr,ig)
             write(iulog,F01)'export: nstep, ig, Faxa_swvdr    = ',nstep,ig,a2x(index_a2x_Faxa_swvdr,ig)
             write(iulog,F01)'export: nstep, ig, Faxa_swndf    = ',nstep,ig,a2x(index_a2x_Faxa_swndf,ig)
             write(iulog,F01)'export: nstep, ig, Faxa_swvdf    = ',nstep,ig,a2x(index_a2x_Faxa_swvdf,ig)
             write(iulog,F01)'export: nstep, ig, Faxa_bcphidry = ',nstep,ig,a2x(index_a2x_Faxa_bcphidry,ig)
             write(iulog,F01)'export: nstep, ig, Faxa_bcphodry = ',nstep,ig,a2x(index_a2x_Faxa_bcphodry,ig)
             write(iulog,F01)'export: nstep, ig, Faxa_bcphiwet = ',nstep,ig,a2x(index_a2x_Faxa_bcphiwet,ig)
             write(iulog,F01)'export: nstep, ig, Faxa_ocphidry = ',nstep,ig,a2x(index_a2x_Faxa_ocphidry,ig)
             write(iulog,F01)'export: nstep, ig, Faxa_ocphodry = ',nstep,ig,a2x(index_a2x_Faxa_ocphodry,ig)
             write(iulog,F01)'export: nstep, ig, Faxa_ocphidry = ',nstep,ig,a2x(index_a2x_Faxa_ocphiwet,ig)
             write(iulog,F01)'export: nstep, ig, Faxa_dstwet1  = ',nstep,ig,a2x(index_a2x_Faxa_dstwet1,ig)
             write(iulog,F01)'export: nstep, ig, Faxa_dstwet1  = ',nstep,ig,a2x(index_a2x_Faxa_dstdry1,ig)
             write(iulog,F01)'export: nstep, ig, Faxa_dstwet1  = ',nstep,ig,a2x(index_a2x_Faxa_dstwet2,ig)
             write(iulog,F01)'export: nstep, ig, Faxa_dstwet1  = ',nstep,ig,a2x(index_a2x_Faxa_dstdry2,ig)
             write(iulog,F01)'export: nstep, ig, Faxa_dstwet1  = ',nstep,ig,a2x(index_a2x_Faxa_dstwet3,ig)
             write(iulog,F01)'export: nstep, ig, Faxa_dstwet1  = ',nstep,ig,a2x(index_a2x_Faxa_dstdry3,ig)
             write(iulog,F01)'export: nstep, ig, Faxa_dstwet1  = ',nstep,ig,a2x(index_a2x_Faxa_dstwet4,ig)
             write(iulog,F01)'export: nstep, ig, Faxa_dstwet1  = ',nstep,ig,a2x(index_a2x_Faxa_dstdry4,ig)
             if (index_a2x_Sa_co2prog /= 0) then
                write(iulog,F01)'export: nstep, ig, Sa_co2prog = ',nstep,ig,a2x(index_a2x_Sa_co2prog,ig)
             end if
             if (index_a2x_Sa_co2diag /= 0) then
                write(iulog,F01)'export: nstep, ig, Sa_co2diag  = ',nstep,ig,a2x(index_a2x_Sa_co2diag,ig)
             end if
             if (index_a2x_Faxa_nhx > 0 ) then
                write(iulog,F01)'export: nstep, ig, Faxa_nhx    = ',nstep,ig,a2x(index_a2x_Faxa_nhx,ig)
             endif
             if (index_a2x_Faxa_noy > 0 ) then
                write(iulog,F01)'export: nstep, ig, Faxa_noy    = ',nstep,ig,a2x(index_a2x_Faxa_noy,ig)
             endif
             ig = ig + 1
          end do
       end do
    end if

  end subroutine atm_export

end module atm_import_export
