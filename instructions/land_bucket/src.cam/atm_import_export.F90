module atm_import_export

  use shr_kind_mod  , only: r8 => shr_kind_r8, cl=>shr_kind_cl
  use time_manager  , only: get_nstep
  use cam_logfile   , only: iulog
  use spmd_utils    , only: masterproc

  implicit none

  integer     ,parameter :: debug = 0 ! internal debug level
  character(*),parameter :: F01 = "('(cam_import_export) ',a, i8,2x,i8,2x,d21.14)"

contains

  subroutine atm_import( x2a, cam_in, cam_out, restart_init ) ! add cam_out - RPF 

    !-----------------------------------------------------------------------
    use cam_cpl_indices
    use camsrfexch        , only: cam_in_t, cam_out_t ! to get cam_out -RPF
    use phys_grid         , only: get_ncols_p, &
                                  get_rlat_p, get_rlon_p 
    use ppgrid            , only: begchunk, endchunk
    use shr_const_mod     , only: shr_const_stebol
    use shr_sys_mod       , only: shr_sys_abort 
    use seq_drydep_mod    , only: n_drydep
    use shr_fire_emis_mod , only: shr_fire_emis_mechcomps_n
    use co2_cycle         , only: c_i, co2_readFlux_ocn, co2_readFlux_fuel
    use co2_cycle         , only: co2_transport, co2_time_interp_ocn, co2_time_interp_fuel
    use co2_cycle         , only: data_flux_ocn, data_flux_fuel
    use physconst         , only: mwco2, gravit, latvap, cpair, cpwv, rair, rh2o, cappa !RPF
    use time_manager      , only: is_first_step, is_first_restart_step, get_curr_time, get_step_size !RPF 

    !Water isotopes:
    use water_tracer_vars, only: wtrc_nsrfvap, wtrc_iasrfvap, wtrc_indices, wtrc_species, &
                                 trace_water
    use water_tracers    , only: wtrc_ratio
    use water_isotopes   , only: isph2o, isph216o, isphdo, isph218o, isph217o, isphto, &
                                  wiso_alpl, wiso_kmol_land ! RPF - alpl added for bucket, and equal fractionation

    !Additional modules for water tags:
    use wv_saturation    , only: qsat

    !
    ! Arguments
    !
    real(r8)      , intent(in)    :: x2a(:,:)
    type(cam_in_t), intent(inout) :: cam_in(begchunk:endchunk)
    type(cam_out_t), intent(inout) :: cam_out(begchunk:endchunk)
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
    ! water isotopes
    real(r8) :: R !water tracer ratio.

    !
    ! Bucket model variables!
    !    
    !adding a water isotope bucket model (added from code from JN)
    !modified on 200818 to add a second bucket. -RPF

    integer  :: dtime                   !time step (in seconds)
    real(r8) :: Re,prcsum,wtprsm        !water tracer evap ratio, precip. sum, water tracers precip. sum
    real(r8) :: dbuckmass,sbuckmass     !isotopic water mass (e.g., col depth in meters) in the deep and shallow buckets
    real(r8) :: dh2obckmass,sh2obckmass !bulk water mass (in meters) in the bucket
    real(r8) :: ovallt, ovalhv          !variables to preserve ratio of upper bucket when it runs dry (original value light/heavy)
    real(r8) :: srunoff, drunoff !mass of water (in meters) that leaves buckets as runoff
    real(r8) :: shalet, deepet   !ET from shallow and deep buckets
    real(r8) :: wtlat            !latitude in degrees
    logical  :: empty_buck       ! logical to recalculate shallow bucket if empty.
    real(r8), parameter :: dbcklim = 0.300_r8 !depth of deep bucket in m 
    real(r8), parameter :: sbcklim = 0.050_r8 !depth of shallow bucket in m 
    integer  :: gwopt = 3        !option for "groundwater" parameterization:
                                 ! 1) JN's GNIP parameterization, 2) from Bowen 2010, 3) from IsoMAP generated on 200720, 4) SMOW
    integer  :: lfopt = 1        ! allow land surface fractionation? 0 = FALSE, 1 = TRUE. applied only to soil evap. 
                                 ! transpiration and interception assumed to be non-fractionating
    real(r8) :: Rf               ! ratio associated w/ soil evaporation
    real(r8) :: nwgt             ! aerodynamic weighting term for scaling kinetic fractionation based on mass of upper bucket
    real(r8) :: aeq              ! equilibrium fractionation
    real(r8) :: akin             ! kinetic fractionation factor, calculated following Benettin et al 2018, Gat 1996, etc.
    real(r8) :: aevap            ! evaporation fractionaion factor (effective) following luz et al 2009
    real(r8) :: R17              ! storage var to hold frac factors for H217O
    real(r8), parameter :: Rmax = 100._r8 
    real(r8), parameter :: Rmin = -50.0_r8   ! hold for min/max 17O bucket isotope ratio.
    real(r8) :: Rp               ! adjust for h2o and h216o
    ! new version of getting kinetic fractionation from wiso_kmol:
    real(r8) :: ustar, rho, zsrf ! parameters that are inputs to wiso_kmol
    real(r8) :: tau, z0m         ! placeholder for total wind stress (kg m /s2)
    real(r8) :: ocnflux   ! ocean flux.
    real(r8) :: lndflux   ! total land evap flux
    real(r8) :: lakeflux  ! total lake evap flux
    real(r8), parameter :: radtodeg = 180.0_r8/SHR_CONST_PI
    real(r8) :: evaprh_ts  ! RH wrt TSRF
    real(r8) :: evaprh_es ! placeholder, saturation vapor pressure.
    real(r8) :: evaprh_qs ! return qsat value from qsat subroutine.
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

          ! RPF 200603 - move x2a fluxes up for species that are useful
          ! to design tags. note: taking reference height q/t instead of
          ! surface for now, but it might be more appropriate in the future
          ! to adjust.
          
          cam_in(c)%landfrac(i)  =  x2a(index_x2a_Sf_lfrac, ig)
          cam_in(c)%icefrac(i)   =  x2a(index_x2a_Sf_ifrac, ig)
          cam_in(c)%ocnfrac(i)   =  x2a(index_x2a_Sf_ofrac, ig)
          cam_in(c)%tref(i)      =  x2a(index_x2a_Sx_tref,  ig)
          cam_in(c)%ts(i)        =  x2a(index_x2a_Sx_t,     ig)
          cam_in(c)%qref(i)      =  x2a(index_x2a_Sx_qref,  ig)
          cam_in(c)%u10(i)       =  x2a(index_x2a_Sx_u10,   ig)
          cam_in(c)%z0m(i)       =  x2a(index_x2a_Sl_z0m,   ig)

          ! RPF calculate required quantities for the bucket model
          prcsum = cam_out(c)%precc(i) + cam_out(c)%precl(i) ! get total precip from cam_out
          dtime  = get_step_size() !determine time step size
          sh2obckmass = cam_out(c)%SbuckH(i) !<-will cause failure until camsrfexch.F90 is updated!
          dh2obckmass = cam_out(c)%DbuckH(i) !<-will cause failure until camsrfexch.F90 is updated!
          wtlat = get_rlat_p(c,i)*radtodeg

          !calculate u* here for alpha_kinetic
          tau = sqrt(cam_in(c)%wsx(i)**2 + cam_in(c)%wsy(i)**2)
          rho = cam_out(c)%rho(i)
          ustar = sqrt(tau/rho)

          ! calculate roughness length from 10m wind speed and ustar.
          z0m = cam_in(c)%z0m(i)

          ! calculate RH wrt surface for kinetic fractionation.
          if (is_first_step() .or. is_first_restart_step()) then 
            evaprh_ts   = 0._r8
          else
            ! calculate RH based using cam_in%qref, cam_in%ts, and cam_out%ps.
            call qsat(cam_in(c)%ts(i),cam_out(c)%pbot(i),evaprh_es,evaprh_qs)
            evaprh_ts = max(0._r8,min(1._r8,cam_out(c)%qbot(i,1)/evaprh_qs))
          end if
          
          ! Iterate over the isotopes that need to go to the surface, and try
          ! to match the ones specified with the surface field names.
          !
          ! NOTE: isph2o is total water, so is the same as Q
          !

          empty_buck = .false. !RPF - start each chunk/col assuming bucket isn't empty.

          if (trace_water) then
             do j = 1, wtrc_nsrfvap
               select case(wtrc_species(wtrc_iasrfvap(j)))
                 case (isph2o)
                   if (j .eq. 1) then !if working on the H2O tracer!
                     !------------------------------------------------
                     ! bucket model for evaporation!
                     !------------------------------------------------
                     !initialize bucket 'flux' terms
                     shalet = 0._r8
                     deepet = 0._r8
                     srunoff = 0._r8
                     drunoff = 0._r8
                     Rf = 0._r8

                     ! calculate ocnflux
                     lndflux = -x2a(index_x2a_Faxx_levap,ig)
                     ocnflux = -x2a(index_x2a_Faxx_evap,ig) - lndflux
                     lakeflux = lndflux - -x2a(index_x2a_Faxx_evpnl,ig)

                     ! zero out ice bucket over land, and land buckets over ice.
                     if (cam_in(c)%icefrac(i) .eq. 1._r8 .or. cam_in(c)%landfrac(i) .eq. 0._r8) then
                       cam_out(c)%SbuckH(i) = 0._r8
                       cam_out(c)%DbuckH(i) = 0._r8
                     end if

                     ! are we on land? if so, do bucket calculations.
                     if (cam_in(c)%ocnfrac(i) .eq. 1._r8) then
                       ! zero out bucket over the ocean, since this could cause issues for areas
                       ! that flip between ice and ocean.
                       cam_out(c)%SbuckH(i) = 0._r8
                       cam_out(c)%DbuckH(i) = 0._r8
                     else
                       !200824 - RPF - IMPORTANT NOTE: the bucket scheme here is *not*
                       !designed to conserve mass; mass conservation in the soil model
                       !is handled in CLM. Instead, the model here is designed to (simply)
                       !contain at least reasonable isotope ratios. Therefore, precipitation
                       !is added to both buckets. The larger bucket will have a longer
                       !turn over time and thus a longer memory/slower change due to precip and evap.
                       if (cam_in(c)%landfrac(i) .gt. 0._r8) then
                        
                         !!!=====================!!!
                         !!! shallow land bucket !!!
                         !!!=====================!!!
                         ! Add precipitation to shallow bucket:
                         sbuckmass = cam_out(c)%SbuckH(i) + prcsum*dtime !m/s*s = m

                         !Check to see if bucket has overflowed.
                         if (sbuckmass .gt. sbcklim) then
                           srunoff = sbuckmass - sbcklim !calculate runoff
                           sbuckmass = sbcklim           !set bucket to full
                         end if

                         sh2obckmass = sbuckmass !save bucket mass in order to calculate ratios

                         !calculate evap flux from shallow bucket using just soil evap.
                         !need to be careful with signs here!
                         if ((-x2a(index_x2a_Faxx_sevap,ig) - lakeflux) .gt. 0._r8) then
                           Re = wtrc_ratio(isph2o,sbuckmass,sh2obckmass)
                           shalet = Re*(-x2a(index_x2a_Faxx_sevap,ig) - lakeflux)
                         else !remove proportional amount of vapor! NO FRACTIONATION!
                           Re = wtrc_ratio(wtrc_species(wtrc_iasrfvap(j)),&
                                    cam_out(c)%qbot(i,wtrc_indices(wtrc_iasrfvap(j))),&
                                    cam_out(c)%qbot(i,wtrc_indices(wtrc_iasrfvap(1))))
                           shalet = Re*(-x2a(index_x2a_Faxx_sevap,ig) - lakeflux) 
                         end if
  
                         !set new bucket water mass value:
                         !NOTE: landevp [kg/m2/s] * dtime [seconds] / density [kg/m3] = m
                         sbuckmass = sbuckmass - shalet*dtime/1000._r8

                         !update shalet if it used too much water
                         if (sbuckmass .lt. 0.02_r8*sbcklim) then
                           !shalet = shalet - sbuckmass*1000._r8/dtime !sbuckmass < 0! 
                           !JN-fill 1/5th of bucket as a default.
                           sbuckmass = 0.2_r8*sbcklim            
                           empty_buck = .true.
                         end if
         
                         !save new bucket mass
                         cam_out(c)%SbuckH(i) = sbuckmass

                         !!!==================!!!
                         !!! deep land bucket !!!
                         !!!==================!!!
  
                         ! Retrieve deep bucket mass
                         dbuckmass = cam_out(c)%DbuckH(i) + prcsum*dtime !m/s*s = !
  
                         if (dbuckmass .gt. dbcklim) then
                           drunoff = dbuckmass - dbcklim
                           dbuckmass = dbcklim
                         end if

                         dh2obckmass = dbuckmass
                         ! calculate loss term from deep bucket. NOTE: intereception is 
                         ! included in this bucket - while conceptually being a bit weird,
                         ! this bucket is for non-frationationg fluxes so seems to make the most sense.
                         
                         if ((-x2a(index_x2a_Faxx_trans,ig) + -x2a(index_x2a_Faxx_inter,ig) + lakeflux) .ge. 0._r8) then 
                           !Calculate deep ET!
                           Re = wtrc_ratio(isph2o,dbuckmass,dh2obckmass)
                           deepet = Re*(-x2a(index_x2a_Faxx_trans,ig) + -x2a(index_x2a_Faxx_inter,ig) + lakeflux)
                         else
                           Re = wtrc_ratio(wtrc_species(wtrc_iasrfvap(j)),&
                                      cam_out(c)%qbot(i,wtrc_indices(wtrc_iasrfvap(j))),&
                                      cam_out(c)%qbot(i,wtrc_indices(wtrc_iasrfvap(1))))
                           deepet = Re*(-x2a(index_x2a_Faxx_trans,ig) + -x2a(index_x2a_Faxx_inter,ig) + lakeflux)
                         end if 

                         ! update deep bucket mass. 
                         dbuckmass = dbuckmass - deepet*dtime/1000._r8 ! [kg/m2/s]*s / density [kg/m3] = m
                            
                         !pull water from 'groundwater' to deep bucket if needed:
                         if (dbuckmass .lt. 0.0_r8) then 
                           deepet = deepet - dbuckmass*1000._r8/dtime !dbuckmass < 0                     
                           !following JN's code, fill to 20%
                           dbuckmass = 0.2_r8*dbcklim
                         end if
    
                         !Save new output bucket masses:
                         cam_out(c)%DbuckH(i) = dbuckmass
                       end if  ! landfrac > 0

                     end if !not ocean!
                   
                   cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)

                 end if ! j = 1
 
                 !-----------------
                 ! H216O
                 !-----------------
                 case (isph216o)
                   if (j .eq. 2) then ! only need to do once, don't need to repeat for j=48
                     !set wtrprsm and evapvars.
                     wtprsm = cam_out(c)%precrl_16O(i)+cam_out(c)%precsl_16O(i)+cam_out(c)%precrc_16O(i)+cam_out(c)%precsc_16O(i) 
                     shalet = 0._r8 !reinitialize shallow land evaporation
                     deepet = 0._r8 !reinitialize deep land evaporation
                     Re     = 0._r8
                     Rf     = 0._r8
                     R      = 0._r8
                     ovallt = 0._r8

                     !Calculate precip ratio as mass fixer:
                     Rp = wtrc_ratio(isph216o,prcsum,wtprsm) 

                     ! zero out ice bucket over land, and land buckets over ice.
                     if (cam_in(c)%icefrac(i) .eq. 1._r8 .or. cam_in(c)%landfrac(i) .eq. 0._r8) then
                       cam_out(c)%Sbuck16(i) = 0._r8
                       cam_out(c)%Dbuck16(i) = 0._r8
                     end if
                     !---------------------------------------------
                     !Apply bucket model for evaporation for H216O:
                     !---------------------------------------------
                     if (cam_in(c)%ocnfrac(i) .eq. 1._r8) then ! zero out buckets over the ocean
                       cam_out(c)%Sbuck16(i) = 0._r8
                       cam_out(c)%Dbuck16(i) = 0._r8
                     else
                       ! Do land buckets if not *entirely* over ocean and lf > 0
                       if (cam_in(c)%landfrac(i) .gt. 0._r8) then ! are we over land?

                         !!!================!!!
                         !!! shallow bucket !!!
                         !!!================!!!
  
                         !Add precipitation to shallow bucket:
                         sbuckmass = cam_out(c)%Sbuck16(i) + Rp*wtprsm*dtime !m/s*s = m

                         !Remove infiltration amount from shallow bucket
                         !RPF - I think these lines below are wrong...they don't achieve what's necessary.
                         Re = wtrc_ratio(isph216o, sbuckmass, sbcklim+srunoff)
                         sbuckmass = sbuckmass - Re*srunoff

                         ! save buckmass for preserving ratio of shallow bucket if empty.
                         ovallt = sbuckmass

                         ! Set fractionation factors
                         Re = wtrc_ratio(isph216o,sbuckmass,sh2obckmass)
                         aeq = wiso_alpl(isph216o,cam_in(c)%ts(i))
                         if (is_first_step() .or. is_first_restart_step()) then !rho, zsrf, u* == 0 @ first step...
                           akin = 1.0_r8 ! set as 1 if first step.
                         else ! call more robust routine                       
                           call wiso_kmol_land(isph216o,rho,zsrf,ustar,z0m,akin) 
                         end if 

                         !Shallow ET (soil evap)
                         ! need to be careful about sign!
                         !nothing to be done here w/ isotope fractioantion - it's h216o
                         if ((-x2a(index_x2a_Faxx_sevap,ig) - lakeflux) .gt. 0.0_r8) then
                           if (abs(evaprh_ts - 1._r8) .gt. 0.02_r8) then 
                             R = wtrc_ratio(isph216o,&
                                            cam_out(c)%qbot(i,wtrc_indices(wtrc_iasrfvap(j))),&
                                            cam_out(c)%qbot(i,wtrc_indices(wtrc_iasrfvap(1))))
                             !Rf = akin*nwgt*(Re/aeq - evaprh_ts*R)/(1.0_r8 - evaprh_ts)
                             Rf = Re/aeq
                           else 
                             Rf = Re/aeq
                           end if
                           shalet = Rf*(-x2a(index_x2a_Faxx_sevap,ig) - lakeflux) 
                         else
                           Re = wtrc_ratio(wtrc_species(wtrc_iasrfvap(j)),&
                                           cam_out(c)%qbot(i,wtrc_indices(wtrc_iasrfvap(j))),&
                                           cam_out(c)%qbot(i,wtrc_indices(wtrc_iasrfvap(1))))
                           shalet = Re*(-x2a(index_x2a_Faxx_sevap,ig) - lakeflux) 
                         end if 
 
                         !update bucket masses
                         sbuckmass = sbuckmass - shalet*dtime/1000._r8
 
                         !adjust fluxes if bucket has become negative
                         if (empty_buck) then
                           ! set to bulk water
                           sbuckmass = cam_out(c)%SbuckH(i)
                         end if 

                         ! save H216O bucket masses.
                         cam_out(c)%Sbuck16(i) = sbuckmass

                         !!!=============!!!
                         !!! deep bucket !!!
                         !!!=============!!!

                         !Retreive deep bucket mass:
                         dbuckmass = cam_out(c)%Dbuck16(i) + Rp*wtprsm*dtime

                         !remove runoff amount from deep bucket
                         Re = wtrc_ratio(isph216o, dbuckmass, dbcklim+drunoff)
                         dbuckmass = dbuckmass - Re*drunoff

                         ! Deep ET (interception and transpiration)
                         if ((-x2a(index_x2a_Faxx_inter,ig) + -x2a(index_x2a_Faxx_trans,ig) + lakeflux) .ge. 0._r8) then
                           Re = wtrc_ratio(isph216o,dbuckmass,dh2obckmass)
                           deepet = Re*(-x2a(index_x2a_Faxx_inter,ig) + -x2a(index_x2a_Faxx_trans,ig) + lakeflux)
                         else
                           Re = wtrc_ratio(wtrc_species(wtrc_iasrfvap(j)),&
                                           cam_out(c)%qbot(i,wtrc_indices(wtrc_iasrfvap(j))),&
                                           cam_out(c)%qbot(i,wtrc_indices(wtrc_iasrfvap(1))))
                           deepet = Re*(-x2a(index_x2a_Faxx_inter,ig) + -x2a(index_x2a_Faxx_trans,ig) + lakeflux)
                         end if

                         !update bucket masses
                         dbuckmass = dbuckmass - deepet*dtime/1000._r8
  
                         !adjust fluxes if bucket has become negative
                         if (dbuckmass .lt. 0._r8) then
                           deepet = deepet + dbuckmass*1000._r8/dtime ! note dbuckmass < 0
                           ! set back to bulk water
                           dbuckmass = cam_out(c)%DbuckH(i) 
                         end if

                         ! save H216O bucket masses.
                         cam_out(c)%Dbuck16(i) = dbuckmass
                       end if  
  
                     end if ! land present

                   end if
                   !assign cflx for hdo
                   cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = shalet + deepet + & ! shalet + deepet + & ! land values
                                                                     -x2a(index_x2a_Faxx_evap_16O,ig)


                 case (isphdo)
                   !reset wtprsm, shalet, deepet vvars for hdo
                   wtprsm = cam_out(c)%precrl_HDO(i)+cam_out(c)%precsl_HDO(i)+cam_out(c)%precrc_HDO(i)+cam_out(c)%precsc_HDO(i) 
                   shalet = 0._r8 !reinitialize shallow land evaporation
                   deepet = 0._r8 !reinitialize deep land evaporation
                   Re     = 0._r8
                   Rf     = 0._r8
                   R      = 0._r8
                   akin   = 0._r8
                   aeq    = 0._r8
                   ovalhv = 0._r8

                   ! zero out ice bucket over land, and land buckets over ice.
                   if (cam_in(c)%icefrac(i) .eq. 1._r8 .or. cam_in(c)%landfrac(i) .eq. 0._r8) then
                     cam_out(c)%SbuckD(i) = 0._r8
                     cam_out(c)%DbuckD(i) = 0._r8
                   end if

                   ! are we on land? if so, do bucket calculations.
                   if (cam_in(c)%ocnfrac(i) .eq. 1._r8) then
                     ! zero out bucket over the ocean, since this could cause issues for areas
                     ! that flip between ice and ocean.
                     cam_out(c)%SbuckD(i) = 0._r8
                     cam_out(c)%DbuckD(i) = 0._r8
                   else
                     !---------------------------------------------
                     !Apply bucket model for evaporation for H216O:
                     !---------------------------------------------
                     if (cam_in(c)%landfrac(i) .gt. 0._r8) then ! are we over land?

                       !!!================!!!
                       !!! shallow bucket !!!
                       !!!================!!!
                       !Add precipitation to shallow bucket:
                       sbuckmass = cam_out(c)%SbuckD(i) + Rp*wtprsm*dtime !m/s*s = m

                       !Remove runoff amount from shallow bucket
                       Re = wtrc_ratio(isphdo, sbuckmass, sbcklim+srunoff)
                       sbuckmass = sbuckmass - Re*srunoff

                       !save bucket mass for fixer at end if bucket goes empty.
                       ovalhv = sbuckmass

                       ! Set fractionation factors. 
                       aeq = wiso_alpl(isphdo,cam_in(c)%ts(i))
                       if (is_first_step() .or. is_first_restart_step()) then !rho, zsrf, u* == 0 @ first step...
                         akin = (1.0_r8/1.0251_r8)**0.5_r8 ! set to merlivat di/d^0.5 for first step.
                       else ! call more robust routine                       
                         call wiso_kmol_land(isphdo,rho,zsrf,ustar,z0m,akin) 
                       end if 

                       R = wtrc_ratio(isphdo,cam_out(c)%qbot(i,wtrc_indices(wtrc_iasrfvap(j))),&
                                      cam_out(c)%qbot(i,wtrc_indices(wtrc_iasrfvap(1))))

                       !Shallow ET (soil evap)
                       !need to be careful w/ sign of evap flux!
                       if ((-x2a(index_x2a_Faxx_sevap,ig) - lakeflux) .gt. 0.0_r8) then
                         Re = wtrc_ratio(isphdo,sbuckmass,sh2obckmass)
                         if (lfopt .eq. 0) then ! if fractionation is turned off...
                           shalet = Re*(-x2a(index_x2a_Faxx_sevap,ig) - lakeflux)
                         else ! make soils fractionate...
                           !per suggestion by DCN, only use alpha for T>0
                           if (cam_in(c)%ts(i) .gt. 273.15_r8) then
                             !simple CG model:
                             if (abs(evaprh_ts - 1._r8) .gt. 0.02_r8) then
                               Rf = akin*(Re/aeq - evaprh_ts*R)/(1._r8 - evaprh_ts)
                             else
                               Rf = Re/aeq
                             end if
                           else ! no fractionation below freezing!
                             Rf = Re
                           end if
                           !Multiply ratio by flux:
                           shalet = Rf*(-x2a(index_x2a_Faxx_sevap,ig) - lakeflux)
                         end if
                       else ! condensation occuring
                         Re = wtrc_ratio(wtrc_species(wtrc_iasrfvap(j)),&
                                   cam_out(c)%qbot(i,wtrc_indices(wtrc_iasrfvap(j))),&
                                   cam_out(c)%qbot(i,wtrc_indices(wtrc_iasrfvap(1))))
                         shalet = Re*(-x2a(index_x2a_Faxx_sevap,ig) - lakeflux)
                       end if

                       !update bucket mass
                       sbuckmass = sbuckmass - shalet*dtime/1000._r8

                       !adjust fluxes if bucket has become negative
                       if (empty_buck) then
                         R = wtrc_ratio(isphdo,ovalhv,ovallt)
                         sbuckmass = R*cam_out(c)%SbuckH(i)
                       end if

                       !save new hdo bucket mass 
                       cam_out(c)%SbuckD(i) = sbuckmass

                       !!!=============!!!
                       !!! deep bucket !!!
                       !!!=============!!!
  
                       !Retreive deep bucket mass:
                       dbuckmass = cam_out(c)%DbuckD(i) + Rp*wtprsm*dtime
                     
                       !remove runoff amount from deep bucket
                       Re = wtrc_ratio(isphdo, dbuckmass, dbcklim+drunoff)
                       dbuckmass = dbuckmass - Re*drunoff

                       !Calculate evaporation flux from deep bucket
                       ! Deep ET (interception and transpiration)
                       if ((-x2a(index_x2a_Faxx_inter,ig) + -x2a(index_x2a_Faxx_trans,ig) + lakeflux) .gt. 0._r8) then
                         Re = wtrc_ratio(isphdo,dbuckmass,dh2obckmass)
                         if (lfopt .eq. 0) then
                           deepet = Re*(-x2a(index_x2a_Faxx_inter,ig) + -x2a(index_x2a_Faxx_trans,ig) + lakeflux)
                         else
                           if (cam_in(c)%ts(i) > 273.15_r8) then
                             if (abs(evaprh_ts - 1.0_r8) .gt. 0.02_r8) then
                               R = wtrc_ratio(isphdo,cam_out(c)%qbot(i,wtrc_indices(wtrc_iasrfvap(j))),&
                                              cam_out(c)%qbot(i,wtrc_indices(wtrc_iasrfvap(1))))
                               Rf = akin*(Re/aeq - evaprh_ts*R)/(1.0_r8 - evaprh_ts)
                             else
                               Rf = Re/aeq
                             end if
                           else ! below freezing.
                             Rf = Re
                           end if
                           deepet = Re*(-x2a(index_x2a_Faxx_inter,ig) + -x2a(index_x2a_Faxx_trans,ig)) + Rf*lakeflux
                         end if ! lfopt.
                       else 
                         Re = wtrc_ratio(wtrc_species(wtrc_iasrfvap(j)),&
                                   cam_out(c)%qbot(i,wtrc_indices(wtrc_iasrfvap(j))),&
                                   cam_out(c)%qbot(i,wtrc_indices(wtrc_iasrfvap(1))))
                         deepet = Re*(-x2a(index_x2a_Faxx_inter,ig) + -x2a(index_x2a_Faxx_trans,ig) + lakeflux)
                       end if
  
                       !update bucket masses
                       dbuckmass = dbuckmass - deepet*dtime/1000._r8
                       
                       !adjust fluxes if bucket has become negative
                       if (dbuckmass .lt. 0._r8) then
                         deepet = deepet + dbuckmass*1000._r8/dtime ! note dbuckmass < 0
                         !groundwater pulled up so that d2H of bucket = predicted d2H
                         if (gwopt .eq. 1) then !JN's GNIP parameterization
                           R = 0.97623_r8 - 0.001976846_r8*abs(wtlat)
                         else if (gwopt .eq. 2) then !Bowen2010 AREPS
                           R = 1._r8 + (0.562*abs(wtlat) - 0.0338*wtlat**2 - 0.0136*cam_out(c)%topo(i))/1000._r8
                         else if (gwopt .eq. 3) then !IsoMAP, July 20, 2020                       
                           R = 1._r8 + (-27.3631_r8 + 1.43709_r8*abs(wtlat) - 0.04349_r8*wtlat**2 &
                                       - 0.0169_r8*cam_out(c)%topo(i))/1000._r8
                         else if (gwopt .eq. 4) then !SMOW - for comparison w/ Nick's 2D model
                           R = 1._r8
                         end if
                         ! set to groundwater isotope ratio
                         dbuckmass = R*cam_out(c)%DbuckH(i)
                      end if

                       ! save HDO deep bucket mass.
                       cam_out(c)%DbuckD(i) = dbuckmass
                     end if
 
                   end if ! land present
                   !assign cflx for hdo
                   cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = shalet + deepet + & ! shalet + deepet + & ! land values
                                                                     -x2a(index_x2a_Faxx_evap_HDO,ig)

                 case (isph218o)
                   !set wtrprsm and evapvars.
                   wtprsm = cam_out(c)%precrl_18O(i)+cam_out(c)%precsl_18O(i)+cam_out(c)%precrc_18O(i)+cam_out(c)%precsc_18O(i) 
                   shalet = 0._r8 !reinitialize shallow land evaporation
                   deepet = 0._r8 !reinitialize deep land evaporation
                   Re     = 0._r8
                   Rf     = 0._r8
                   akin   = 0._r8
                   aeq    = 0._r8
                   
                   ! zero out ice bucket over land, and land buckets over ice. 
                   if (cam_in(c)%icefrac(i) .eq. 1._r8 .or. cam_in(c)%landfrac(i) .eq. 0._r8) then
                     cam_out(c)%Sbuck18(i) = 0._r8
                     cam_out(c)%Dbuck18(i) = 0._r8
                   end if

                   ! are we on land? if so, do bucket calculations.
                   if (cam_in(c)%ocnfrac(i) .eq. 1._r8) then
                     ! zero out bucket over the ocean, since this could cause issues for areas
                     ! that flip between ice and ocean.
                     cam_out(c)%Sbuck18(i) = 0._r8
                     cam_out(c)%Dbuck18(i) = 0._r8
                   else
                     !---------------------------------------------
                     !Apply bucket model for evaporation for H218O:
                     !---------------------------------------------
                     if (cam_in(c)%landfrac(i) .gt. 0._r8) then ! are we over land?

                       !Add precipitation to shallow bucket:
                       sbuckmass = cam_out(c)%Sbuck18(i) + Rp*wtprsm*dtime !m/s*s = m

                       !Remove infiltration amount from shallow bucket
                       Re = wtrc_ratio(isph218o, sbuckmass, sbcklim+srunoff)
                       sbuckmass = sbuckmass - Re*srunoff
                     
                       !save bucket mass for fixer at end if bucket goes empty.
                       ovalhv = sbuckmass

                       ! Set fractionation factors:                     
                       aeq = wiso_alpl(isph218o,cam_in(c)%ts(i))
                       if (is_first_step() .or. is_first_restart_step()) then !rho, zsrf, u* == 0 @ first step...
                         akin = (1.0_r8/1.0283_r8)**0.5_r8 ! set to merlivat Di/D^0.5 for first step.
                       else ! call more robust routine                       
                         call wiso_kmol_land(isph218o,rho,zsrf,ustar,z0m,akin) 
                       end if 

                       R = wtrc_ratio(isph218o,cam_out(c)%qbot(i,wtrc_indices(wtrc_iasrfvap(j))),&
                                      cam_out(c)%qbot(i,wtrc_indices(wtrc_iasrfvap(1))))

                       !Shallow ET (soil evap)
                       if ((-x2a(index_x2a_Faxx_sevap,ig) - lakeflux) .gt. 0.0_r8) then
                         Re = wtrc_ratio(isph218o,sbuckmass,sh2obckmass)
                         if (lfopt .eq. 0) then ! if fractionation is turned off...
                           shalet = Re*-x2a(index_x2a_Faxx_sevap,ig)
                         else ! make soils fractionate...
                           !per suggestion by DCN, only use alpha for T>0
                           if (cam_in(c)%ts(i) .gt. 273.15_r8) then
                             !simple CG model:
                             if (abs(evaprh_ts - 1._r8) .gt. 0.02_r8) then 
                               Rf = akin*(Re/aeq - evaprh_ts*R)/(1._r8-evaprh_ts)
                             else
                               Rf = Re/aeq
                             end if
                           else ! no fractionation below freezing!
                             Rf = Re
                           end if
                           shalet = Rf*(-x2a(index_x2a_Faxx_sevap,ig) - lakeflux)
                         end if
                       else !condensation occuring!
                         Re = wtrc_ratio(wtrc_species(wtrc_iasrfvap(j)),&
                                   cam_out(c)%qbot(i,wtrc_indices(wtrc_iasrfvap(j))),&
                                   cam_out(c)%qbot(i,wtrc_indices(wtrc_iasrfvap(1))))
                         shalet = Re*(-x2a(index_x2a_Faxx_sevap,ig) - lakeflux)
                       end if

                       !update bucket masses
                       sbuckmass = sbuckmass - shalet*dtime/1000._r8
 
                       !adjust fluxes if bucket has become negative
                       if (empty_buck) then
                         R = wtrc_ratio(isph218o,ovalhv,ovallt)
                         sbuckmass = R*cam_out(c)%SbuckH(i)
                       end if
                       
                       ! save H218O bucket mass.
                       cam_out(c)%Sbuck18(i) = sbuckmass

                       !!!=============!!!
                       !!! deep bucket !!!
                       !!!=============!!!

                       !Retreive deep bucket mass:
                       dbuckmass = cam_out(c)%Dbuck18(i) + Rp*wtprsm*dtime

                       !remove runoff amount from deep bucket
                       Re = wtrc_ratio(isph218o, dbuckmass, dbcklim+drunoff)
                       dbuckmass = dbuckmass - Re*drunoff

                       ! Deep ET (interception and transpiration)
                       if ((-x2a(index_x2a_Faxx_inter,ig) + -x2a(index_x2a_Faxx_trans,ig) + lakeflux) .gt. 0.0_r8) then 
                         Re = wtrc_ratio(isph218o,dbuckmass,dh2obckmass)
                         if (lfopt .eq. 0) then
                           deepet = Re*(-x2a(index_x2a_Faxx_inter,ig) + -x2a(index_x2a_Faxx_trans,ig) + lakeflux)
                         else
                           if (cam_in(c)%ts(i) > 273.15_r8) then
                             if (abs(evaprh_ts - 1.0_r8) .gt. 0.02_r8) then
                               R = wtrc_ratio(isph218o,cam_out(c)%qbot(i,wtrc_indices(wtrc_iasrfvap(j))),&
                                              cam_out(c)%qbot(i,wtrc_indices(wtrc_iasrfvap(1))))
                               Rf = akin*(Re/aeq - evaprh_ts*R)/(1.0_r8 - evaprh_ts)
                             else
                               Rf = Re/aeq
                             end if
                           else ! below freezing.
                             Rf = Re
                           end if
                           deepet = Re*(-x2a(index_x2a_Faxx_inter,ig) + -x2a(index_x2a_Faxx_trans,ig)) + Rf*lakeflux
                         end if ! lfopt.
                       else 
                         Re = wtrc_ratio(wtrc_species(wtrc_iasrfvap(j)),&
                                   cam_out(c)%qbot(i,wtrc_indices(wtrc_iasrfvap(j))),&
                                   cam_out(c)%qbot(i,wtrc_indices(wtrc_iasrfvap(1))))
                         deepet = Re*(-x2a(index_x2a_Faxx_inter,ig) + -x2a(index_x2a_Faxx_trans,ig) + lakeflux)
                       end if

                       !update bucket masses
                       dbuckmass = dbuckmass - deepet*dtime/1000._r8
                       
                       !adjust fluxes if bucket has become negative
                       if (dbuckmass .lt. 0._r8) then
                         deepet = deepet + dbuckmass*1000._r8/dtime ! note dbuckmass < 0
                        !groundwater pulled up so that d18o of bucket = predicted d18o
                         if (gwopt .eq. 1) then !JN's GNIP parameterization
                           R = 0.99607_r8 - 0.0002494177_r8*abs(wtlat)
                         else if (gwopt .eq. 2) then !Bowen2010 AREPS
                           R = 1._r8 + (0.078*abs(wtlat) - 0.00428*wtlat**2 - 0.00194*cam_out(c)%topo(i))/1000._r8
                         else if (gwopt .eq. 3) then !IsoMAP, July 20, 2020                       
                           R = 1._r8 + (-4.73842_r8 + 0.17409_r8*abs(wtlat) - 0.00524525_r8*wtlat**2 &
                                        - 0.00218674_r8*cam_out(c)%topo(i))/1000._r8
                         else if (gwopt .eq. 4) then !SMOW - for comparison w/ Nick's 2D model
                           R = 1._r8
                         end if
                         ! set to groundwater isotope ratio
                         dbuckmass = R*cam_out(c)%DbuckH(i)
                       end if

                       ! save H218O bucket masses.
                       cam_out(c)%Dbuck18(i) = dbuckmass
                     end if 

                   end if ! land present
                   !assign cflx for h218o
                   cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = shalet + deepet + -x2a(index_x2a_Faxx_evap_18O,ig)

                 case (isph217o)
                  cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap_17O,ig)
                  ! RPF - no bucket here for 17O

               case (isphto)
                 cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap_HTO,ig) 
                 !RPF - no bucket for HTO, but should be off anyway!

               end select
             end do
          end if

          cam_in(c)%asdir(i)     =  x2a(index_x2a_Sx_avsdr, ig)
          cam_in(c)%aldir(i)     =  x2a(index_x2a_Sx_anidr, ig)
          cam_in(c)%asdif(i)     =  x2a(index_x2a_Sx_avsdf, ig)
          cam_in(c)%aldif(i)     =  x2a(index_x2a_Sx_anidf, ig)
          cam_in(c)%sst(i)       =  x2a(index_x2a_So_t,     ig)
          cam_in(c)%snowhland(i) =  x2a(index_x2a_Sl_snowh, ig)
          cam_in(c)%snowhice(i)  =  x2a(index_x2a_Si_snowh, ig)

          if ( associated(cam_in(c)%ram1) ) &
               cam_in(c)%ram1(i) =  x2a(index_x2a_Sl_ram1 , ig)
          if ( associated(cam_in(c)%fv) ) &
               cam_in(c)%fv(i)   =  x2a(index_x2a_Sl_fv   , ig)
         ! if ( associated(cam_in(c)%z0m) ) & !RPF _ roughness length
         !      cam_in(c)%z0m(i)   =  x2a(index_x2a_Sl_z0m   , ig)
         !      write(*,*) 'roughness length from clm:', cam_in(c)%z0m(i)
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
    use water_isotopes   , only: isph2o, isph216o, isphdo, isph218o, isph217o, isphto

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
    logical :: pass16, passD, pass18, pass17, passT !logicals that prevent the passing of water tag infromation to iCLM4.
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
            a2x(index_a2x_Sa_shum_17O   ,ig) = 0._r8
            a2x(index_a2x_Sa_shum_HTO   ,ig) = 0._r8

           !logical to prevent surface vapor from tags being passed on. -JN
            pass16 = .true.
            passD  = .true.
            pass18 = .true.
            pass17 = .true.
            passT  = .true. 

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
                case (isph217o)
                  if(pass17) then !pass on H217O?
                    a2x(index_a2x_Sa_shum_17O   ,ig) = cam_out(c)%qbot(i,wtrc_indices(wtrc_iasrfvap(j)))
                    pass17 = .false.
                  end if
                case (isphto)
                  if(passT) then !pass on HTO?
                    a2x(index_a2x_Sa_shum_HTO   ,ig) = cam_out(c)%qbot(i,wtrc_indices(wtrc_iasrfvap(j)))
                    passT = .false.
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
            a2x(index_a2x_Faxa_rainl_17O,ig)=cam_out(c)%precrl_17O(i)*1000._r8
            a2x(index_a2x_Faxa_snowl_17O,ig)=cam_out(c)%precsl_17O(i)*1000._r8
            a2x(index_a2x_Faxa_rainc_17O,ig)=cam_out(c)%precrc_17O(i)*1000._r8
            a2x(index_a2x_Faxa_snowc_17O,ig)=cam_out(c)%precsc_17O(i)*1000._r8
            a2x(index_a2x_Faxa_rainl_HTO,ig)=cam_out(c)%precrl_HTO(i)*1000._r8
            a2x(index_a2x_Faxa_snowl_HTO,ig)=cam_out(c)%precsl_HTO(i)*1000._r8
            a2x(index_a2x_Faxa_rainc_HTO,ig)=cam_out(c)%precrc_HTO(i)*1000._r8
            a2x(index_a2x_Faxa_snowc_HTO,ig)=cam_out(c)%precsc_HTO(i)*1000._r8
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
