!
! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!  ocn_forward_mode
!
!> \brief Main driver for MPAS ocean core
!> \author Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date   September 2011
!> \details
!> This module contains initialization and timestep drivers for the MPAS ocean core.
!
MODULE ocn_forward_mode

   USE mpas_kind_types
   USE mpas_derived_types
   USE mpas_pool_routines
   USE mpas_stream_manager
   USE mpas_timekeeping
   USE mpas_dmpar
   USE mpas_timer
   USE mpas_log
   USE mpas_decomp
   USE mpas_log

   USE ocn_analysis_driver
   USE ocn_init_routines

   USE ocn_time_integration
   USE ocn_time_integration_split
   USE ocn_tendency
   USE ocn_diagnostics
   USE ocn_test

   USE ocn_thick_hadv
   USE ocn_thick_vadv
   USE ocn_thick_ale
   USE ocn_thick_surface_flux

   USE ocn_vel_pressure_grad
   USE ocn_vel_vadv
   USE ocn_vel_hmix
   USE ocn_vel_forcing
   USE ocn_vel_coriolis
   USE ocn_vel_forcing_surface_stress
   USE ocn_surface_bulk_forcing
   USE ocn_surface_land_ice_fluxes
   USE ocn_frazil_forcing

   USE ocn_tracer_hmix
   USE ocn_tracer_hmix_redi
   USE ocn_tracer_surface_flux_to_tend
   USE ocn_tracer_short_wave_absorption
   USE ocn_tracer_short_wave_absorption_variable
   USE ocn_tracer_nonlocalflux
   USE ocn_tracer_advection
   USE ocn_tracer_ecosys
   USE ocn_tracer_DMS
   USE ocn_tracer_MacroMolecules
   USE ocn_tracer_surface_restoring
   USE ocn_gm

   USE ocn_high_freq_thickness_hmix_del2

   USE ocn_equation_of_state

   USE ocn_vmix

   USE ocn_forcing
   USE ocn_sea_ice

   USE ocn_constants

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: ocn_forward_mode_init, ocn_forward_mode_run, ocn_forward_mode_finalize, &
             ocn_forward_mode_setup_clock

   CONTAINS
   !
   !  FUNCTION ocn_forward_mode_init
   !
   !> \brief   Initialize MPAS-Ocean core
   !> \author  Doug Jacobsen, Mark Petersen, Todd Ringler
   !> \date    September 2011
   !> \details
   !> This function calls all initializations required to begin a
   !> simulation with MPAS-Ocean
   !
   FUNCTION ocn_forward_mode_init(domain, startTimeStamp) RESULT(ierr)
      
      IMPLICIT NONE

      TYPE(domain_type), INTENT(INOUT) :: domain
      CHARACTER(LEN=*), intent(out) :: startTimeStamp
      INTEGER :: ierr
      
      ! Local variables      
      CHARACTER(LEN=StrKIND), POINTER &
      :: xtime, simulationStartTime, config_vert_coord_movement, &
         config_pressure_gradient_type, &
         config_surface_salinity_monthly_restoring_compute_interval
      INTEGER :: err_tmp
      INTEGER, POINTER :: nVertLevels      
      LOGICAL, POINTER :: config_do_restart, config_read_nearest_restart, &
                          config_filter_btr_mode, config_conduct_tests, &
                          config_use_surface_salinity_monthly_restoring            
      REAL(KIND=RKIND) :: dt, maxDensity, maxDensity_global
      REAL(KIND=RKIND), POINTER :: daysSinceStartOfSim, config_maxMeshDensity
      REAL(KIND=RKIND), DIMENSION(:), POINTER :: meshDensity
      TYPE(mpas_pool_type), POINTER :: meshPool, diagnosticsPool
      TYPE(block_type), POINTER :: block
      TYPE(MPAS_Time_type) :: xtime_timeType, simulationStartTime_timeType, startTime, &
                              alarmTime
      TYPE(MPAS_TimeInterval_type) :: timeStep, alarmTimeStep

      ierr = 0
      !
      ! Set startTimeStamp based on the start time of the simulation clock
      !
      startTime = mpas_get_clock_time(domain % clock, MPAS_START_TIME, err_tmp)
      CALL mpas_get_time(startTime, dateTimeString=startTimeStamp)
      ierr = ior(ierr, err_tmp)

      ! Setup ocean config pool
      CALL ocn_constants_init(domain % configs, domain % packages)

      CALL mpas_pool_get_config(domain % configs, 'config_do_restart', config_do_restart)
      CALL mpas_pool_get_config(domain % configs, 'config_read_nearest_restart', &
                                config_read_nearest_restart)
      CALL mpas_pool_get_config(domain % configs, 'config_vert_coord_movement', &
                                config_vert_coord_movement)
      CALL mpas_pool_get_config(domain % configs, 'config_pressure_gradient_type', &
                                config_pressure_gradient_type)
      CALL mpas_pool_get_config(domain % configs, 'config_filter_btr_mode', &
                                config_filter_btr_mode)
      CALL mpas_pool_get_config(domain % configs, 'config_maxMeshDensity', &
                                config_maxMeshDensity)
      CALL mpas_pool_get_config(domain % configs, 'config_conduct_tests', &
                                config_conduct_tests)
      CALL mpas_pool_get_config(domain % configs, &
                                'config_use_surface_salinity_monthly_restoring', &
                                config_use_surface_salinity_monthly_restoring)
      CALL mpas_pool_get_config &
      (domain % configs, 'config_surface_salinity_monthly_restoring_compute_interval', &
       config_surface_salinity_monthly_restoring_compute_interval)
      !
      ! Read input data for model
      !
      CALL mpas_timer_start('io_read')
      CALL mpas_stream_mgr_read(domain % streamManager, streamID='mesh', &
                                whence=MPAS_STREAM_NEAREST, ierr=err_tmp)

      IF (config_do_restart) THEN
         IF (config_read_nearest_restart) THEN
            CALL mpas_stream_mgr_read(domain % streamManager, streamID='restart', &
                                      whence=MPAS_STREAM_NEAREST, ierr=err_tmp)
         ELSE
            CALL mpas_stream_mgr_read(domain % streamManager, streamID='restart', &                           
                                      ierr=err_tmp)
         END IF
      ELSE
         CALL mpas_stream_mgr_read(domain % streamManager, streamID='input', &
                                   ierr=err_tmp)
      END IF

      CALL ocn_analysis_bootstrap(domain, err=err_tmp)

      CALL mpas_timer_stop('io_read')
      CALL mpas_timer_start('reset_io_alarms')
      CALL mpas_stream_mgr_reset_alarms(domain % streamManager, streamID='input', &                     
                                        ierr=err_tmp)
      CALL mpas_stream_mgr_reset_alarms(domain % streamManager, streamID='restart', &
                                        ierr=err_tmp)
      CALL mpas_stream_mgr_reset_alarms(domain % streamManager, &
                                        direction=MPAS_STREAM_OUTPUT, ierr=err_tmp)
      CALL mpas_timer_stop('reset_io_alarms')

      ! Read the remaining input streams
      CALL mpas_timer_start('io_read')
      CALL mpas_stream_mgr_read(domain % streamManager, ierr=err_tmp)
      ierr = ior(ierr, err_tmp)
      CALL mpas_timer_stop('io_read')
      CALL mpas_timer_start('reset_io_alarms')
      CALL mpas_stream_mgr_reset_alarms(domain % streamManager, &
                                        direction=MPAS_STREAM_INPUT, ierr=err_tmp)
      ierr = ior(ierr, err_tmp)
      CALL mpas_timer_stop('reset_io_alarms')

      ! Initialize submodules before initializing blocks.
      CALL ocn_timestep_init(ierr)

      CALL ocn_thick_hadv_init(err_tmp)
      ierr = ior(ierr, err_tmp)
      CALL ocn_thick_vadv_init(err_tmp)
      ierr = ior(ierr, err_tmp)
      CALL ocn_thick_surface_flux_init(err_tmp)
      ierr = ior(ierr, err_tmp)
      CALL ocn_thick_ale_init(err_tmp)
      ierr = ior(ierr,err_tmp)

      CALL ocn_vel_coriolis_init(err_tmp)
      ierr = ior(ierr, err_tmp)
      CALL ocn_vel_pressure_grad_init(err_tmp)
      ierr = ior(ierr, err_tmp)
      CALL ocn_vel_vadv_init(err_tmp)
      ierr = ior(ierr, err_tmp)
      CALL ocn_vel_hmix_init(err_tmp)
      ierr = ior(ierr, err_tmp)
      CALL ocn_vel_forcing_init(err_tmp)
      ierr = ior(ierr, err_tmp)
      CALL ocn_vel_forcing_surface_stress_init(err_tmp)
      ierr = ior(ierr, err_tmp)
      CALL ocn_surface_bulk_forcing_init(err_tmp)
      ierr = ior(ierr, err_tmp)
      CALL ocn_surface_land_ice_fluxes_init(err_tmp)
      ierr = ior(ierr, err_tmp)
      CALL ocn_frazil_forcing_init(err_tmp)
      ierr = ior(ierr, err_tmp)

      CALL ocn_tracer_hmix_init(err_tmp)
      ierr = ior(ierr, err_tmp)
      CALL ocn_tracer_hmix_redi_init(err_tmp)
      ierr = ior(ierr, err_tmp)
      CALL ocn_tracer_surface_flux_init(err_tmp)
      ierr = ior(ierr, err_tmp)
      CALL ocn_tracer_advection_init(err_tmp)
      ierr = ior(ierr,err_tmp)
      CALL ocn_tracer_short_wave_absorption_init(domain,err_tmp)
      ierr = ior(ierr,err_tmp)
      CALL ocn_gm_init(err_tmp)
      ierr = ior(ierr,err_tmp)
      CALL ocn_tracer_nonlocalflux_init(err_tmp)
      ierr = ior(ierr,err_tmp)
      CALL ocn_tracer_ecosys_init(domain, err_tmp)
      ierr = ior(ierr,err_tmp)
      CALL ocn_tracer_DMS_init(domain, err_tmp)
      ierr = ior(ierr,err_tmp)
      CALL ocn_tracer_MacroMolecules_init(domain, err_tmp)
      ierr = ior(ierr,err_tmp)

      CALL ocn_vmix_init(domain, err_tmp)
      ierr = ior(ierr, err_tmp)

      CALL ocn_equation_of_state_init(err_tmp)
      ierr = ior(ierr, err_tmp)

      CALL ocn_tendency_init(err_tmp)
      ierr = ior(ierr,err_tmp)
      CALL ocn_diagnostics_init(err_tmp)
      ierr = ior(ierr,err_tmp)

      CALL ocn_forcing_init(err_tmp)
      ierr = ior(ierr,err_tmp)

      CALL ocn_high_freq_thickness_hmix_del2_init(err_tmp)
      ierr = ior(ierr,err_tmp)

      CALL mpas_pool_get_DIMENSION(domain % blocklist % dimensions, 'nVertLevels', &
                                   nVertLevels)
      CALL ocn_sea_ice_init(nVertLevels, err_tmp)
      ierr = ior(ierr, err_tmp)

      IF (ierr == 1) THEN
          CALL mpas_log_write &
          ('An error was encountered while initializing the MPAS-Ocean forward mode', &
           MPAS_LOG_CRIT)
      END IF

      CALL ocn_init_metadata(domain)

      CALL ocn_init_routines_vert_coord(domain)

      CALL ocn_init_routines_compute_max_level(domain)

      CALL ocn_time_integration_split_init(domain)

      CALL mpas_log_write(' Vertical coordinate movement is: ' // &
                          TRIM(config_vert_coord_movement))

      IF (config_vert_coord_movement /= 'fixed' .AND. &
          config_vert_coord_movement /= 'uniform_stretching' .AND. &
          config_vert_coord_movement /= 'impermeable_interfaces' .AND. &
          config_vert_coord_movement /= 'user_specified') THEN
         CALL mpas_log_write(' Incorrect choice of config_vert_coord_movement.')
         CALL mpas_log_write('Incorrect choice of config_vert_coord_movement.', &
                             MPAS_LOG_CRIT)
      END IF

      IF (config_vert_coord_movement /= 'impermeable_interfaces' &
          .AND. config_pressure_gradient_type == 'MontgomeryPotential') THEN
         CALL mpas_log_write( &
         'Incorrect combination of config_vert_coord_movement and config_pressure_gradient_type', MPAS_LOG_CRIT)
      END IF

      IF (config_filter_btr_mode.AND. &
          config_vert_coord_movement /= 'fixed')then
         CALL mpas_log_write('filter_btr_mode has only been tested with config_vert_coord_movement=fixed.', MPAS_LOG_CRIT)
      END IF

      ! find the maximum value of the meshDensity
      IF (config_maxMeshDensity < 0.0_RKIND) THEN
        maxDensity = -1
        block => domain % blocklist
        DO WHILE (ASSOCIATED(block))
          CALL mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
          CALL mpas_pool_get_array(meshPool, 'meshDensity', meshDensity)
          maxDensity = max(maxDensity, maxval(meshDensity))
          block => block % next
        END DO
        CALL mpas_dmpar_max_real(domain % dminfo, maxDensity, maxDensity_global)
        config_maxMeshDensity = maxDensity_global
      END IF
      !
      ! Initialize core
      !
      timeStep = mpas_get_clock_timestep(domain % clock, ierr=err_tmp)
      CALL mpas_get_timeInterval(timeStep, dt=dt)

      block => domain % blocklist
      DO WHILE (ASSOCIATED(block))
      
         CALL ocn_init_routines_block(block, dt, ierr)

         IF (ierr==1) THEN
             CALL mpas_log_write('An error was encountered in ocn_init_routines_block', &
                                 MPAS_LOG_CRIT)
         END IF

         CALL mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         CALL mpas_pool_get_array(diagnosticsPool, 'xtime', xtime)
         xtime = startTimeStamp

         ! Set simulationStartTime only if that variable is not read from the restart
         ! file.
         CALL mpas_pool_get_array(diagnosticsPool, 'simulationStartTime', &
                                  simulationStartTime)
         IF (TRIM(simulationStartTime)=="no_date_available") THEN
             simulationStartTime = startTimeStamp
         END IF

         ! compute time since start of simulation, in days
         CALL mpas_pool_get_array(diagnosticsPool, &
                                  'daysSinceStartOfSim',daysSinceStartOfSim)
         CALL mpas_set_time(xtime_timeType, dateTimeString=xtime)
         CALL mpas_set_time(simulationStartTime_timeType, &
                            dateTimeString=simulationStartTime)
         CALL mpas_get_timeInterval(xtime_timeType - simulationStartTime_timeType, &
                                    dt=daysSinceStartOfSim)
         daysSinceStartOfSim = daysSinceStartOfSim*days_per_second

         block => block % next
         
      END DO

      IF (config_conduct_tests) THEN
         CALL mpas_timer_start("test suite")
         CALL ocn_test_suite(domain,ierr)
         CALL mpas_timer_stop("test suite")
      END IF

      CALL ocn_analysis_init(domain, err_tmp)
      ierr = ior(ierr, err_tmp)
      IF (ierr==1) THEN
          CALL mpas_log_write('An error was encountered while initializing ' &
                              // 'the analysis members in the MPAS-Ocean forward mode', &
                              MPAS_LOG_CRIT)
      END IF

      IF (config_use_surface_salinity_monthly_restoring) THEN
      
         ! Initialize the alarm for reading salinity data for restoring
         IF (config_surface_salinity_monthly_restoring_compute_interval == 'dt') THEN
            alarmTimeStep = mpas_get_clock_timestep(domain % clock, err_tmp)
            CALL mpas_get_timeInterval &
            (alarmTimeStep, timeString = &
             config_surface_salinity_monthly_restoring_compute_interval, ierr=err_tmp)
         END IF
         alarmTime = mpas_get_clock_time(domain % clock, MPAS_START_TIME, ierr=err_tmp)
   
         CALL mpas_set_timeInterval &
         (alarmTimeStep, timeString= &
          config_surface_salinity_monthly_restoring_compute_interval, ierr=err_tmp)
   
         CALL mpas_add_clock_alarm &
         (domain % clock, 'salinityDataReadAlarm', alarmTime, &
          alarmTimeInterval=alarmTimeStep, ierr=err_tmp)
         CALL mpas_reset_clock_alarm(domain % clock, 'salinityDataReadAlarm', &
                                     ierr=err_tmp)

      END IF

   END FUNCTION ocn_forward_mode_init
   !------------------------------------------------------------------------------------!
   !                                                                                    !
   !------------------------------------------------------------------------------------!
   !
   !  FUNCTION ocn_forward_mode_setup_clock
   !
   !> \brief   Initialize timer variables
   !> \author  Doug Jacobsen, Mark Petersen, Todd Ringler
   !> \date    September 2011
   !> \details
   !> This routine initializes all timer variables
   !
   FUNCTION ocn_forward_mode_setup_clock(core_clock, configs) RESULT(ierr)

      IMPLICIT NONE

      TYPE(MPAS_Clock_type), INTENT(INOUT) :: core_clock
      TYPE(mpas_pool_type), INTENT(INOUT) :: configs
      INTEGER :: ierr
      
      ! Local variables
      CHARACTER(LEN=StrKIND) :: restartTimeStamp
      CHARACTER(LEN=StrKIND), POINTER :: config_start_time, config_stop_time, &
                                         config_run_duration, config_dt, &
                                         config_restart_timestamp_name
      INTEGER :: err_tmp
      TYPE(MPAS_Time_Type) :: startTime, stopTime, alarmStartTime
      TYPE(MPAS_TimeInterval_type) :: runDuration, timeStep, alarmTimeStep

      ierr = 0

      CALL mpas_pool_get_config(configs, 'config_dt', config_dt)
      CALL mpas_pool_get_config(configs, 'config_start_time', config_start_time)
      CALL mpas_pool_get_config(configs, 'config_stop_time', config_stop_time)
      CALL mpas_pool_get_config(configs, 'config_run_duration', config_run_duration)
      CALL mpas_pool_get_config(configs, 'config_restart_timestamp_name', &
                                config_restart_timestamp_name)

      IF (TRIM(config_start_time) == "file") THEN
         OPEN (UNIT=22, FILE = config_restart_timestamp_name, FORM='FORMATTED', &
               STATUS='OLD')
         READ (22,*) restartTimeStamp
         CLOSE (UNIT=22)
         CALL mpas_set_time(curr_time=startTime, dateTimeString=restartTimeStamp, &                  
                            ierr=ierr)
      ELSE
         CALL mpas_set_time(curr_time=startTime, dateTimeString=config_start_time, &
                            ierr=err_tmp)
      END IF

      CALL mpas_set_timeInterval(timeStep, timeString=config_dt, ierr=err_tmp)
      IF (TRIM(config_run_duration) /= "none") THEN
         CALL mpas_set_timeInterval(runDuration, timeString=config_run_duration, &
                                    ierr=err_tmp)
         CALL mpas_create_clock(core_clock, startTime=startTime, timeStep=timeStep, &
                                runDuration=runDuration, ierr=err_tmp)

         IF (TRIM(config_stop_time) /= "none") THEN
            CALL mpas_set_time(curr_time=stopTime, dateTimeString=config_stop_time, &
                               ierr=err_tmp)
            IF (startTime + runduration /= stopTime) THEN
               CALL mpas_log_write('Warning: config_run_duration and config_stop_time are inconsitent: using config_run_duration.')
            END IF
         END IF
      ELSE IF (TRIM(config_stop_time) /= "none") THEN
         CALL mpas_set_time(curr_time=stopTime, dateTimeString=config_stop_time, &
                            ierr=err_tmp)
         CALL mpas_create_clock(core_clock, startTime=startTime, timeStep=timeStep, &
                                stopTime=stopTime, ierr=err_tmp)
      ELSE
         CALL mpas_log_write('Error: Neither config_run_duration nor config_stop_time were specified.')
          ierr = 1
      END IF

   END FUNCTION ocn_forward_mode_setup_clock
   !------------------------------------------------------------------------------------!
   !                                                                                    !
   !------------------------------------------------------------------------------------!
   !
   !  FUNCTION ocn_forward_mode_run
   !
   !> \brief   Main driver for MPAS-Ocean time-stepping
   !> \author  Doug Jacobsen, Mark Petersen, Todd Ringler
   !> \date    September 2011
   !> \details
   !> This function includes the time-stepping loop, and calls routines to write output
   !> and restart files.
   !
   FUNCTION ocn_forward_mode_run(domain) RESULT(ierr)
   
      IMPLICIT NONE 

      TYPE(domain_type), INTENT(INOUT) :: domain
      INTEGER :: ierr
      
      ! Local variables
      CHARACTER(LEN=StrKIND) :: timeStamp
      CHARACTER(LEN=StrKIND), POINTER :: config_restart_timestamp_name, &
                                         config_sw_absorption_type
      INTEGER :: itimestep, err, err_tmp
      LOGICAL, POINTER :: config_write_output_on_startup, config_use_ecosysTracers, &                  
                          config_use_activeTracers_surface_restoring, &
                          config_use_surface_salinity_monthly_restoring
      REAL(KIND=RKIND) :: dt
      TYPE(block_type), POINTER :: block_ptr
      TYPE(MPAS_Time_Type) :: currTime
      TYPE(mpas_pool_type), POINTER :: averagePool, meshPool, statePool, forcingPool, &
                                       diagnosticsPool, scratchPool
      TYPE(MPAS_timeInterval_type) :: timeStep

      ierr = 0

      CALL mpas_pool_get_config(domain % configs, 'config_write_output_on_startup', &
                                config_write_output_on_startup)
      CALL mpas_pool_get_config(domain % configs, 'config_restart_timestamp_name', &
                                config_restart_timestamp_name)
      CALL mpas_pool_get_config(domain % configs, 'config_sw_absorption_type', &
                                config_sw_absorption_type)
      CALL mpas_pool_get_config(domain % configs, 'config_use_ecosysTracers', &
                                config_use_ecosysTracers)
      CALL mpas_pool_get_config(domain % configs, &
                                'config_use_activeTracers_surface_restoring',  &
                                config_use_activeTracers_surface_restoring)
      CALL mpas_pool_get_config(domain % configs, &
                                'config_use_surface_salinity_monthly_restoring',  &
                                config_use_surface_salinity_monthly_restoring)

      ! Eventually, dt should be domain specific
      timeStep = mpas_get_clock_timestep(domain % clock, ierr=ierr)
      CALL mpas_get_timeInterval(timeStep, dt=dt)

      currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr)
      CALL mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)
#ifdef MPAS_DEBUG
      CALL mpas_log_write('Initial time ' // TRIM(timeStamp))
#endif

      CALL ocn_analysis_compute_startup(domain, err)

      IF (config_write_output_on_startup) THEN
          CALL mpas_timer_start('io_write')
          CALL mpas_stream_mgr_write(domain % streamManager, 'output', &
                                     forceWriteNow=.TRUE., ierr=ierr)
          CALL mpas_timer_stop('io_write')
      END IF

      ! read initial data required for variable shortwave
      CALL mpas_timer_start('io_shortwave',.FALSE.)
      CALL ocn_get_shortWaveData(domain % streamManager, domain, domain % clock, .TRUE.)
      CALL mpas_timer_stop('io_shortwave')

      ! read initial data required for ecosys forcing
      IF (config_use_ecosysTracers) THEN
        CALL mpas_timer_start('io_ecosys',.FALSE.)
        CALL ocn_get_ecosysData(domain % streamManager, domain, domain % clock, .TRUE.)
        CALL mpas_timer_stop('io_ecosys')
      END IF

      ! read initial data required for monthly surface salinity restoring
      ! always execute this call as initial data needed regardless of alarm
      IF (config_use_activeTracers_surface_restoring .AND. &
          config_use_surface_salinity_monthly_restoring) THEN
        CALL mpas_timer_start('io_monthly_surface_salinity', .FALSE.)
        CALL ocn_get_surfaceSalinityData(domain % streamManager, domain, &
                                         domain % clock, .TRUE.)
        CALL mpas_timer_stop('io_monthly_surface_salinity')
      END IF

      ! During integration, time level 1 stores the model state at the beginning of the
      ! time step, and time level 2 stores the state advanced dt in time by timestep(...)
      itimestep = 0

      DO WHILE (.not. mpas_is_clock_stop_time(domain % clock))
         CALL mpas_timer_start('io_read')
         CALL mpas_stream_mgr_read(domain % streamManager, ierr=ierr)
         CALL mpas_timer_stop('io_read')
         CALL mpas_timer_start('reset_io_alarms')
         CALL mpas_stream_mgr_reset_alarms(domain % streamManager, &
                                           direction=MPAS_STREAM_INPUT, ierr=ierr)

         ! Also restart all block_* streams, if any are defined.
         IF (mpas_stream_mgr_stream_exists(domain % streamManager, 'block_.*')) THEN
            CALL mpas_stream_mgr_reset_alarms(domain % streamManager, &
                                              streamID='block_.*', ierr=ierr)
         END IF
         CALL mpas_timer_stop('reset_io_alarms')

         itimestep = itimestep + 1
         CALL mpas_advance_clock(domain % clock)

         currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr)
         CALL mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)
         CALL mpas_log_write('Doing timestep ' // TRIM(timeStamp))

#ifdef MPAS_DEBUG
         CALL mpas_log_write( '   Computing surface flux arrays')
#endif
         block_ptr => domain % blocklist
         DO WHILE(ASSOCIATED(block_ptr))
           CALL mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
           CALL mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
           CALL mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
           CALL mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', &
                                      diagnosticsPool)
           CALL mpas_pool_get_subpool(block_ptr % structs, 'scratch', scratchPool)
           CALL ocn_forcing_build_fraction_absorbed_array(meshPool, statePool, &
                                                          diagnosticsPool, forcingPool, &
                                                          ierr, 1)
           CALL mpas_timer_start("land_ice_build_arrays")
           CALL ocn_surface_land_ice_fluxes_build_arrays(meshPool, diagnosticsPool, &
                                                         forcingPool, scratchPool, &
                                                         statePool, dt, err)
           CALL mpas_timer_stop("land_ice_build_arrays")

           CALL ocn_frazil_forcing_build_arrays(domain, meshPool, forcingPool, &
                                                diagnosticsPool, statePool, err)
           block_ptr => block_ptr % next
         END DO

         CALL mpas_timer_start("time integration")

#ifdef MPAS_DEBUG
         CALL mpas_log_write( '   Computing forward time step')
#endif
         !$omp parallel default(firstprivate) shared(domain, dt, timeStamp)

         CALL ocn_timestep(domain, dt, timeStamp)

         !$omp end parallel

         CALL mpas_timer_stop("time integration")

         ! Move time level 2 fields back into time level 1 for next time step
#ifdef MPAS_DEBUG
         CALL mpas_log_write( '   Shifting time levels')
#endif
         CALL mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)
         CALL mpas_pool_shift_time_levels(statePool)

#ifdef MPAS_DEBUG
         CALL mpas_log_write( '   Handling analysis members')
#endif
         CALL ocn_analysis_compute(domain, err)
         CALL ocn_analysis_restart(domain, err)
         CALL ocn_analysis_write(domain, err)

#ifdef MPAS_DEBUG
         CALL mpas_log_write( '   Performing I/O')
#endif
         CALL mpas_timer_start('io_write')
         CALL mpas_stream_mgr_write(domain % streamManager, streamID='output', ierr=ierr)
         CALL mpas_timer_stop('io_write')
         CALL mpas_timer_start('reset_io_alarms')
         CALL mpas_stream_mgr_reset_alarms(domain % streamManager, streamID='output', &
                                           ierr=ierr)
         CALL mpas_timer_stop('reset_io_alarms')

         IF (mpas_stream_mgr_ringing_alarms(domain % streamManager, streamID='restart', &
                                            direction=MPAS_STREAM_OUTPUT, ierr=ierr)) &
         THEN
#ifdef MPAS_DEBUG
            CALL mpas_log_write( '   Writing restart timestamp file')
#endif
            IF (domain % dminfo % my_proc_id == 0) THEN
               OPEN(UNIT=22, FILE=config_restart_timestamp_name, FORM='FORMATTED', &
                    STATUS='REPLACE')
               WRITE (22, *) TRIM(timeStamp)
               CLOSE (UNIT=22)
            END IF

            IF (TRIM(config_sw_absorption_type)=='ohlmann00') &
            CALL ocn_shortwave_forcing_write_restart(domain)

            IF (config_use_ecosysTracers) CALL ocn_ecosys_forcing_write_restart(domain)

            IF (config_use_activeTracers_surface_restoring .AND. &
                config_use_surface_salinity_monthly_restoring)  &
               CALL ocn_salinity_restoring_forcing_write_restart(domain)

         END IF

         CALL mpas_timer_start('io_write')
         CALL mpas_stream_mgr_write(domain % streamManager, streamID='restart', &
                                    ierr=ierr)
         CALL mpas_timer_stop('io_write')


         CALL mpas_timer_start('reset_io_alarms')
         CALL mpas_stream_mgr_reset_alarms(domain % streamManager, streamID='restart', &
                                           ierr=ierr)
         CALL mpas_timer_stop('reset_io_alarms')

         CALL mpas_timer_start('io_write')
         CALL mpas_stream_mgr_write(domain % streamManager, ierr=ierr)
         CALL mpas_timer_stop('io_write')
         CALL mpas_timer_start('reset_io_alarms')
         CALL mpas_stream_mgr_reset_alarms(domain % streamManager, &
                                           direction=MPAS_STREAM_OUTPUT, ierr=ierr)
         CALL mpas_timer_stop('reset_io_alarms')

         ! read next time level data required for variable shortwave
         CALL mpas_timer_start('io_shortwave',.FALSE.)
         CALL ocn_get_shortWaveData(domain % streamManager, domain, domain % clock, &
                                    .FALSE.)
         CALL mpas_timer_stop('io_shortwave')

         ! read next time level data required for ecosys forcing
         IF (config_use_ecosysTracers) THEN
           CALL mpas_timer_start('io_ecosys',.FALSE.)
           CALL ocn_get_ecosysData(domain % streamManager, domain, domain % clock, &
                                   .FALSE.)
           CALL mpas_timer_stop('io_ecosys')
         END IF

         ! read next time level data required for monthly surface salinity restoring
         IF (config_use_activeTracers_surface_restoring .AND. &
             config_use_surface_salinity_monthly_restoring) THEN
           IF (mpas_is_alarm_ringing(domain % clock, 'salinityDataReadAlarm', &
               ierr=err_tmp)) THEN
              CALL mpas_reset_clock_alarm(domain % clock, 'salinityDataReadAlarm', &
                                          ierr=err_tmp)
              CALL mpas_timer_start('io_monthly_surface_salinity',.FALSE.)
              CALL ocn_get_surfaceSalinityData(domain % streamManager, domain, &
                                               domain % clock, .FALSE.)
              CALL mpas_timer_stop('io_monthly_surface_salinity')
           END IF
         END IF

         ! Validate that the state is OK to run with for the next timestep.
         CALL ocn_validate_state(domain, timeLevel=1)

      END DO

   END FUNCTION ocn_forward_mode_run
   !------------------------------------------------------------------------------------!
   !                                                                                    !
   !------------------------------------------------------------------------------------!
   !
   !  FUNCTION ocn_forward_mode_finalize
   !
   !> \brief   Finalize MPAS-Ocean Forward Mode
   !> \author  Doug Jacobsen, Mark Petersen, Todd Ringler
   !> \date    September 2011
   !> \details
   !> This function finalizes the MPAS-Ocean core in forward mode.
   !
   FUNCTION ocn_forward_mode_finalize(domain) RESULT(iErr)
      
      IMPLICIT NONE 
      
      TYPE(domain_type), INTENT(INOUT) :: domain
      INTEGER :: ierr

      CALL ocn_analysis_finalize(domain, ierr)
      CALL mpas_destroy_clock(domain % clock, ierr)
      CALL mpas_decomp_destroy_decomp_list(domain % decompositions)

   END FUNCTION ocn_forward_mode_finalize

END MODULE ocn_forward_mode
