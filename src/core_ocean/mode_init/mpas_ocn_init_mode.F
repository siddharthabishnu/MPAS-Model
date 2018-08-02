!
! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!  ocn_init_mode
!
!> \brief Main driver for MPAS ocean core
!> \author Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date   September 2011
!> \details
!> This module contains initialization and timestep drivers for the MPAS ocean core.
!
MODULE ocn_init_mode

   USE mpas_kind_types
   USE mpas_derived_types
   USE mpas_pool_routines
   USE mpas_stream_manager
   USE mpas_timekeeping
   USE mpas_dmpar
   USE mpas_timer
   USE mpas_io_units
   USE mpas_constants
   USE mpas_decomp
   USE ocn_init_routines
   USE ocn_equation_of_state
   USE ocn_constants
   USE ocn_init_spherical_utils
!  USE ocn_init_TEMPLATE
   USE ocn_init_baroclinic_channel
   USE ocn_init_lock_exchange
   USE ocn_init_internal_waves
   USE ocn_init_overflow
   USE ocn_init_global_ocean
   USE ocn_init_cvmix_WSwSBF
   USE ocn_init_iso
   USE ocn_init_soma
   USE ocn_init_ziso
   USE ocn_init_sub_ice_shelf_2D
   USE ocn_init_periodic_planar
   USE ocn_init_ecosys_column
   USE ocn_init_sea_mount
   USE ocn_init_isomip
   USE ocn_init_isomip_plus
   USE ocn_init_surface_gravity_wave

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: ocn_init_mode_init, ocn_init_mode_run, ocn_init_mode_finalize
   PUBLIC :: ocn_init_mode_setup_clock, ocn_init_mode_validate_configuration

   CONTAINS
!
!  FUNCTION ocn_init_mode_init
!
!> \brief   Initialize MPAS-Ocean core in init mode
!> \author  Doug Jacobsen
!> \date    06/15/2015
!> \details
!> This function calls all initializations required to start MPAS-Ocean in init mode.
!
   FUNCTION ocn_init_mode_init(domain, startTimeStamp) RESULT(iErr)

      TYPE(domain_type), INTENT(INOUT) :: domain
      CHARACTER(LEN=*), INTENT(OUT) :: startTimeStamp
      INTEGER :: iErr
      
      ! Local variables
      CHARACTER(LEN=StrKIND), POINTER :: xtime
      CHARACTER(LEN=StrKIND), POINTER :: config_vert_coord_movement, &
                                         config_pressure_gradient_type      
      INTEGER :: err_tmp
      INTEGER, POINTER :: nVertLevels      
      LOGICAL, POINTER :: config_do_restart, config_filter_btr_mode, &
                          config_conduct_tests, config_write_stats_on_startup            
      REAL(KIND=RKIND) :: dt, maxDensity, maxDensity_global
      REAL(KIND=RKIND), POINTER :: config_maxMeshDensity
      REAL(KIND=RKIND), DIMENSION(:), POINTER :: meshDensity
      TYPE(block_type), POINTER :: block
      TYPE(mpas_pool_type), POINTER :: meshPool, diagnosticsPool 
      TYPE(MPAS_Time_Type) :: startTime
      TYPE(MPAS_TimeInterval_type) :: timeStep

      iErr = 0
      !
      ! Set startTimeStamp based on the start time of the simulation clock
      !
      startTime = mpas_get_clock_time(domain % clock, MPAS_START_TIME, err_tmp)
      CALL mpas_get_time(startTime, dateTimeString=startTimeStamp)
      iErr = ior(iErr, err_tmp)

      ! Setup ocean config pool
      CALL ocn_constants_init(domain % configs, domain % packages)
      !
      ! Read input data for model
      !
      CALL mpas_timer_start('io_read', .false.)
      CALL mpas_stream_mgr_read(domain % streamManager, streamID='input_init', &
                                iErr=err_tmp)
      CALL mpas_timer_stop('io_read')

      CALL mpas_timer_start('reset_io_alarms', .false.)
      CALL mpas_stream_mgr_reset_alarms(domain % streamManager, streamID='input_init', &
                                        iErr=err_tmp)
!     CALL mpas_stream_mgr_reset_alarms(domain % streamManager, &
!                                       direction=MPAS_STREAM_OUTPUT, iErr=err_tmp)
      CALL mpas_timer_stop('reset_io_alarms')

      ! Read the remaining input streams
      CALL mpas_timer_start('io_read', .false.)
      CALL mpas_stream_mgr_read(domain % streamManager, iErr=err_tmp)
      iErr = ior(iErr, err_tmp)
      CALL mpas_timer_stop('io_read')
      CALL mpas_timer_start('reset_io_alarms', .false.)
      CALL mpas_stream_mgr_reset_alarms(domain % streamManager, &
                                        direction=MPAS_STREAM_INPUT, iErr=err_tmp)
      iErr = ior(iErr, err_tmp)
      CALL mpas_timer_stop('reset_io_alarms')

      ! Initialize submodules before initializing blocks.
      CALL ocn_equation_of_state_init(err_tmp)
      iErr = ior(iErr, err_tmp)
      IF (iErr.eq.1) THEN
          CALL mpas_log_write('An error was encountered while initializing the &
                               &MPAS-Ocean init mode', MPAS_LOG_CRIT)
      END IF
      !
      ! Initialize core
      !
      timeStep = mpas_get_clock_timestep(domain % clock, iErr=err_tmp)
      CALL mpas_get_timeInterval(timeStep, dt=dt)

      block => domain % blocklist
      DO WHILE (ASSOCIATED(block))
         CALL mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         CALL mpas_pool_get_array(diagnosticsPool, 'xtime', xtime)
         xtime = startTimeStamp
         block => block % next
      END DO

      ! Expand sphere if it needs to be expanded
      CALL ocn_init_expand_sphere(domain, domain % streamManager, a, iErr)

   END FUNCTION ocn_init_mode_init
   !------------------------------------------------------------------------------------!
   !                                                                                    !
   !------------------------------------------------------------------------------------!
!
!  FUNCTION ocn_init_mode_setup_clock
!
!> \brief   Setup MPAS-Ocean clock
!> \author  Doug Jacobsen
!> \date    06/15/2015
!> \details
!>  This function initializes the MPAS-Ocean clock for the init mode.
!
   FUNCTION ocn_init_mode_setup_clock(core_clock, configs) RESULT(iErr)

      IMPLICIT NONE

      TYPE(MPAS_Clock_type), INTENT(INOUT) :: core_clock
      TYPE(mpas_pool_type), INTENT(INOUT) :: configs
      INTEGER :: iErr
      
      ! Local variables
      TYPE(MPAS_Time_Type) :: startTime, stopTime, alarmStartTime
      TYPE(MPAS_TimeInterval_type) :: runDuration, timeStep, alarmTimeStep
      CHARACTER(LEN=StrKIND) :: restartTimeStamp
      CHARACTER(LEN=StrKIND), POINTER :: config_start_time, config_stop_time, &
                                         config_run_duration, config_dt, &
                                         config_restart_timestamp_name
      INTEGER :: err_tmp

      iErr = 0

      CALL mpas_pool_get_config(configs, 'config_dt', config_dt)
      CALL mpas_pool_get_config(configs, 'config_start_time', config_start_time)
      CALL mpas_pool_get_config(configs, 'config_stop_time', config_stop_time)
      CALL mpas_pool_get_config(configs, 'config_run_duration', config_run_duration)
      CALL mpas_pool_get_config(configs, 'config_restart_timestamp_name', &
                                config_restart_timestamp_name)

      CALL mpas_set_time(startTime, dateTimeString=config_start_time, iErr=err_tmp)
      CALL mpas_set_timeInterval(timeStep, timeString=config_dt, iErr=err_tmp)
      IF (TRIM(config_run_duration) /= "none") THEN
         CALL mpas_set_timeInterval(runDuration, timeString=config_run_duration, &
                                    iErr=err_tmp)
         CALL mpas_create_clock(core_clock, startTime=startTime, timeStep=timeStep, &
                                runDuration=runDuration, iErr=err_tmp)

         IF (TRIM(config_stop_time) /= "none") THEN
            CALL mpas_set_time(curr_time=stopTime, dateTimeString=config_stop_time, &
                               iErr=err_tmp)
            IF (startTime + runduration /= stopTime) THEN
               CALL mpas_log_write('Warning: config_run_duration and config_stop_time &
                                    &are inconsitent: using config_run_duration.')
            END IF
         END IF
      else IF (TRIM(config_stop_time) /= "none") THEN
         CALL mpas_set_time(curr_time=stopTime, dateTimeString=config_stop_time, &
                            iErr=err_tmp)
         CALL mpas_create_clock(core_clock, startTime=startTime, timeStep=timeStep, &
                                stopTime=stopTime, iErr=err_tmp)
      else
         CALL mpas_log_write(' Warning: config_run_duration and config_start_time were &
                              &"none", setting run duration to 1 second.')
         CALL mpas_set_timeInterval(runDuration, timeString="0000_00:00:01", &
                                    iErr=err_tmp)
         CALL mpas_create_clock(core_clock, startTime=startTime, timeStep=timeStep, &
                                runDuration=runDuration, iErr=err_tmp)
      END IF

   END FUNCTION ocn_init_mode_setup_clock
   !------------------------------------------------------------------------------------!
   !                                                                                    !
   !------------------------------------------------------------------------------------!
!
!  FUNCTION ocn_init_mode_run
!
!> \brief   MPAS-Ocean init mode run step
!> \author  Doug Jacobsen
!> \date    06/15/2015
!> \details
!>  This function sets up the initial configuration using the MPAS-Ocean init
!>  mode.
!
   FUNCTION ocn_init_mode_run(domain) RESULT(iErr)

      TYPE(domain_type), INTENT(INOUT) :: domain
      INTEGER :: iErr
      
      ! Local variables
      CHARACTER(LEN=StrKIND) :: timeStamp
      CHARACTER(LEN=StrKIND), POINTER :: config_init_configuration
      INTEGER :: itimestep
      REAL(KIND=RKIND) :: dt
      TYPE(block_type), POINTER :: block_ptr
      TYPE(MPAS_Time_Type) :: currTime
      TYPE(mpas_pool_type), POINTER :: averagePool, meshPool, statePool, forcingPool
      TYPE(MPAS_timeInterval_type) :: timeStep

      iErr = 0

      ! Eventually, dt should be domain specific
      timeStep = mpas_get_clock_timestep(domain % clock, iErr=iErr)
      CALL mpas_get_timeInterval(timeStep, dt=dt)

      currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, iErr)
      CALL mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, iErr=iErr)

      CALL mpas_pool_get_config(domain % configs, 'config_init_configuration', &
                                config_init_configuration)
      CALL mpas_log_write(' Generating configuration: ' // &
                          TRIM(config_init_configuration))

      CALL ocn_init_setup_baroclinic_channel(domain, iErr)
      CALL ocn_init_setup_lock_exchange(domain, iErr)
      CALL ocn_init_setup_internal_waves(domain, iErr)
      CALL ocn_init_setup_overflow(domain, iErr)
      CALL ocn_init_setup_global_ocean(domain, iErr)
      CALL ocn_init_setup_cvmix_WSwSBF(domain, iErr)
      CALL ocn_init_setup_iso(domain, iErr)
      CALL ocn_init_setup_soma(domain, iErr)
      CALL ocn_init_setup_ziso(domain, iErr)
      CALL ocn_init_setup_sub_ice_shelf_2D(domain, iErr)
      CALL ocn_init_setup_periodic_planar(domain, iErr)
      CALL ocn_init_setup_ecosys_column(domain, iErr)
      CALL ocn_init_setup_sea_mount(domain, iErr)
      CALL ocn_init_setup_isomip(domain, iErr)
      CALL ocn_init_setup_isomip_plus(domain, iErr)
      CALL ocn_init_setup_surface_gravity_wave(domain, iErr)
!     CALL ocn_init_setup_TEMPLATE(domain, iErr)

      CALL mpas_log_write(' Completed setup of: ' // TRIM(config_init_configuration))
      CALL mpas_timer_start('io_write', .false.)
      CALL mpas_stream_mgr_write(domain % streamManager, iErr=iErr)
      CALL mpas_timer_stop('io_write')
      CALL mpas_timer_start('reset_io_alarms', .false.)
      CALL mpas_stream_mgr_reset_alarms(domain % streamManager, &
                                        direction=MPAS_STREAM_OUTPUT, iErr=iErr)
      CALL mpas_timer_stop('reset_io_alarms')

   END FUNCTION ocn_init_mode_run
   !------------------------------------------------------------------------------------!
   !                                                                                    !
   !------------------------------------------------------------------------------------!
!
!  function ocn_init_mode_finalize
!
!> \brief   MPAS-Ocean init mode run step
!> \author  Doug Jacobsen
!> \date    06/15/2015
!> \details
!>  This function sets up the initial configuration using the MPAS-Ocean init
!>  mode.
!
   FUNCTION ocn_init_mode_finalize(domain) RESULT(iErr)

      TYPE(domain_type), INTENT(INOUT) :: domain
      INTEGER :: iErr

      iErr = 0

      CALL mpas_destroy_clock(domain % clock, iErr)

      CALL mpas_decomp_destroy_decomp_list(domain % decompositions)

   END FUNCTION ocn_init_mode_finalize
   !------------------------------------------------------------------------------------!
   !                                                                                    !
   !------------------------------------------------------------------------------------!
!
!  routine ocn_init_validate_configuration
!
!> \brief   Configuration validation routine
!> \author  Doug Jacobsen
!> \date    03/20/2015
!> \details
!>  This routine is used to validate the namelist options against the
!>  configuration definition. It will call the validate routines for each of the
!>  configurations to ensure namelist options are set in a valid way.
!
   SUBROUTINE ocn_init_mode_validate_configuration(configPool, packagePool, iocontext, &
                                                   iErr)
                                                   
      IMPLICIT NONE                                                 

      TYPE(mpas_pool_type), INTENT(INOUT) :: configPool 
      !< Input: Pool with namelist options
      TYPE(mpas_pool_type), INTENT(INOUT) :: packagePool !< Input: Pool with packages
      TYPE(mpas_io_context_type), INTENT(INOUT) :: iocontext
      INTEGER, INTENT(OUT) :: iErr !< Output: Error core


      ! Local variables
      INTEGER :: err_tmp
      LOGICAL, POINTER :: cullCellsActive, config_write_cull_cell_mask
      
      iErr = 0

      CALL mpas_pool_get_config(configPool, 'config_write_cull_cell_mask', &
                                config_write_cull_cell_mask)
      CALL mpas_pool_get_package(packagePool, 'cullCellsActive', cullCellsActive)

      IF (config_write_cull_cell_mask) THEN
         cullCellsActive = .true.
      END IF

      CALL ocn_init_validate_baroclinic_channel(configPool, packagePool, iocontext, &
                                                iErr=err_tmp)
      iErr = ior(iErr, err_tmp)
      CALL ocn_init_validate_lock_exchange(configPool, packagePool, iocontext, &
                                           iErr=err_tmp)
      iErr = ior(iErr, err_tmp)
      CALL ocn_init_validate_internal_waves(configPool, packagePool, iocontext, &
                                            iErr=err_tmp)
      iErr = ior(iErr, err_tmp)
      CALL ocn_init_validate_overflow(configPool, packagePool, iocontext, iErr=err_tmp)
      iErr = ior(iErr, err_tmp)
      CALL ocn_init_validate_global_ocean(configPool, packagePool, iocontext, &
                                          iErr=err_tmp)
      iErr = ior(iErr, err_tmp)
      CALL ocn_init_validate_cvmix_WSwSBF(configPool, packagePool, iocontext, &
                                          iErr=err_tmp)
      iErr = ior(iErr, err_tmp)
      CALL ocn_init_validate_iso(configPool, packagePool, iocontext, iErr=err_tmp)
      iErr = ior(iErr, err_tmp)
      CALL ocn_init_validate_soma(configPool, packagePool, iocontext, iErr=err_tmp)
      iErr = ior(iErr, err_tmp)
      CALL ocn_init_validate_ziso(configPool, packagePool, iErr=err_tmp)
      iErr = ior(iErr, err_tmp)
      CALL ocn_init_validate_sub_ice_shelf_2D(configPool, packagePool,iocontext, &
                                              iErr=err_tmp)
      iErr = ior(iErr, err_tmp)
      CALL ocn_init_validate_periodic_planar(configPool, packagePool, iocontext, &
                                             iErr=err_tmp)
      iErr = ior(iErr, err_tmp)
      CALL ocn_init_validate_ecosys_column(configPool, packagePool, iocontext, &
                                           iErr=err_tmp)
      iErr = ior(iErr, err_tmp)
      CALL ocn_init_validate_sea_mount(configPool, packagePool, iocontext, iErr=err_tmp)
      iErr = ior(iErr, err_tmp)
      CALL ocn_init_validate_isomip(configPool, packagePool, iocontext, iErr=err_tmp)
      iErr = ior(iErr, err_tmp)
      CALL ocn_init_validate_isomip_plus(configPool, packagePool, iocontext, &
                                         iErr=err_tmp)
      iErr = ior(iErr, err_tmp)
      CALL ocn_init_validate_surface_gravity_wave(configPool, packagePool, iocontext, &
                                                  iErr=err_tmp)
      iErr = ior(iErr, err_tmp)
!     CALL ocn_init_validate_TEMPLATE(configPool, packagePool, iocontext, iErr=err_tmp)
!     iErr = ior(iErr, err_tmp)
      
   END SUBROUTINE ocn_init_mode_validate_configuration

END MODULE ocn_init_mode