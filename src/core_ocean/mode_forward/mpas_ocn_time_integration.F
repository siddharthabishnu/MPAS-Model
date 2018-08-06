!
! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!  ocn_time_integration
!
!> \brief MPAS ocean time integration driver
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   September 2011
!> \details
!>  This module contains the main driver routine for calling
!>  the time integration scheme
!
MODULE ocn_time_integration

   USE mpas_derived_types
   USE mpas_pool_routines
   USE mpas_constants
   USE mpas_timekeeping
   USE mpas_dmpar
   USE mpas_vector_reconstruction
   USE mpas_spline_interpolation
   USE mpas_timer
   USE mpas_log

   USE ocn_constants
   USE ocn_time_integration_rk4
   USE ocn_time_integration_split

   IMPLICIT NONE
   PRIVATE
   SAVE

   PUBLIC :: ocn_timestep, ocn_timestep_init

   !------------------------------------------------------------------------------------!
   !                                                                  
   ! Public parameters
   !
   !------------------------------------------------------------------------------------!

   !------------------------------------------------------------------------------------!
   !
   ! Public member functions
   !
   !------------------------------------------------------------------------------------!

   !------------------------------------------------------------------------------------!
   !
   ! Private module variables
   !
   !------------------------------------------------------------------------------------!

   LOGICAL :: rk4On, splitOn

   CONTAINS
   !
   !  ocn_timestep
   !
   !> \brief MPAS ocean time integration driver
   !> \author Mark Petersen, Doug Jacobsen, Todd Ringler
   !> \date   September 2011
   !> \details
   !> This routine handles a single timestep for the ocean. It determines the time 
   !> integrator that will be used for the run, and calls the appropriate one.
   !
      SUBROUTINE ocn_timestep(domain, dt, timeStamp)
      !
      ! Advance model state forward in time by the specified time step
      !
      ! Input: domain - current model state in time level 1 
      !                 (e.g., time_levs(1) % state % h(:,:)) plus mesh meta-data
      ! Output: domain - upon exit, time level 2 
      !                  (e.g., time_levs(2) % state % h(:,:)) contains model state 
      !                  advanced forward in time by dt seconds
      !
      IMPLICIT NONE

      TYPE(domain_type), INTENT(INOUT) :: domain
      REAL(KIND=RKIND), INTENT(IN) :: dt
      CHARACTER(LEN=*), INTENT(IN) :: timeStamp
      
      ! Local variables
      CHARACTER(LEN=StrKIND), POINTER :: xtime, simulationStartTime
      REAL(KIND=RKIND), POINTER :: daysSinceStartOfSim
      TYPE(dm_info) :: dminfo
      TYPE(block_type), POINTER :: block
      TYPE(mpas_pool_type), POINTER :: diagnosticsPool, statePool, meshPool
      TYPE(MPAS_Time_type) :: xtime_timeType, simulationStartTime_timeType

      IF (rk4On) THEN
         CALL ocn_time_integrator_rk4(domain, dt)
      ELSE IF (splitOn) THEN
         CALL ocn_time_integrator_split(domain, dt)
      END IF

      block => domain % blocklist
     
      DO WHILE (ASSOCIATED(block))
     
         CALL mpas_pool_get_subpool(block % structs, 'state', statePool)
         CALL mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         CALL mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

         CALL mpas_pool_get_array(diagnosticsPool, 'xtime', xtime)

         xtime = timeStamp

         ! compute time since start of simulation, in days
         CALL mpas_pool_get_array(diagnosticsPool, 'simulationStartTime', &
                                  simulationStartTime)
         CALL mpas_pool_get_array(diagnosticsPool, 'daysSinceStartOfSim', &
                                  daysSinceStartOfSim)
         CALL mpas_set_time(xtime_timeType, dateTimeString=xtime)
         CALL mpas_set_time(simulationStartTime_timeType, &
                            dateTimeString=simulationStartTime)
         CALL mpas_get_timeInterval(xtime_timeType - simulationStartTime_timeType, &
                                    dt=daysSinceStartOfSim)

         !$omp single
         daysSinceStartOfSim = daysSinceStartOfSim*days_per_second
         !$omp end single

         block => block % next
        
      END DO

      END SUBROUTINE ocn_timestep
      !---------------------------------------------------------------------------------!
      !                                                                                 !
      !---------------------------------------------------------------------------------!
      SUBROUTINE ocn_timestep_init(err)
      
      IMPLICIT NONE

      INTEGER, INTENT(OUT) :: err
      
      ! Local variables
      CHARACTER(LEN=StrKIND), POINTER :: config_time_integrator

      err = 0

      CALL mpas_pool_get_config(ocnConfigs, 'config_time_integrator', &
                                config_time_integrator)

      rk4On = .FALSE.
      splitOn = .FALSE.

      IF (TRIM(config_time_integrator) == 'RK4') THEN
          rk4On = .TRUE.
      ELSE IF (TRIM(config_time_integrator) == 'split_explicit' &
               .OR. TRIM(config_time_integrator) == 'unsplit_explicit') THEN
          splitOn = .TRUE.
      ELSE
          err = 1
          CALL mpas_log_write('Incorrect choice for config_time_integrator:' &
                              // TRIM(config_time_integrator) // &
                              '   choices are: RK4, split_explicit, unsplit_explicit', &
                              MPAS_LOG_CRIT)
      END IF

      END SUBROUTINE ocn_timestep_init

END MODULE ocn_time_integration