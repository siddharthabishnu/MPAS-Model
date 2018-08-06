! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!  ocn_time_integration_rk4
!
!> \brief MPAS ocean RK4 Time integration scheme
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   September 2011
!> \details
!> This module contains the RK4 time integration routine.
!
MODULE ocn_time_integration_rk4

   USE mpas_derived_types
   USE mpas_pool_routines
   USE mpas_constants
   USE mpas_dmpar
   USE mpas_threading
   USE mpas_vector_reconstruction
   USE mpas_spline_interpolation
   USE mpas_timer

   USE ocn_constants
   USE ocn_tendency
   USE ocn_diagnostics
   USE ocn_gm

   USE ocn_equation_of_state
   USE ocn_vmix
   USE ocn_time_average_coupled

   USE ocn_effective_density_in_land_ice

   IMPLICIT NONE
   PRIVATE
   SAVE

   !-----------------------------------------------------------------------------------!
   !
   ! Public parameters
   !
   !-----------------------------------------------------------------------------------!

   !-----------------------------------------------------------------------------------!
   !
   ! Public member functions
   !
   !-----------------------------------------------------------------------------------!

   PUBLIC :: ocn_time_integrator_rk4

   CONTAINS
   !
   !  ocn_time_integrator_rk4
   !
   !> \brief MPAS ocean RK4 Time integration scheme
   !> \author Mark Petersen, Doug Jacobsen, Todd Ringler
   !> \date   September 2011
   !> \details
   !> This routine integrates one timestep (dt) using an RK4 time integrator.
   !
      SUBROUTINE ocn_time_integrator_rk4(domain, dt)
      !
      ! Advance model state forward in time by the specified time step using
      ! 4th order Runge-Kutta
      !
      ! Input: domain - current model state in time level 1 
      !                 (e.g., time_levs(1) % state % h(:,:)) plus mesh meta-data
      ! Output: domain - upon exit, time level 2 (e.g., time_levs(2) % state % h(:,:)) 
      !                  contains model state advanced forward in time by dt seconds
      !
      IMPLICIT NONE

      TYPE(domain_type), INTENT(INOUT) :: domain !< Input/Output: domain information
      REAL(KIND=RKIND), INTENT(IN) :: dt !< Input: timestep
      
      ! Local variables
      INTEGER :: iCell, iEdge, k, i, err, rk_step
      REAL(KIND=RKIND) :: coef
      REAL(KIND=RKIND), DIMENSION(:,:), POINTER :: vertViscTopOfEdge, vertDiffTopOfCell
      REAL(KIND=RKIND), DIMENSION(:), ALLOCATABLE :: A,C,uTemp
      REAL(KIND=RKIND), DIMENSION(:,:), ALLOCATABLE :: tracersTemp
      REAL(KIND=RKIND), DIMENSION(1:4) :: rk_weights, rk_substep_weights
      TYPE(block_type), POINTER :: block
      TYPE(mpas_pool_type), POINTER :: tendPool, tracersTendPool, statePool, &
                                       tracersPool, meshPool, provisStatePool, &
                                       provisTracersPool, diagnosticsPool, &
                                       verticalMeshPool, forcingPool, scratchPool, &
                                       swForcingPool, nextProvisPool, prevProvisPool
   
      ! Dimensions
      INTEGER, POINTER :: nCells, nEdges, nVertLevels, num_tracers

      ! Config options
      CHARACTER(LEN=StrKIND), POINTER :: config_land_ice_flux_mode
      LOGICAL, POINTER :: config_prescribe_velocity, config_prescribe_thickness, &
                          config_filter_btr_mode, config_use_freq_filtered_thickness, &
                          config_use_standardGM, config_use_cvmix_kpp, &
                          config_use_tracerGroup
      REAL(KIND=RKIND), POINTER :: config_mom_del4

      ! State indices
      INTEGER, POINTER :: indexTemperature, indexSalinity

      ! Diagnostics Indices
      INTEGER, POINTER :: indexSurfaceVelocityZonal, indexSurfaceVelocityMeridional, &
                          indexSSHGradientZonal, indexSSHGradientMeridional

      ! Mesh array pointers
      INTEGER, DIMENSION(:), POINTER :: maxLevelCell, maxLevelEdgeTop

      ! Provis Array Pointers
      REAL(KIND=RKIND), DIMENSION(:,:), POINTER &
      :: normalVelocityProvis, layerThicknessProvis, highFreqThicknessProvis, &
         lowFreqDivergenceProvis
      REAL(KIND=RKIND), DIMENSION(:,:,:), POINTER :: tracersGroupProvis

      ! Tend Array Pointers
      REAL(KIND=RKIND), DIMENSION(:,:), POINTER &
      :: highFreqThicknessTend, lowFreqDivergenceTend, normalVelocityTend, &
         layerThicknessTend
      REAL(KIND=RKIND), DIMENSION(:,:,:), POINTER :: tracersGroupTend

      ! Diagnostics Array Pointers
      REAL(KIND=RKIND), DIMENSION(:), POINTER :: gradSSH, gradSSHX, gradSSHY, &
                                                 gradSSHZ, gradSSHZonal, &
                                                 gradSSHMeridional
      REAL(KIND=RKIND), DIMENSION(:,:), POINTER &
      :: layerThicknessEdge, vertAleTransportTop, normalTransportVelocity, &
         normalGMBolusVelocity, velocityX, velocityY, velocityZ, velocityZonal, &
         velocityMeridional, surfaceVelocity, sshGradient

      ! State Array Pointers
      REAL(KIND=RKIND), DIMENSION(:,:), POINTER :: normalVelocityCur, normalVelocityNew
      REAL(KIND=RKIND), DIMENSION(:,:), POINTER :: layerThicknessCur, layerThicknessNew
      REAL(KIND=RKIND), DIMENSION(:,:), POINTER :: highFreqThicknessCur, &
      highFreqThicknessNew
      REAL(KIND=RKIND), DIMENSION(:,:), POINTER :: lowFreqDivergenceCur, &
                                                   lowFreqDivergenceNew
      REAL(KIND=RKIND), DIMENSION(:), POINTER :: sshCur, sshNew

      REAL(KIND=RKIND), DIMENSION(:,:,:), POINTER :: tracerGroup, tracersCur, tracersNew

      ! Diagnostics Field Pointers
      TYPE(field1DReal), POINTER :: boundaryLayerDepthField, effectiveDensityField
      TYPE(field2DReal), POINTER :: normalizedRelativeVorticityEdgeField, &
                                    divergenceField, relativeVorticityField

      ! State/Tend Field Pointers
      TYPE(field2DReal), POINTER :: normalVelocityField, layerThicknessField
      TYPE(field3DReal), POINTER :: tracersGroupField

      ! Tracer Group Iteartion
      TYPE(mpas_pool_iterator_type) :: groupItr
      CHARACTER(LEN=StrKIND) :: modifiedGroupName
      CHARACTER(LEN=StrKIND) :: configName

      ! Get config options
      CALL mpas_pool_get_config(domain % configs, 'config_mom_del4', config_mom_del4)
      CALL mpas_pool_get_config(domain % configs, 'config_filter_btr_mode', &
                                config_filter_btr_mode)
      CALL mpas_pool_get_config(domain % configs, 'config_prescribe_velocity', &
                                config_prescribe_velocity)
      CALL mpas_pool_get_config(domain % configs, 'config_prescribe_thickness', &
                                config_prescribe_thickness)
      CALL mpas_pool_get_config(domain % configs, &
                                'config_use_freq_filtered_thickness', &
                                config_use_freq_filtered_thickness)
      CALL mpas_pool_get_config(domain % configs, 'config_use_standardGM', &
                                config_use_standardGM)
      CALL mpas_pool_get_config(domain % configs, 'config_use_cvmix_kpp', &
                                config_use_cvmix_kpp)
      CALL mpas_pool_get_config(domain % configs, 'config_land_ice_flux_mode', &
                                config_land_ice_flux_mode)
      !
      ! Initialize time_levs(2) with state at current time
      ! Initialize first RK state
      ! Couple tracers time_levs(2) with layerThickness in time-levels
      ! Initialize RK weights
      !
      block => domain % blocklist
      
      DO WHILE (ASSOCIATED(block))
      
         CALL mpas_pool_get_subpool(block % structs, 'state', statePool)
         CALL mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         CALL mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

         CALL mpas_pool_create_pool(provisStatePool)

         CALL mpas_pool_clone_pool(statePool, provisStatePool, 1)
         CALL mpas_pool_add_subpool(block % structs, 'provis_state', provisStatePool)
         CALL mpas_threading_barrier()

         CALL mpas_pool_get_DIMENSION(block % dimensions, 'nCells', nCells)
         CALL mpas_pool_get_DIMENSION(block % dimensions, 'nEdges', nEdges)

         CALL mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
         CALL mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
         CALL mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
         CALL mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)

         CALL mpas_pool_get_array(statePool, 'highFreqThickness', &
                                  highFreqThicknessCur, 1)
         CALL mpas_pool_get_array(statePool, 'highFreqThickness', &
                                  highFreqThicknessNew, 2)
         CALL mpas_pool_get_array(statePool, 'lowFreqDivergence', &      
                                  lowFreqDivergenceCur, 1)
         CALL mpas_pool_get_array(statePool, 'lowFreqDivergence', &
                                  lowFreqDivergenceNew, 2)

         CALL mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
         CALL mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)

         !$OMP DO SCHEDULE(runtime) PRIVATE(k)
         DO iEdge = 1, nEdges
            DO k = 1, maxLevelEdgeTop(iEdge)
               normalVelocityNew(k, iEdge) = normalVelocityCur(k, iEdge)
            END DO
         END DO
         !$OMP END DO

         !$OMP DO SCHEDULE(runtime) PRIVATE(k)
         DO iCell = 1, nCells
            DO k = 1, maxLevelCell(iCell)
               layerThicknessNew(k, iCell) = layerThicknessCur(k, iCell)
            END DO
         END DO
         !$OMP END DO

         CALL mpas_pool_begin_iteration(tracersPool)
         DO WHILE (mpas_pool_get_next_member(tracersPool, groupItr))

            IF (groupItr % memberType == MPAS_POOL_FIELD) THEN

               CALL mpas_pool_get_array(tracersPool, TRIM(groupItr % memberName), &
                                        tracersCur, 1)
               CALL mpas_pool_get_array(tracersPool, TRIM(groupItr % memberName), &
                                        tracersNew, 2)

               IF (ASSOCIATED(tracersCur) .AND. ASSOCIATED(tracersNew)) THEN
                  !$OMP DO SCHEDULE(runtime) PRIVATE(k)
                  DO iCell = 1, nCells  ! couple tracers to thickness
                     DO k = 1, maxLevelCell(iCell)
                        tracersNew(:, k, iCell) &
                        = tracersCur(:, k, iCell)*layerThicknessCur(k, iCell)
                     END DO
                  END DO
                  !$OMP END DO
               END IF
            END IF
         END DO

         IF (ASSOCIATED(highFreqThicknessCur)) THEN
              !$OMP DO SCHEDULE(runtime)
              DO iCell = 1, nCells
                 highFreqThicknessNew(:, iCell) = highFreqThicknessCur(:, iCell)
              END DO
              !$OMP END DO
         END IF

         IF (ASSOCIATED(lowFreqDivergenceCur)) THEN
              !$OMP DO SCHEDULE(runtime)
              DO iCell = 1, nCells
                 lowFreqDivergenceNew(:, iCell) = lowFreqDivergenceCur(:, iCell)
              END DO
              !$OMP END DO
         END IF

         block => block % next
      END DO

      block => domain % blocklist
      
      DO WHILE(ASSOCIATED(block))
      
         IF (ASSOCIATED(block % prev)) THEN
            CALL mpas_pool_get_subpool(block % prev % structs, 'provis_state', &
                                       prevProvisPool)
         ELSE
            nullify(prevProvisPool)
         END IF

         IF (ASSOCIATED(block % next)) THEN
            CALL mpas_pool_get_subpool(block % next % structs, 'provis_state', &
                                       nextProvisPool)
         ELSE
            nullify(nextProvisPool)
         END IF

         CALL mpas_pool_get_subpool(block % structs, 'provis_state', provisStatePool)

         IF (ASSOCIATED(prevProvisPool) .AND. ASSOCIATED(nextProvisPool)) THEN
            CALL mpas_pool_link_pools(provisStatePool, prevProvisPool, nextProvisPool)
         ELSE IF (ASSOCIATED(prevProvisPool)) THEN
            CALL mpas_pool_link_pools(provisStatePool, prevProvisPool)
         ELSE IF (ASSOCIATED(nextProvisPool)) THEN
            CALL mpas_pool_link_pools(provisStatePool, nextPool=nextProvisPool)
         ELSE
            CALL mpas_pool_link_pools(provisStatePool)
         END IF

         CALL mpas_pool_link_parinfo(block, provisStatePool)

         block => block % next
         
      END DO

      CALL mpas_threading_barrier()

      ! Fourth-order Runge-Kutta, solving dy/dt = f(t,y) is typically written as follows
      ! where h = delta t is the large time step.  Here f(t,y) is the right hand side,
      ! called the tendencies in the code below.
      ! k_1 = h f(t_n        , y_n)
      ! k_2 = h f(t_n + 1/2 h, y_n + 1/2 k_1)
      ! k_3 = h f(t_n + 1/2 h, y_n + 1/2 k_2)
      ! k_4 = h f(t_n +     h, y_n +     k_3)
      ! y_{n+1} = y_n + 1/6 k_1 + 1/3 k_2 + 1/3 k_3 + 1/6 k_4

      ! in index notation:
      ! k_{j+1} = h f(t_n + a_j h, y_n + a_j k_j)
      ! y_{n+1} = y_n + sum (b_j k_j)

      ! The coefficients of k_j are b_j = (1/6, 1/3, 1/3, 1/6) and are
      ! initialized here as delta t*b_j:

      rk_weights(1) = dt/6.
      rk_weights(2) = dt/3.
      rk_weights(3) = dt/3.
      rk_weights(4) = dt/6.

      ! The a_j coefficients of h in the computation of k_j are typically written 
      ! (0, 1/2, 1/2, 1). However, in the algorithm below we pre-compute the state for 
      ! the tendency one iteration early. That is, on j=1 (rk_step=1, below) we pre-
      ! compute y_n + 1/2 k_1 and save it in provis_state. Then we compute 1/6 k_1 and 
      ! add it to state % time_levs(2). That is why the coefficients of h are one index
      ! early in the following, i.e. a = (1/2, 1/2, 1)

      rk_substep_weights(1) = dt/2.
      rk_substep_weights(2) = dt/2.
      rk_substep_weights(3) = dt
      rk_substep_weights(4) = dt ! a_4 only used for ALE step, otherwise it is skipped.

      CALL mpas_timer_start("RK4-main loop")
      !
      ! Begin RK Loop
      !
      DO rk_step = 1, 4
      
        CALL mpas_pool_get_subpool(domain % blocklist % structs, 'diagnostics', &
                                   diagnosticsPool)

        ! Update halos for diagnostic variables.
        IF (config_use_cvmix_kpp) THEN
           CALL mpas_timer_start("RK4-boundary layer depth halo update")
           CALL mpas_dmpar_field_halo_exch(domain, 'boundaryLayerDepth')
           CALL mpas_timer_stop("RK4-boundary layer depth halo update")
        END IF


        CALL mpas_timer_start("RK4-diagnostic halo update")

        CALL mpas_dmpar_field_halo_exch(domain, 'normalizedRelativeVorticityEdge')
        IF (config_mom_del4 > 0.0_RKIND) THEN
           CALL mpas_dmpar_field_halo_exch(domain, 'divergence')
           CALL mpas_dmpar_field_halo_exch(domain, 'relativeVorticity')
        END IF
        CALL mpas_timer_stop("RK4-diagnostic halo update")
        CALL mpas_threading_barrier()

        ! Compute tendencies for high frequency thickness
        ! In RK4 notation, we are computing the right hand side f(t,y),
        ! which is the same as k_j / h.

        IF (config_use_freq_filtered_thickness) THEN
        
           CALL mpas_timer_start("RK4-tendency computations")

           block => domain % blocklist
           
           DO WHILE (ASSOCIATED(block))
              CALL mpas_pool_get_subpool(block % structs, 'tend', tendPool)
              CALL mpas_pool_get_subpool(block % structs, 'diagnostics', &
                                         diagnosticsPool)
              CALL mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
              CALL mpas_pool_get_subpool(block % structs, 'state', statePool)
              CALL mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
              CALL mpas_pool_get_subpool(block % structs, 'provis_state', &
                                         provisStatePool)
              CALL ocn_tend_freq_filtered_thickness(tendPool, provisStatePool, &
                                                    diagnosticsPool, meshPool, 1)
              CALL mpas_threading_barrier()
              block => block % next
           END DO

           CALL mpas_timer_stop("RK4-tendency computations")

           CALL mpas_timer_start("RK4-prognostic halo update")

           CALL mpas_dmpar_field_halo_exch(domain, 'tendHighFreqThickness')
           CALL mpas_dmpar_field_halo_exch(domain, 'tendLowFreqDivergence')

           CALL mpas_timer_stop("RK4-prognostic halo update")
           CALL mpas_threading_barrier()

           ! Compute next substep state for high frequency thickness.
           ! In RK4 notation, we are computing y_n + a_j k_j.
           block => domain % blocklist
           
           DO WHILE (ASSOCIATED(block))
              CALL mpas_pool_get_subpool(block % structs, 'state', statePool)
              CALL mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
              CALL mpas_pool_get_subpool(block % structs, 'tend', tendPool)
              CALL mpas_pool_get_subpool(block % structs, 'provis_state', &
                                         provisStatePool)
              CALL mpas_pool_get_array(statePool, 'highFreqThickness', &
                                       highFreqThicknessCur, 1)
              CALL mpas_pool_get_array(provisStatePool, 'highFreqThickness', &
                                       highFreqThicknessProvis, 1)
              CALL mpas_pool_get_array(tendPool, 'highFreqThickness', &
                                       highFreqThicknessTend)
              CALL mpas_pool_get_DIMENSION(block % dimensions, 'nCells', nCells)
              !$OMP DO SCHEDULE(runtime)
              DO iCell = 1, nCells
                 highFreqThicknessProvis(:, iCell) &
                 = highFreqThicknessCur(:, iCell) + rk_substep_weights(rk_step) &
                   *highFreqThicknessTend(:, iCell)
              END DO
              !$OMP END DO
              CALL mpas_threading_barrier()
              block => block % next
           END DO

        END IF 

        ! Compute tendencies for velocity, thickness, and tracers.
        ! In RK4 notation, we are computing the right hand side f(t,y),
        ! which is the same as k_j / h.
        CALL mpas_timer_start("RK4 vel/thick tendency computations")

        block => domain % blocklist
        DO WHILE (ASSOCIATED(block))
           CALL &
           ocn_time_integrator_rk4_compute_vel_tends(block, dt, &
                                                     rk_substep_weights(rk_step), err)
           CALL &
           ocn_time_integrator_rk4_compute_thick_tends(block, dt, &
                                                       rk_substep_weights(rk_step), err)
           block => block % next
        END DO

        CALL mpas_timer_stop("RK4 vel/thick tendency computations")

        ! Update halos for prognostic variables.

        CALL mpas_timer_start("RK4 vel/thick prognostic halo update")

        CALL mpas_dmpar_field_halo_exch(domain, 'tendNormalVelocity')
        CALL mpas_dmpar_field_halo_exch(domain, 'tendLayerThickness')

        CALL mpas_timer_stop("RK4 vel/thick prognostic halo update")

        CALL mpas_timer_start("RK4 tracer tendency computations")

        block => domain % blocklist
        
        DO WHILE (ASSOCIATED(block))
           CALL &
           ocn_time_integrator_rk4_compute_tracer_tends(block, dt, &
                                                        rk_substep_weights(rk_step), &
                                                        err)
           block => block % next
        END DO

        CALL mpas_timer_stop("RK4 tracer tendency computations")

        CALL mpas_timer_start("RK4 tracer prognostic halo update")

        CALL mpas_pool_get_subpool(domain % blocklist % structs, 'tend', tendPool)
        CALL mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)

        CALL mpas_pool_begin_iteration(tracersTendPool)
        DO WHILE (mpas_pool_get_next_member(tracersTendPool, groupItr))
           IF (groupItr % memberType == MPAS_POOL_FIELD) THEN
              CALL mpas_dmpar_field_halo_exch(domain, TRIM(groupItr % memberName))
           END IF
        END DO

        CALL mpas_timer_stop("RK4 tracer prognostic halo update")
        CALL mpas_threading_barrier()

        ! Compute next substep state for velocity, thickness, and tracers.
        ! In RK4 notation, we are computing y_n + a_j k_j.

        CALL mpas_timer_start("RK4-update diagnostic variables")

        IF (rk_step < 4) THEN
           block => domain % blocklist
           DO WHILE (ASSOCIATED(block))
              CALL &
              ocn_time_integrator_rk4_diagnostic_update(block, dt, &
                                                        rk_substep_weights(rk_step), &
                                                        err)
              block => block % next
           END DO
        END IF

        CALL mpas_timer_stop("RK4-update diagnostic variables")
        CALL mpas_threading_barrier()

        ! Accumulate update.
        ! In RK4 notation, we are computing b_j k_j and adding it to an accumulating sum
        ! so that we have
        !    y_{n+1} = y_n + sum (b_j k_j)
        ! after the fourth iteration.

        CALL mpas_timer_start("RK4-RK4 accumulate update")

        block => domain % blocklist
        
        DO WHILE (ASSOCIATED(block))
           CALL ocn_time_integrator_rk4_accumulate_update(block, rk_weights(rk_step), &            
                                                          err)
           block => block % next
        END DO

        CALL mpas_timer_stop("RK4-RK4 accumulate update")
        CALL mpas_threading_barrier()

      END DO
      !
      ! End RK Loop
      !
      CALL mpas_timer_stop("RK4-main loop")
      CALL mpas_threading_barrier()
      !
      ! A little clean up at the end: rescale tracer fields and compute diagnostics for
      ! new state
      !
      CALL mpas_timer_start("RK4-cleaup phase")

      ! Rescale tracers
      block => domain % blocklist
      DO WHILE(ASSOCIATED(block))
        CALL ocn_time_integrator_rk4_cleanup(block, dt, err)

        block => block % next
      END DO

      CALL mpas_timer_start("RK4-implicit vert mix")
      ! Update halo on u and tracers, which were just updated for implicit vertical 
      ! mixing. If not done, this leads to lack of volume conservation. It is required 
      ! because halo updates in RK4 are only conducted on tendencies, not on the
      ! velocity and tracer fields.  So this update is required to communicate the 
      ! change due to implicit vertical mixing across the boundary.
      CALL mpas_timer_start("RK4-implicit vert mix halos")

      CALL mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)
      CALL mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

      CALL mpas_dmpar_field_halo_exch(domain, 'normalVelocity', timeLevel=2)

      CALL mpas_pool_begin_iteration(tracersPool)
      DO WHILE (mpas_pool_get_next_member(tracersPool, groupItr))
         IF (groupItr % memberType == MPAS_POOL_FIELD) THEN
            CALL mpas_dmpar_field_halo_exch(domain, groupItr % memberName, timeLevel=2)
         END IF
      END DO

      CALL mpas_timer_stop("RK4-implicit vert mix halos")

      CALL mpas_timer_stop("RK4-implicit vert mix")
      CALL mpas_threading_barrier()

      block => domain % blocklist
      
      DO WHILE (ASSOCIATED(block))
      
         CALL mpas_pool_get_subpool(block % structs, 'state', statePool)
         CALL mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         CALL mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
         CALL mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         CALL mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         CALL mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)

         CALL mpas_pool_get_DIMENSION(meshPool, 'nCells', nCells)
         CALL mpas_pool_get_DIMENSION(meshPool, 'nEdges', nEdges)

         CALL mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
         CALL mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
         CALL mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
         CALL mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)

         CALL mpas_pool_get_DIMENSION(diagnosticsPool, 'index_surfaceVelocityZonal', &
                                      indexSurfaceVelocityZonal)
         CALL mpas_pool_get_DIMENSION(diagnosticsPool, &
                                      'index_surfaceVelocityMeridional', &
                                      indexSurfaceVelocityMeridional)
         CALL mpas_pool_get_DIMENSION(diagnosticsPool, 'index_SSHGradientZonal', &
                                      indexSSHGradientZonal)
         CALL mpas_pool_get_DIMENSION(diagnosticsPool, 'index_SSHGradientMeridional', &
                                      indexSSHGradientMeridional)
         CALL mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', &
                                  normalTransportVelocity)
         CALL mpas_pool_get_array(diagnosticsPool, 'normalGMBolusVelocity', &
                                  normalGMBolusVelocity)
         CALL mpas_pool_get_array(diagnosticsPool, 'velocityX', velocityX)
         CALL mpas_pool_get_array(diagnosticsPool, 'velocityY', velocityY)
         CALL mpas_pool_get_array(diagnosticsPool, 'velocityZ', velocityZ)
         CALL mpas_pool_get_array(diagnosticsPool, 'velocityZonal', velocityZonal)
         CALL mpas_pool_get_array(diagnosticsPool, 'velocityMeridional', &
                                  velocityMeridional)
         CALL mpas_pool_get_array(diagnosticsPool, 'gradSSH', gradSSH)
         CALL mpas_pool_get_array(diagnosticsPool, 'gradSSHX', gradSSHX)
         CALL mpas_pool_get_array(diagnosticsPool, 'gradSSHY', gradSSHY)
         CALL mpas_pool_get_array(diagnosticsPool, 'gradSSHZ', gradSSHZ)
         CALL mpas_pool_get_array(diagnosticsPool, 'gradSSHZonal', gradSSHZonal)
         CALL mpas_pool_get_array(diagnosticsPool, 'gradSSHMeridional', &
                                  gradSSHMeridional)
         CALL mpas_pool_get_array(diagnosticsPool, 'surfaceVelocity', surfaceVelocity)
         CALL mpas_pool_get_array(diagnosticsPool, 'SSHGradient', SSHGradient)

         IF (config_prescribe_velocity) THEN
            !$OMP DO SCHEDULE(runtime)
            DO iEdge = 1, nEdges
               normalVelocityNew(:, iEdge) = normalVelocityCur(:, iEdge)
            END DO
            !$OMP END DO
         END IF

         IF (config_prescribe_thickness) THEN
            !$OMP DO SCHEDULE(runtime)
            DO iCell = 1, nCells
               layerThicknessNew(:, iCell) = layerThicknessCur(:, iCell)
            END DO
            !$OMP END DO
         END IF

         CALL ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool, &
                                   diagnosticsPool, scratchPool, tracersPool, 2)
         CALL mpas_threading_barrier()

         ! Update the effective desnity in land ice if we're coupling to land ice
         CALL ocn_effective_density_in_land_ice_update(meshPool, forcingPool, &
                                                       statePool, scratchPool, err)

         !-----------------------------------------------------------------------------!
         ! Accumulating various parameterizations of the transport velocity
         !-----------------------------------------------------------------------------!
         
         !$OMP DO SCHEDULE(runtime)
         DO iEdge = 1, nEdges
            normalTransportVelocity(:, iEdge) = normalVelocityNew(:, iEdge)
         END DO
         !$OMP END DO
         
         CALL mpas_threading_barrier()

         ! Compute normalGMBolusVelocity and the tracer transport velocity
         IF (config_use_standardGM) THEN
             CALL ocn_gm_compute_Bolus_velocity(diagnosticsPool, meshPool, scratchPool)
         END IF
         CALL mpas_threading_barrier()

         IF (config_use_standardGM) THEN
            !$OMP DO SCHEDULE(runtime)
            DO iEdge = 1, nEdges
               normalTransportVelocity(:, iEdge) &
               = normalTransportVelocity(:, iEdge) + normalGMBolusVelocity(:, iEdge)
            END DO
            !$OMP END DO
         END IF
         
         !-----------------------------------------------------------------------------!
         ! End: Accumulating various parameterizations of the transport velocity
         !-----------------------------------------------------------------------------!

         CALL mpas_reconstruct(meshPool, normalVelocityNew, &
                          velocityX, velocityY, velocityZ,   &
                          velocityZonal, velocityMeridional, &
                          includeHalos = .TRUE.)

         CALL mpas_reconstruct(meshPool, gradSSH,          &
                          gradSSHX, gradSSHY, gradSSHZ,    &
                          gradSSHZonal, gradSSHMeridional, &
                          includeHalos = .TRUE.)
                          
         CALL mpas_threading_barrier()

         !$OMP DO SCHEDULE(runtime)
         DO iCell = 1, nCells
            surfaceVelocity(indexSurfaceVelocityZonal, iCell) = velocityZonal(1, iCell)
            surfaceVelocity(indexSurfaceVelocityMeridional, iCell) &
            = velocityMeridional(1, iCell)
            SSHGradient(indexSSHGradientZonal, iCell) = gradSSHZonal(iCell)
            SSHGradient(indexSSHGradientMeridional, iCell) = gradSSHMeridional(iCell)
         END DO
         !$OMP END DO

         CALL ocn_time_average_coupled_accumulate(diagnosticsPool, statePool, &
                                                  forcingPool, 2)

         IF (config_use_standardGM) THEN
            CALL ocn_reconstruct_gm_vectors(diagnosticsPool, meshPool)
         END IF
         
         CALL mpas_threading_barrier()

         block => block % next
         
      END DO

      IF (TRIM(config_land_ice_flux_mode) == 'coupled') THEN
         CALL mpas_timer_start("RK4-effective density halo")
         CALL mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)
         CALL mpas_pool_get_field(statePool, 'effectiveDensityInLandIce', &
                                  effectiveDensityField, 2)
         CALL mpas_dmpar_exch_halo_field(effectiveDensityField)
         CALL mpas_timer_stop("RK4-effective density halo")
      END IF

      CALL mpas_timer_stop("RK4-cleaup phase")

      CALL mpas_threading_barrier()

      block => domain % blocklist
      
      DO WHILE(ASSOCIATED(block))
         CALL mpas_pool_get_subpool(block % structs, 'provis_state', provisStatePool)

         CALL mpas_pool_destroy_pool(provisStatePool)

         CALL mpas_pool_remove_subpool(block % structs, 'provis_state')
         block => block % next
      END DO
      
      CALL mpas_threading_barrier()

   END SUBROUTINE ocn_time_integrator_rk4
   !-----------------------------------------------------------------------------------!
   !                                                                                   !
   !-----------------------------------------------------------------------------------!
   SUBROUTINE ocn_time_integrator_rk4_compute_vel_tends(block, dt, rkWeight, err)
   
      IMPLICIT NONE
   
      TYPE(block_type), INTENT(IN) :: block
      REAL(KIND=RKIND), INTENT(IN) :: dt
      REAL(KIND=RKIND), INTENT(IN) :: rkWeight
      INTEGER, INTENT(OUT) :: err
      
      ! Local variables
      LOGICAL, POINTER :: config_filter_btr_mode
      REAL(KIND=RKIND), DIMENSION(:), POINTER :: sshCur
      REAL(KIND=RKIND), DIMENSION(:,:), POINTER &
      :: layerThicknessCur, normalVelocityCur, layerThicknessEdge, &
         vertAleTransportTop, normalTransportVelocity, normalVelocityProvis, &
         highFreqThicknessProvis
      TYPE(mpas_pool_type), POINTER :: meshPool, verticalMeshPool, statePool, &
                                       diagnosticsPool, forcingPool, scratchPool, &
                                       tendPool, provisStatePool, tracersPool

      err = 0

      CALL mpas_pool_get_config(block % configs, 'config_filter_btr_mode', &
                                config_filter_btr_mode)

      CALL mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      CALL mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)
      CALL mpas_pool_get_subpool(block % structs, 'state', statePool)
      CALL mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
      CALL mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
      CALL mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
      CALL mpas_pool_get_subpool(block % structs, 'tend', tendPool)
      CALL mpas_pool_get_subpool(block % structs, 'provis_state', provisStatePool)

      CALL mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

      CALL mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
      CALL mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
      CALL mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)

      CALL mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', &
                               layerThicknessEdge)
      CALL mpas_pool_get_array(diagnosticsPool, 'vertAleTransportTop', &
                               vertAleTransportTop)
      CALL mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', &
                               normalTransportVelocity)

      CALL mpas_pool_get_array(provisStatePool, 'normalVelocity', &
                               normalVelocityProvis, 1)
      CALL mpas_pool_get_array(provisStatePool, 'highFreqThickness', &
                               highFreqThicknessProvis, 1)

      ! advection of u uses u, while advection of layerThickness and tracers use 
      ! normalTransportVelocity.
      IF (ASSOCIATED(highFreqThicknessProvis)) THEN
         CALL ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, scratchPool, &
            layerThicknessCur,layerThicknessEdge, normalVelocityProvis, &
            sshCur, rkWeight, &
            vertAleTransportTop, err, highFreqThicknessProvis)
      ELSE
         CALL ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, scratchPool, &
            layerThicknessCur,layerThicknessEdge, normalVelocityProvis, &
            sshCur, rkWeight, &
            vertAleTransportTop, err)
      END IF 
      CALL mpas_threading_barrier()

      CALL ocn_tend_vel(tendPool, provisStatePool, forcingPool, diagnosticsPool, &
                        meshPool, scratchPool, 1, dt)
      CALL mpas_threading_barrier()

   END SUBROUTINE ocn_time_integrator_rk4_compute_vel_tends
   !-----------------------------------------------------------------------------------!
   !                                                                                   !
   !-----------------------------------------------------------------------------------!
   SUBROUTINE ocn_time_integrator_rk4_compute_thick_tends(block, dt, rkWeight, err)
   
      IMPLICIT NONE 
   
      TYPE(block_type), INTENT(IN) :: block
      REAL(KIND=RKIND), INTENT(IN) :: dt
      REAL(KIND=RKIND), INTENT(IN) :: rkWeight
      INTEGER, INTENT(OUT) :: err
   
      ! Local variables
      LOGICAL, POINTER :: config_filter_btr_mode
      REAL(KIND=RKIND), DIMENSION(:), POINTER :: sshCur
      REAL(KIND=RKIND), DIMENSION(:,:), POINTER &
      :: layerThicknessCur, normalVelocityCur, layerThicknessEdge, &
         vertAleTransportTop, normalTransportVelocity, normalVelocityProvis, &
         highFreqThicknessProvis
      TYPE(mpas_pool_type), POINTER :: meshPool, verticalMeshPool, statePool, &
                                       diagnosticsPool, forcingPool, scratchPool, &
                                       tendPool, provisStatePool, tracersPool

      err = 0

      CALL mpas_pool_get_config(block % configs, 'config_filter_btr_mode', &
                              config_filter_btr_mode)

      CALL mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      CALL mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)
      CALL mpas_pool_get_subpool(block % structs, 'state', statePool)
      CALL mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
      CALL mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
      CALL mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
      CALL mpas_pool_get_subpool(block % structs, 'tend', tendPool)
      CALL mpas_pool_get_subpool(block % structs, 'provis_state', provisStatePool)

      CALL mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

      CALL mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
      CALL mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
      CALL mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)

      CALL mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', &
                               layerThicknessEdge)
      CALL mpas_pool_get_array(diagnosticsPool, 'vertAleTransportTop', &
                               vertAleTransportTop)
      CALL mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', &
                               normalTransportVelocity)

      CALL mpas_pool_get_array(provisStatePool, 'normalVelocity', &
                               normalVelocityProvis, 1)
      CALL mpas_pool_get_array(provisStatePool, 'highFreqThickness', &
                               highFreqThicknessProvis, 1)

      ! Advection of u uses u, while advection of layerThickness and tracers use 
      ! normalTransportVelocity.
      IF (ASSOCIATED(highFreqThicknessProvis)) THEN
         CALL ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, scratchPool, &
            layerThicknessCur, layerThicknessEdge, normalTransportVelocity, &
            sshCur, rkWeight, &
            vertAleTransportTop, err, highFreqThicknessProvis)
      ELSE
         CALL ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, scratchPool, &
            layerThicknessCur, layerThicknessEdge, normalTransportVelocity, &
            sshCur, rkWeight, &
            vertAleTransportTop, err)
      END IF 
      
      CALL mpas_threading_barrier()

      CALL ocn_tend_thick(tendPool, forcingPool, diagnosticsPool, meshPool)

      CALL mpas_threading_barrier()

   END SUBROUTINE ocn_time_integrator_rk4_compute_thick_tends
   !-----------------------------------------------------------------------------------!
   !                                                                                   !
   !-----------------------------------------------------------------------------------!
   SUBROUTINE ocn_time_integrator_rk4_compute_tracer_tends(block, dt, rkWeight, err)
   
      IMPLICIT NONE 
   
      TYPE(block_type), INTENT(IN) :: block
      REAL(KIND=RKIND), INTENT(IN) :: dt
      REAL(KIND=RKIND), INTENT(IN) :: rkWeight
      INTEGER, INTENT(OUT) :: err
      
      ! Local variables
      LOGICAL, POINTER :: config_filter_btr_mode
      REAL(KIND=RKIND), DIMENSION(:), POINTER :: sshCur
      REAL(KIND=RKIND), DIMENSION(:,:), POINTER &
      :: layerThicknessCur, normalVelocityCur, layerThicknessEdge, &
         vertAleTransportTop, normalTransportVelocity, normalVelocityProvis, &
         highFreqThicknessProvis
      TYPE(mpas_pool_type), POINTER :: meshPool, verticalMeshPool, statePool, &
                                       diagnosticsPool, forcingPool, scratchPool, &
                                       tendPool, provisStatePool, swForcingPool, &
                                       tracersPool
                                       
      err = 0

      CALL mpas_pool_get_config(block % configs, 'config_filter_btr_mode', &
                                config_filter_btr_mode)

      CALL mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      CALL mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)
      CALL mpas_pool_get_subpool(block % structs, 'state', statePool)
      CALL mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
      CALL mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
      CALL mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
      CALL mpas_pool_get_subpool(block % structs, 'tend', tendPool)
      CALL mpas_pool_get_subpool(block % structs, 'provis_state', provisStatePool)
      CALL mpas_pool_get_subpool(block % structs, 'shortwave', swForcingPool)

      CALL mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

      CALL mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
      CALL mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
      CALL mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)

      CALL mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', &
                               layerThicknessEdge)
      CALL mpas_pool_get_array(diagnosticsPool, 'vertAleTransportTop', &
                               vertAleTransportTop)
      CALL mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', &
                               normalTransportVelocity)

      CALL mpas_pool_get_array(provisStatePool, 'normalVelocity', &
                               normalVelocityProvis, 1)
      CALL mpas_pool_get_array(provisStatePool, 'highFreqThickness', &
                               highFreqThicknessProvis, 1)

      ! Advection of u uses u, while advection of layerThickness and tracers use 
      ! normalTransportVelocity.
      IF (ASSOCIATED(highFreqThicknessProvis)) THEN
         CALL ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, scratchPool, &
            layerThicknessCur, layerThicknessEdge, normalTransportVelocity, &
            sshCur, rkWeight, &
            vertAleTransportTop, err, highFreqThicknessProvis)
      ELSE
         CALL ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, scratchPool, &
            layerThicknessCur, layerThicknessEdge, normalTransportVelocity, &
            sshCur, rkWeight, &
            vertAleTransportTop, err)
      END IF 
      
      CALL mpas_threading_barrier()

      IF (config_filter_btr_mode) THEN
          CALL ocn_filter_btr_mode_tend_vel(tendPool, provisStatePool, &
                                            diagnosticsPool, meshPool, 1)
      END IF 
      CALL mpas_threading_barrier()

      CALL ocn_tend_tracer(tendPool, provisStatePool, forcingPool, diagnosticsPool, &
                           meshPool, swForcingPool, scratchPool, dt, &
                           activeTracersOnlyIn=.FALSE., timeLevelIn=1)
                           
      CALL mpas_threading_barrier()

   END SUBROUTINE ocn_time_integrator_rk4_compute_tracer_tends
   !-----------------------------------------------------------------------------------!
   !                                                                                   !
   !-----------------------------------------------------------------------------------!
   SUBROUTINE ocn_time_integrator_rk4_compute_tends(block, dt, rkWeight, err)
   
      IMPLICIT NONE 
   
      TYPE(block_type), INTENT(IN) :: block
      REAL(KIND=RKIND), INTENT(IN) :: dt
      REAL(KIND=RKIND), INTENT(IN) :: rkWeight
      INTEGER, INTENT(OUT) :: err
      
      ! Local variables
      LOGICAL, POINTER :: config_filter_btr_mode
      REAL(KIND=RKIND), DIMENSION(:), POINTER :: sshCur
      REAL(KIND=RKIND), DIMENSION(:,:), POINTER :: layerThicknessCur, normalVelocityCur
      REAL(KIND=RKIND), DIMENSION(:,:), POINTER :: layerThicknessEdge, &
                                                   vertAleTransportTop
      REAL(KIND=RKIND), DIMENSION(:,:), POINTER :: normalTransportVelocity
      REAL(KIND=RKIND), DIMENSION(:,:), POINTER ::  normalVelocityProvis, &
                                                    highFreqThicknessProvis                                                    
      TYPE(mpas_pool_type), POINTER :: meshPool, verticalMeshPool, statePool, &
                                       diagnosticsPool, forcingPool, scratchPool, &
                                       tendPool, provisStatePool, swForcingPool, &
                                       tracersPool
                                       
      err = 0

      CALL mpas_pool_get_config(block % configs, 'config_filter_btr_mode', &
                                config_filter_btr_mode)

      CALL mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      CALL mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)
      CALL mpas_pool_get_subpool(block % structs, 'state', statePool)
      CALL mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
      CALL mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
      CALL mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
      CALL mpas_pool_get_subpool(block % structs, 'tend', tendPool)
      CALL mpas_pool_get_subpool(block % structs, 'provis_state', provisStatePool)
      CALL mpas_pool_get_subpool(block % structs, 'shortwave', swForcingPool)

      CALL mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

      CALL mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
      CALL mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
      CALL mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)

      CALL mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', &
                               layerThicknessEdge)
      CALL mpas_pool_get_array(diagnosticsPool, 'vertAleTransportTop', &
                               vertAleTransportTop)
      CALL mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', &
                               normalTransportVelocity)

      CALL mpas_pool_get_array(provisStatePool, 'normalVelocity', &
                               normalVelocityProvis, 1)
      CALL mpas_pool_get_array(provisStatePool, 'highFreqThickness', &
                               highFreqThicknessProvis, 1)

      ! Advection of u uses u, while advection of layerThickness and tracers use 
      ! normalTransportVelocity.
      IF (ASSOCIATED(highFreqThicknessProvis)) THEN
         CALL ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, scratchPool, &
            layerThicknessCur,layerThicknessEdge, normalVelocityProvis, &
            sshCur, rkWeight, &
            vertAleTransportTop, err, highFreqThicknessProvis)
      ELSE
         CALL ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, scratchPool, &
            layerThicknessCur,layerThicknessEdge, normalVelocityProvis, &
            sshCur, rkWeight, &
            vertAleTransportTop, err)
      END IF 
      
      CALL mpas_threading_barrier()

      CALL ocn_tend_vel(tendPool, provisStatePool, forcingPool, diagnosticsPool, &
                        meshPool, scratchPool, 1, dt)
      CALL mpas_threading_barrier()

      IF (ASSOCIATED(highFreqThicknessProvis)) THEN
         CALL ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, scratchPool, &
            layerThicknessCur, layerThicknessEdge, normalTransportVelocity, &
            sshCur, rkWeight, &
            vertAleTransportTop, err, highFreqThicknessProvis)
      ELSE
         CALL ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, scratchPool, &
            layerThicknessCur, layerThicknessEdge, normalTransportVelocity, &
            sshCur, rkWeight, &
            vertAleTransportTop, err)
      END IF 
      
      CALL mpas_threading_barrier()

      CALL ocn_tend_thick(tendPool, forcingPool, diagnosticsPool, meshPool)

      IF (config_filter_btr_mode) THEN
          CALL ocn_filter_btr_mode_tend_vel(tendPool, provisStatePool, &
                                            diagnosticsPool, meshPool, 1)
      END IF 
      
      CALL mpas_threading_barrier()

      CALL ocn_tend_tracer(tendPool, provisStatePool, forcingPool, diagnosticsPool, &
                           meshPool, swForcingPool, scratchPool, dt, &
                           activeTracersOnlyIn=.FALSE., timeLevelIn=1)
                           
      CALL mpas_threading_barrier()

   END SUBROUTINE ocn_time_integrator_rk4_compute_tends
   !-----------------------------------------------------------------------------------!
   !                                                                                   !
   !-----------------------------------------------------------------------------------!
   SUBROUTINE ocn_time_integrator_rk4_diagnostic_update(block, dt, rkWeight, err)
   
      IMPLICIT NONE 
   
      TYPE(block_type), INTENT(IN) :: block
      REAL(KIND=RKIND), INTENT(IN) :: dt
      REAL(KIND=RKIND), INTENT(IN) :: rkWeight
      INTEGER, INTENT(OUT) :: err
      
      ! Local variables
      CHARACTER(LEN=StrKIND) :: modifiedGroupName, configName
      INTEGER :: iCell, iEdge, k
      INTEGER, POINTER :: nCells, nEdges
      INTEGER, DIMENSION(:), POINTER :: maxLevelCell, maxLevelEdgeTop
      LOGICAL, POINTER :: config_prescribe_velocity, config_prescribe_thickness, &
                          config_use_standardGM, config_use_tracerGroup
      REAL(KIND=RKIND), DIMENSION(:,:), POINTER &
      :: normalVelocityCur, normalVelocityProvis, normalVelocityTend, &
         layerThicknessCur, layerThicknessProvis, layerThicknessTend, &
         lowFreqDivergenceCur, lowFreqDivergenceProvis, lowFreqDivergenceTend, &
         normalTransportVelocity, normalGMBolusVelocity
      REAL(KIND=RKIND), DIMENSION(:,:,:), POINTER &
      :: tracersGroupCur, tracersGroupProvis, tracersGroupTend
      TYPE(mpas_pool_iterator_type) :: groupItr
      TYPE(mpas_pool_type), POINTER :: statePool, tendPool, meshPool, scratchPool, &
                                       diagnosticsPool, provisStatePool, forcingPool, &
                                       tracersPool, tracersTendPool, provisTracersPool      

      err = 0

      CALL mpas_pool_get_config(block % configs, 'config_prescribe_velocity', &
                                config_prescribe_velocity)
      CALL mpas_pool_get_config(block % configs, 'config_prescribe_thickness', &
                                config_prescribe_thickness)
      CALL mpas_pool_get_config(block % configs, 'config_use_standardGM', &
                                config_use_standardGM)

      CALL mpas_pool_get_DIMENSION(block % dimensions, 'nCells', nCells)
      CALL mpas_pool_get_DIMENSION(block % dimensions, 'nEdges', nEdges)

      CALL mpas_pool_get_subpool(block % structs, 'state', statePool)
      CALL mpas_pool_get_subpool(block % structs, 'tend', tendPool)
      CALL mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      CALL mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
      CALL mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
      CALL mpas_pool_get_subpool(block % structs, 'provis_state', provisStatePool)
      CALL mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)

      CALL mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
      CALL mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)

      CALL mpas_pool_get_subpool(provisStatePool, 'tracers', provisTracersPool)

      CALL mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
      CALL mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
      CALL mpas_pool_get_array(statePool, 'lowFreqDivergence', lowFreqDivergenceCur, 1)

      CALL mpas_pool_get_array(provisStatePool, 'normalVelocity', &
                               normalVelocityProvis, 1)
      CALL mpas_pool_get_array(provisStatePool, 'layerThickness', &
                               layerThicknessProvis, 1)
      CALL mpas_pool_get_array(provisStatePool, 'lowFreqDivergence', &
                               lowFreqDivergenceProvis, 1)

      CALL mpas_pool_get_array(tendPool, 'normalVelocity', normalVelocityTend)
      CALL mpas_pool_get_array(tendPool, 'layerThickness', layerThicknessTend)
      CALL mpas_pool_get_array(tendPool, 'lowFreqDivergence', lowFreqDivergenceTend)

      CALL mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      CALL mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)

      CALL mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', &
                               normalTransportVelocity)
      CALL mpas_pool_get_array(diagnosticsPool, 'normalGMBolusVelocity', &
                               normalGMBolusVelocity)

      CALL mpas_threading_barrier()

      !$OMP DO SCHEDULE(runtime) PRIVATE(k)
      DO iEdge = 1, nEdges
         DO k = 1, maxLevelEdgeTop(iEdge)
            normalVelocityProvis(k, iEdge) = normalVelocityCur(k, iEdge) + rkWeight &
                                             *normalVelocityTend(k, iEdge)
         END DO
      END DO
      !$OMP END DO

      !$OMP DO SCHEDULE(runtime) PRIVATE(k)
      DO iCell = 1, nCells
         DO k = 1, maxLevelCell(iCell)
            layerThicknessProvis(k, iCell) = layerThicknessCur(k, iCell) + rkWeight &
                                             *layerThicknessTend(k, iCell)
         END DO
      END DO
      !$OMP END DO

      CALL mpas_pool_begin_iteration(tracersPool)
      DO WHILE (mpas_pool_get_next_member(tracersPool, groupItr))
         IF (groupItr % memberType == MPAS_POOL_FIELD) THEN
            configName = 'config_use_' // TRIM(groupItr % memberName)
            CALL mpas_pool_get_config(block % configs, configName, &
                                      config_use_tracerGroup)

            IF (config_use_tracerGroup) THEN
               CALL mpas_pool_get_array(tracersPool, groupItr % memberName, &
                                        tracersGroupCur, 1)
               CALL mpas_pool_get_array(provisTracersPool, groupItr % memberName, &
                                        tracersGroupProvis, 1)

               modifiedGroupName = TRIM(groupItr % memberName) // 'Tend'
               CALL mpas_pool_get_array(tracersTendPool, modifiedGroupName, &
                                        tracersGroupTend)
               IF (ASSOCIATED(tracersGroupProvis) .AND. &
                   ASSOCIATED(tracersGroupCur) .AND. ASSOCIATED(tracersGroupTend)) THEN
                  !$OMP DO SCHEDULE(runtime) PRIVATE(k)
                  DO iCell = 1, nCells
                     DO k = 1, maxLevelCell(iCell)
                        tracersGroupProvis(:, k, iCell) &
                        = (layerThicknessCur(k, iCell)*tracersGroupCur(:, k, iCell)  &
                           + rkWeight*tracersGroupTend(:, k, iCell)) &
                          /layerThicknessProvis(k, iCell)
                     END DO

                  END DO
                  !$OMP END DO
               END IF
            END IF
         END IF
      END DO

      IF (ASSOCIATED(lowFreqDivergenceCur)) THEN
         !$OMP DO SCHEDULE(runtime)
         DO iCell = 1, nCells
            lowFreqDivergenceProvis(:, iCell) &
            = lowFreqDivergenceCur(:, iCell) + rkWeight*lowFreqDivergenceTend(:, iCell)
         END DO
         !$OMP END DO
      END IF

      IF (config_prescribe_velocity) THEN
         !$OMP DO SCHEDULE(runtime)
         DO iEdge = 1, nEdges
            normalVelocityProvis(:, iEdge) = normalVelocityCur(:, iEdge)
         END DO
         !$OMP END DO
      END IF

      IF (config_prescribe_thickness) THEN
         !$OMP DO SCHEDULE(runtime)
         DO iCell = 1, nCells
            layerThicknessProvis(:, iCell) = layerThicknessCur(:, iCell)
         END DO
         !$OMP END DO
      END IF
      CALL mpas_threading_barrier()

      CALL ocn_diagnostic_solve(dt, provisStatePool, forcingPool, meshPool, &
                                diagnosticsPool, scratchPool, tracersPool, 1)
      CALL mpas_threading_barrier()

      ! -------------------------------------------------------------------------------!
      ! Accumulating various parametrizations of the transport velocity
      ! -------------------------------------------------------------------------------!
      
      !$OMP DO SCHEDULE(runtime)
      DO iEdge = 1, nEdges
         normalTransportVelocity(:, iEdge) = normalVelocityProvis(:, iEdge)
      END DO
      !$OMP END DO
      
      CALL mpas_threading_barrier()

      ! Compute normalGMBolusVelocity, relativeSlope and RediDiffVertCoef if respective 
      ! flags are turned on
      IF (config_use_standardGM) THEN
         CALL ocn_gm_compute_Bolus_velocity(diagnosticsPool, meshPool, scratchPool)
      END IF
      
      CALL mpas_threading_barrier()

      IF (config_use_standardGM) THEN
         !$OMP DO SCHEDULE(runtime)
         DO iEdge = 1, nEdges
            normalTransportVelocity(:, iEdge) &
            = normalTransportVelocity(:, iEdge) + normalGMBolusVelocity(:,iEdge)
         END DO
         !$OMP END DO
      END IF
      
      CALL mpas_threading_barrier()
      
      ! -------------------------------------------------------------------------------!
      ! End: Accumulating various parametrizations of the transport velocity
      ! -------------------------------------------------------------------------------!

   END SUBROUTINE ocn_time_integrator_rk4_diagnostic_update
   !-----------------------------------------------------------------------------------!
   !                                                                                   !
   !-----------------------------------------------------------------------------------!
   SUBROUTINE ocn_time_integrator_rk4_accumulate_update(block, rkWeight, err)
   
      IMPLICIT NONE
   
      TYPE(block_type), INTENT(IN) :: block
      REAL(KIND=RKIND), INTENT(IN) :: rkWeight
      INTEGER, INTENT(OUT) :: err
      
      ! Local variables
      CHARACTER(LEN=StrKIND) :: modifiedGroupName, configName
      INTEGER, POINTER :: nCells, nEdges
      INTEGER :: iCell, iEdge, k
      INTEGER, DIMENSION(:), POINTER :: maxLevelCell
      LOGICAL, POINTER :: config_use_tracerGroup
      REAL(KIND=RKIND), DIMENSION(:,:), POINTER &
      :: normalVelocityNew, normalVelocityTend, layerThicknessNew, layerThicknessTend, &
         highFreqThicknessNew, highFreqThicknessTend, lowFreqDivergenceNew, &
         lowFreqDivergenceTend
      REAL(KIND=RKIND), DIMENSION(:,:,:), POINTER :: tracersGroupNew, tracersGroupTend
      TYPE(mpas_pool_iterator_type) :: groupItr
      TYPE(mpas_pool_type), POINTER :: statePool, tendPool, meshPool, tracersPool, &
                                       tracersTendPool

      err = 0

      CALL mpas_pool_get_DIMENSION(block % dimensions, 'nCells', nCells)
      CALL mpas_pool_get_DIMENSION(block % dimensions, 'nEdges', nEdges)

      CALL mpas_pool_get_subpool(block % structs, 'state', statePool)
      CALL mpas_pool_get_subpool(block % structs, 'tend', tendPool)
      CALL mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

      CALL mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
      CALL mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)

!     CALL mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
!     CALL mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
!     CALL mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessCur, 1)
!     CALL mpas_pool_get_array(statePool, 'lowFreqDivergence', lowFreqDivergenceCur, 1)

      CALL mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
      CALL mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)
      CALL mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessNew, 2)
      CALL mpas_pool_get_array(statePool, 'lowFreqDivergence', lowFreqDivergenceNew, 2)

      CALL mpas_pool_get_array(tendPool, 'normalVelocity', normalVelocityTend)
      CALL mpas_pool_get_array(tendPool, 'layerThickness', layerThicknessTend)

      CALL mpas_pool_get_array(tendPool, 'highFreqThickness', highFreqThicknessTend)
      CALL mpas_pool_get_array(tendPool, 'lowFreqDivergence', lowFreqDivergenceTend)

      CALL mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

      !$OMP DO SCHEDULE(runtime) PRIVATE(k)
      DO iCell = 1, nCells
         DO k = 1, maxLevelCell(iCell)
            layerThicknessNew(k, iCell) = layerThicknessNew(k, iCell) &
                                          + rkWeight*layerThicknessTend(k, iCell)
         END DO
      END DO
      !$OMP END DO

      !$OMP DO SCHEDULE(runtime)
      DO iEdge = 1, nEdges
         normalVelocityNew(:, iEdge) &
         = normalVelocityNew(:, iEdge) + rkWeight*normalVelocityTend(:, iEdge)
      END DO
      !$OMP END DO

      CALL mpas_pool_begin_iteration(tracersPool)
      DO WHILE (mpas_pool_get_next_member(tracersPool, groupItr))
         IF (groupItr % memberType == MPAS_POOL_FIELD) THEN
            configName = 'config_use_' // TRIM(groupItr % memberName)
            CALL mpas_pool_get_config(block % configs, configName, &
                                      config_use_tracerGroup)
            IF (config_use_tracerGroup) THEN
               CALL mpas_pool_get_array(tracersPool, groupItr % memberName, &
                                        tracersGroupNew, 2)
               modifiedGroupName = TRIM(groupItr % memberName) // 'Tend'
               CALL mpas_pool_get_array(tracersTendPool, modifiedGroupName, &
                                        tracersGroupTend)
               IF (ASSOCIATED(tracersGroupNew) .AND. ASSOCIATED(tracersGroupTend)) THEN
                  !$OMP DO SCHEDULE(runtime) PRIVATE(k)
                  DO iCell = 1, nCells
                     DO k = 1, maxLevelCell(iCell)
                        tracersGroupNew(:, k, iCell) &
                        = tracersGroupNew(:, k, iCell) &
                          + rkWeight*tracersGroupTend(:, k, iCell)
                     END DO
                  END DO
                  !$OMP END DO
               END IF
            END IF
         END IF
      END DO

      IF (ASSOCIATED(highFreqThicknessNew)) THEN
         !$OMP DO SCHEDULE(runtime)
         DO iCell = 1, nCells
            highFreqThicknessNew(:, iCell) = highFreqThicknessNew(:, iCell) &
                                             + rkWeight*highFreqThicknessTend(:, iCell)
         END DO
         !$OMP END DO
      END IF

      IF (ASSOCIATED(lowFreqDivergenceNew)) THEN
         !$OMP DO SCHEDULE(runtime)
         DO iCell = 1, nCells
            lowFreqDivergenceNew(:, iCell) = lowFreqDivergenceNew(:, iCell) &
                                             + rkWeight*lowFreqDivergenceTend(:, iCell)
         END DO
         !$OMP END DO
      END IF

   END SUBROUTINE ocn_time_integrator_rk4_accumulate_update
   !-----------------------------------------------------------------------------------!
   !                                                                                   !
   !-----------------------------------------------------------------------------------!
   SUBROUTINE ocn_time_integrator_rk4_cleanup(block, dt, err)
   
      IMPLICIT NONE 
   
      TYPE(block_type), INTENT(IN) :: block
      REAL(KIND=RKIND), INTENT(IN) :: dt
      INTEGER, INTENT(OUT) :: err
      
      ! Local variables
      CHARACTER(LEN=StrKIND) :: modifiedGroupName, configName
      INTEGER :: iCell, iEdge, k
      INTEGER, POINTER :: nCells, nEdges, indexTemperature, indexSalinity
      INTEGER, DIMENSION(:), POINTER :: maxLevelCell     
      LOGICAL, POINTER :: config_use_tracerGroup, config_use_standardGM
      REAL(KIND=RKIND), DIMENSION(:,:), POINTER &
      :: layerThicknessNew, normalVelocityNew,  normalTransportVelocity, &
         normalGMBolusVelocity
      REAL(KIND=RKIND), DIMENSION(:,:,:), POINTER :: tracersGroupNew
      TYPE(mpas_pool_type), POINTER :: statePool, meshPool, forcingPool, &
                                       diagnosticsPool, scratchPool, tracersPool
      TYPE(mpas_pool_iterator_type) :: groupItr

      err = 0

      CALL mpas_pool_get_config(block % configs, 'config_use_standardGM', &
                                config_use_standardGM)

      CALL mpas_pool_get_DIMENSION(block % dimensions, 'nCells', nCells)
      CALL mpas_pool_get_DIMENSION(block % dimensions, 'nEdges', nEdges)

      CALL mpas_pool_get_subpool(block % structs, 'state', statePool)
      CALL mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      CALL mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
      CALL mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
      CALL mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)

      CALL mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

      CALL mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)
      CALL mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)

      CALL mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', &
                               normalTransportVelocity)
      CALL mpas_pool_get_array(diagnosticsPool, 'normalGMBolusVelocity', &
                               normalGMBolusVelocity)

      CALL mpas_pool_get_DIMENSION(tracersPool, 'index_temperature', indexTemperature)
      CALL mpas_pool_get_DIMENSION(tracersPool, 'index_salinity', indexSalinity)

      CALL mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

      CALL mpas_pool_begin_iteration(tracersPool)
      DO WHILE (mpas_pool_get_next_member(tracersPool, groupItr))
         IF (groupItr % memberType == MPAS_POOL_FIELD) THEN
            CALL mpas_pool_get_array(tracersPool, groupItr % memberName, &
                                     tracersGroupNew, 2)
            IF (ASSOCIATED(tracersGroupNew)) THEN
               !$OMP DO SCHEDULE(runtime) PRIVATE(k)
               DO iCell = 1, nCells
                 DO k = 1, maxLevelCell(iCell)
                   tracersGroupNew(:, k, iCell) = tracersGroupNew(:, k, iCell) &
                                                  /layerThicknessNew(k, iCell)
                 END DO
               END DO
               !$OMP END DO
            END IF
         END IF
      END DO

      CALL ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool, diagnosticsPool, &
                                scratchPool, tracersPool, 2)
      CALL mpas_threading_barrier()

      CALL ocn_vmix_implicit(dt, meshPool, diagnosticsPool, statePool, forcingPool, &
                             scratchPool, err, 2)
      CALL mpas_threading_barrier()

   END SUBROUTINE ocn_time_integrator_rk4_cleanup

END MODULE ocn_time_integration_rk4