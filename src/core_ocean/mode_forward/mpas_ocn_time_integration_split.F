! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!  ocn_time_integration_split
!
!> \brief MPAS ocean split explicit time integration scheme
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   September 2011
!> \details
!> This module contains the routine for the split explicit time integration scheme.
!
MODULE ocn_time_integration_split

   USE mpas_derived_types
   USE mpas_pool_routines
   USE mpas_constants
   USE mpas_dmpar
   USE mpas_vector_reconstruction
   USE mpas_spline_interpolation
   USE mpas_timer
   USE mpas_threading
   USE mpas_timekeeping
   USE mpas_log
   USE mpas_io_units

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

   PUBLIC :: ocn_time_integrator_split, ocn_time_integration_split_init

   CHARACTER(LEN=*), PARAMETER :: subcycleGroupName = 'subcycleFields'
   CHARACTER(LEN=*), PARAMETER :: finalBtrGroupName = 'finalBtrFields'
   INTEGER :: nBtrSubcycles

   CONTAINS
   !
   !  ocn_time_integration_split
   !
   !> \brief MPAS ocean split explicit time integration scheme
   !> \author Mark Petersen, Doug Jacobsen, Todd Ringler
   !> \date   September 2011
   !> \details
   !> This routine integrates a master time step (dt) using a split explicit time 
   !  integrator.
   !
   SUBROUTINE ocn_time_integrator_split(domain, dt)
   !
   ! Advance model state forward in time by the specified time step using
   ! Split_Explicit timestepping scheme
   !
   ! Input: domain - current model state in time level 1 (e.g., time_levs(1) % state &
   ! % h(:,:)) plus mesh meta-data
   ! 
   ! Output: domain - upon exit, time level 2 (e.g., time_levs(2) % state % h(:,:)) 
   ! contains model state advanced forward in time by dt seconds
   ! 
      IMPLICIT NONE

      TYPE(domain_type), INTENT(INOUT) :: domain
      REAL(KIND=RKIND), INTENT(IN) :: dt
      
      ! Local variables
      INTEGER :: iCell, i,k,j, iEdge, cell1, cell2, split_explicit_step, split, &
                 eoe, oldBtrSubcycleTime, newBtrSubcycleTime, uPerpTime, BtrCorIter, &
                 stage1_tend_time, useVelocityCorrection, err, tsIter, &
                 edgeHaloComputeCounter, cellHaloComputeCounter
      INTEGER, DIMENSION(:), ALLOCATABLE :: n_bcl_iter
      LOGICAL :: activeTracersOnly ! if true only compute tendencies for active tracers 
      REAL(KIND=RKIND) :: normalThicknessFluxSum, thicknessSum, flux, sshEdge, hEdge1, &
                          CoriolisTerm, normalVelocityCorrection, temp, temp_h, coef, &
                          barotropicThicknessFlux_coeff, sshCell1, sshCell2
      REAL(KIND=RKIND), DIMENSION(:), POINTER :: btrvel_temp                    
      REAL(KIND=RKIND), DIMENSION(:,:), POINTER :: vertViscTopOfEdge, vertDiffTopOfCell
      REAL(KIND=RKIND), DIMENSION(:,:,:), POINTER :: tracersGroup
      REAL(KIND=RKIND), DIMENSION(:), ALLOCATABLE:: uTemp
      REAL(KIND=RKIND), DIMENSION(:,:), ALLOCATABLE:: tracersTemp
      TYPE(dm_info) :: dminfo
      TYPE(block_type), POINTER :: BLOCK
      TYPE(field1DReal), POINTER :: btrvel_tempField
      TYPE(mpas_pool_type), POINTER :: statePool, tracersPool, meshPool, &
                                       verticalMeshPool, diagnosticsPool, tendPool, &
                                       tracersTendPool, forcingPool, scratchPool, &
                                       swForcingPool

      ! Config options
      CHARACTER(LEN=StrKIND), POINTER :: config_time_integrator, &
                                         config_debug_split_explicit_filename, &
                                         config_land_ice_flux_mode
      INTEGER :: debugSplitExplicitUnit
      INTEGER, POINTER :: config_n_bcl_iter_mid, config_n_bcl_iter_beg, &
                          config_n_bcl_iter_end, config_n_ts_iter, &
                          config_btr_subcycle_loop_factor, config_n_btr_cor_iter, &
                          config_num_halos
      LOGICAL, POINTER :: config_use_standardGM, config_use_freq_filtered_thickness, &
                          config_btr_solve_SSH2, config_filter_btr_mode, &
                          config_vel_correction, config_prescribe_velocity, &
                          config_prescribe_thickness, config_use_cvmix_kpp, &
                          config_use_tracerGroup, &
                          config_compute_active_tracer_budgets, &
                          config_debug_split_explicit
      REAL(KIND=RKIND), POINTER :: config_mom_del4, config_btr_gam1_velWt1, &
                                   config_btr_gam2_SSHWt1, config_btr_gam3_velWt2

      ! Dimensions
      INTEGER :: nCells, nEdges
      INTEGER, POINTER :: nCellsPtr, nEdgesPtr, nVertLevels, num_tracersGroup, &
                          startIndex, endIndex, indexTemperature, indexSalinity, &
                          indexSurfaceVelocityZonal, indexSurfaceVelocityMeridional, &
                          indexSSHGradientZonal, indexSSHGradientMeridional
      INTEGER, DIMENSION(:), POINTER :: nCellsArray, nEdgesArray

      ! Mesh array pointers
      INTEGER, DIMENSION(:), POINTER :: maxLevelCell, maxLevelEdgeTop, nEdgesOnEdge, &
                                        nEdgesOnCell
      INTEGER, DIMENSION(:,:), POINTER :: cellsOnEdge, edgeMask, edgesOnEdge, &
                                          edgesOnCell, edgeSignOnCell

      REAL(KIND=RKIND), DIMENSION(:), POINTER &
      :: dcEdge, fEdge, bottomDepth, refBottomDepthTopOfCell, dvEdge, areaCell
      REAL(KIND=RKIND), DIMENSION(:,:), POINTER :: weightsOnEdge

      ! State Array Pointers
      REAL(KIND=RKIND), DIMENSION(:), POINTER &
      :: sshSubcycleCur, sshSubcycleNew, normalBarotropicVelocitySubcycleCur, &
         normalBarotropicVelocitySubcycleNew, sshCur, sshNew, &
         normalBarotropicVelocityCur, normalBarotropicVelocityNew
      
      REAL(KIND=RKIND), DIMENSION(:,:), POINTER &
      :: normalBaroclinicVelocityCur, normalBaroclinicVelocityNew, normalVelocityCur, &
         normalVelocityNew, layerThicknessCur, layerThicknessNew, &
         highFreqThicknessCur, highFreqThicknessNew, lowFreqDivergenceCur, &
         lowFreqDivergenceNew
      REAL(KIND=RKIND), DIMENSION(:,:,:), POINTER :: tracersGroupCur, tracersGroupNew

      ! Tend Array Pointers
      REAL(KIND=RKIND), DIMENSION(:), POINTER :: sshTend
      REAL(KIND=RKIND), DIMENSION(:,:), POINTER &
      :: highFreqThicknessTend, lowFreqDivergenceTend, normalVelocityTend, &
         layerThicknessTend
      REAL(KIND=RKIND), DIMENSION(:,:,:), POINTER :: tracersGroupTend, activeTracersTend

      ! Diagnostics Array Pointers
      REAL(KIND=RKIND), DIMENSION(:), POINTER &
      :: barotropicForcing, barotropicThicknessFlux, gradSSH, gradSSHX, gradSSHY, &
         gradSSHZ, gradSSHZonal, gradSSHMeridional
      REAL(KIND=RKIND), DIMENSION(:,:), POINTER &
      :: layerThicknessEdge, normalTransportVelocity, normalGMBolusVelocity, &
         vertAleTransportTop, velocityX, velocityY, velocityZ, velocityZonal, &
         velocityMeridional, surfaceVelocity, SSHGradient

      ! Diagnostics Field Pointers
      TYPE(field1DReal), POINTER :: barotropicThicknessFluxField, &
                                    boundaryLayerDepthField, effectiveDensityField
      TYPE(field2DReal), POINTER :: normalizedRelativeVorticityEdgeField, &
                                    divergenceField, relativeVorticityField
                                    
      ! Tracer tendencies brought in here to normalize by new layer thickness
      REAL(KIND=RKIND), DIMENSION(:,:), POINTER :: temperatureShortWaveTendency
      REAL(KIND=RKIND), DIMENSION(:,:,:), POINTER &
      :: activeTracerHorizontalAdvectionTendency, &
         activeTracerVerticalAdvectionTendency, activeTracerSurfaceFluxTendency, &
         activeTracerNonLocalTendency
        
      ! State/Tend Field Pointers
      TYPE(field1DReal), POINTER :: normalBarotropicVelocitySubcycleField, &      
                                    sshSubcycleField
      TYPE(field2DReal), POINTER :: highFreqThicknessField, lowFreqDivergenceField
      TYPE(field2DReal), POINTER :: normalBaroclinicVelocityField, layerThicknessField
      TYPE(field2DReal), POINTER :: normalVelocityField
      TYPE(field3DReal), POINTER :: tracersGroupField

      ! Tracer iterators
      CHARACTER(LEN=StrKIND) :: modifiedGroupName, configName
      INTEGER :: threadNum, temp_mask
      TYPE(mpas_pool_iterator_type) :: groupItr

      CALL mpas_timer_start("se timestep")

      CALL mpas_pool_get_config(domain % configs, 'config_n_bcl_iter_beg', &
                                config_n_bcl_iter_beg)
      CALL mpas_pool_get_config(domain % configs, 'config_n_bcl_iter_mid', &
                                config_n_bcl_iter_mid)
      CALL mpas_pool_get_config(domain % configs, 'config_n_bcl_iter_end', &
                                config_n_bcl_iter_end)
      CALL mpas_pool_get_config(domain % configs, 'config_n_ts_iter', config_n_ts_iter)
      CALL mpas_pool_get_config(domain % configs, 'config_btr_subcycle_loop_factor', &
                                config_btr_subcycle_loop_factor)
      CALL mpas_pool_get_config(domain % configs, 'config_btr_gam1_velWt1', &
                                config_btr_gam1_velWt1)
      CALL mpas_pool_get_config(domain % configs, 'config_btr_gam3_velWt2', &
                                config_btr_gam3_velWt2)
      CALL mpas_pool_get_config(domain % configs, 'config_btr_solve_SSH2', &
                                config_btr_solve_SSH2)
      CALL mpas_pool_get_config(domain % configs, 'config_n_btr_cor_iter', &
                                config_n_btr_cor_iter)
      CALL mpas_pool_get_config(domain % configs, 'config_btr_gam2_SSHWt1', &
                                config_btr_gam2_SSHWt1)
      CALL mpas_pool_get_config(domain % configs, 'config_filter_btr_mode', &
                                config_filter_btr_mode)
      CALL mpas_pool_get_config(domain % configs, 'config_mom_del4', config_mom_del4)
      CALL mpas_pool_get_config(domain % configs, &
                                'config_use_freq_filtered_thickness', &
                                config_use_freq_filtered_thickness)
      CALL mpas_pool_get_config(domain % configs, 'config_time_integrator', &
                                config_time_integrator)
      CALL mpas_pool_get_config(domain % configs, 'config_vel_correction', &
                                config_vel_correction)

      CALL mpas_pool_get_config(domain % configs, 'config_prescribe_velocity', &
                                config_prescribe_velocity)
      CALL mpas_pool_get_config(domain % configs, 'config_prescribe_thickness', &
                                config_prescribe_thickness)

      CALL mpas_pool_get_config(domain % configs, 'config_prescribe_velocity', &
                                config_prescribe_velocity)
      CALL mpas_pool_get_config(domain % configs, 'config_prescribe_thickness', &
                                config_prescribe_thickness)

      CALL mpas_pool_get_config(domain % configs, 'config_use_standardGM', &
                                config_use_standardGM)
      CALL mpas_pool_get_config(domain % configs, 'config_use_cvmix_kpp', &
                                config_use_cvmix_kpp)
      CALL mpas_pool_get_config(domain % configs, 'config_land_ice_flux_mode', &
                                config_land_ice_flux_mode)

      CALL mpas_pool_get_config(domain % configs, 'config_num_halos', config_num_halos)

      CALL mpas_pool_get_config(domain % configs, &
                                'config_compute_active_tracer_budgets', &
                                config_compute_active_tracer_budgets)
      CALL mpas_pool_get_config(domain % configs, 'config_debug_split_explicit', &
                                config_debug_split_explicit)
      CALL mpas_pool_get_config(domain % configs, &
                                'config_debug_split_explicit_filename', &
                                config_debug_split_explicit_filename)                                
                                
      ALLOCATE(n_bcl_iter(config_n_ts_iter))
      !
      ! Prep variables before first iteration
      !
      CALL mpas_timer_start("se prep")      
      block => domain % blocklist
      
      DO WHILE (ASSOCIATED(block))
      
         CALL mpas_pool_get_DIMENSION(block % dimensions, 'nCells', nCellsPtr)
         CALL mpas_pool_get_DIMENSION(block % dimensions, 'nEdges', nEdgesPtr)
         CALL mpas_pool_get_DIMENSION(block % dimensions, 'nCellsArray', nCellsArray)
         CALL mpas_pool_get_DIMENSION(block % dimensions, 'nEdgesArray', nEdgesArray)
         CALL mpas_pool_get_DIMENSION(block % dimensions, 'nVertLevels', nVertLevels)

         CALL mpas_pool_get_subpool(block % structs, 'state', statePool)
         CALL mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         CALL mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         CALL mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

         CALL mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', &
                                  normalBaroclinicVelocityCur, 1)
         CALL mpas_pool_get_array(statePool, 'normalBarotropicVelocity', &
                                  normalBarotropicVelocityCur, 1)
         CALL mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)

         CALL mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', &
                                  normalBaroclinicVelocityNew, 2)
         CALL mpas_pool_get_array(statePool, 'normalBarotropicVelocity', &
                                  normalBarotropicVelocityNew, 2)
         CALL mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)

         CALL mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
         CALL mpas_pool_get_array(statePool, 'ssh', sshNew, 2)

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

         CALL mpas_pool_get_array(diagnosticsPool, 'vertAleTransportTop', &
                                  vertAleTransportTop)

         CALL mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
         CALL mpas_pool_get_DIMENSION(tracersPool, 'index_salinity', indexSalinity)

         nCells = nCellsPtr
         nEdges = nEdgesPtr

         ! Initialize*variables that are used to compute baroclinic tendencies below.

         !$OMP DO SCHEDULE(runtime) PRIVATE(k)
         DO iEdge = 1, nEdges
            DO k = 1, nVertLevels !maxLevelEdgeTop % array(iEdge)

               ! The baroclinic velocity needs be recomputed at the beginning of a
               ! timestep because the implicit vertical mixing is conducted on the
               ! total u.  We keep normalBarotropicVelocity from the previous timestep.
               ! Note that normalBaroclinicVelocity may now include a barotropic 
               ! component, because the weights layerThickness have changed. That is OK,
               ! because the barotropicForcing variable subtracts out the barotropic 
               ! component from the baroclinic.
               
               normalBaroclinicVelocityCur(k,iEdge) &
               = normalVelocityCur(k,iEdge) - normalBarotropicVelocityCur(iEdge)

               normalVelocityNew(k,iEdge) = normalVelocityCur(k,iEdge)

               normalBaroclinicVelocityNew(k,iEdge) &
               = normalBaroclinicVelocityCur(k,iEdge)
               
            END DO
         END DO
         !$OMP END DO

         !$OMP DO SCHEDULE(runtime) PRIVATE(k)
         DO iCell = 1, nCells
            sshNew(iCell) = sshCur(iCell)
            DO k = 1, maxLevelCell(iCell)
               layerThicknessNew(k,iCell) = layerThicknessCur(k,iCell)
               ! Set vertAleTransportTop to zero for stage 1 velocity tendency, first 
               ! time through.
               vertAleTransportTop(k,iCell) = 0.0_RKIND
            END DO
         END DO
         !$OMP END DO

         CALL mpas_pool_begin_iteration(tracersPool)
         DO WHILE (mpas_pool_get_next_member(tracersPool, groupItr))
            IF (groupItr % memberType == MPAS_POOL_FIELD) THEN
               CALL mpas_pool_get_array(tracersPool, groupItr % memberName, &
                                        tracersGroupCur, 1)
               CALL mpas_pool_get_array(tracersPool, groupItr % memberName, &
                                        tracersGroupNew, 2)
               IF (ASSOCIATED(tracersGroupCur) .AND. ASSOCIATED(tracersGroupNew)) THEN
                  !$OMP DO SCHEDULE(runtime) PRIVATE(k)
                  DO iCell = 1, nCells
                     DO k = 1, maxLevelCell(iCell)
                        tracersGroupNew(:,k,iCell) = tracersGroupCur(:,k,iCell)
                     END DO
                  END DO
                  !$OMP END DO
               END IF
            END IF
         END DO

         IF (ASSOCIATED(highFreqThicknessNew)) THEN
            !$OMP DO SCHEDULE(runtime)
            DO iCell = 1, nCells
               highFreqThicknessNew(:, iCell) = highFreqThicknessCur(:, iCell)
            END DO
            !$OMP END DO
         END IF

         IF (ASSOCIATED(lowFreqDivergenceNew)) THEN
            !$OMP DO SCHEDULE(runtime)
            DO iCell = 1, nCells
               lowFreqDivergenceNew(:, iCell) = lowFreqDivergenceCur(:, iCell)
            END DO
            !$OMP END DO
         END IF

         block => block % next
         
      END DO

      CALL mpas_timer_stop("se prep")
      ! 
      ! Begin large iteration loop
      ! 
      n_bcl_iter = config_n_bcl_iter_mid
      n_bcl_iter(1) = config_n_bcl_iter_beg
      n_bcl_iter(config_n_ts_iter) = config_n_bcl_iter_end

      DO split_explicit_step = 1, config_n_ts_iter
      
         CALL mpas_timer_start('se loop')

         stage1_tend_time = MIN(split_explicit_step,2)

         CALL mpas_pool_get_subpool(domain % blocklist % structs, 'diagnostics', &
                                    diagnosticsPool)

         CALL mpas_threading_barrier()
         
         ! --- Update halos for diagnostic ocean boundayr layer depth
         IF (config_use_cvmix_kpp) THEN
            CALL mpas_timer_start("se halo diag obd")
            CALL mpas_dmpar_field_halo_exch(domain, 'boundaryLayerDepth')
            CALL mpas_timer_stop("se halo diag obd")
         END IF

         ! --- Update halos for diagnostic variables
         CALL mpas_timer_start("se halo diag")

         CALL mpas_dmpar_field_halo_exch(domain, 'normalizedRelativeVorticityEdge')
         IF (config_mom_del4 > 0.0_RKIND) THEN
           CALL mpas_dmpar_field_halo_exch(domain, 'divergence')
           CALL mpas_dmpar_field_halo_exch(domain, 'relativeVorticity')
         END IF
         CALL mpas_timer_stop("se halo diag")
         !
         ! Stage 1: Baroclinic velocity (3D) prediction, explicit with long timestep
         !
         IF (config_use_freq_filtered_thickness) THEN
         
            CALL mpas_timer_start("se freq-filtered-thick computations")            
            block => domain % blocklist
            
            DO WHILE (ASSOCIATED(block))
               CALL mpas_pool_get_subpool(block % structs, 'tend', tendPool)
               CALL mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
               CALL mpas_pool_get_subpool(block % structs, 'state', statepool)
               CALL mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
               CALL mpas_pool_get_subpool(block % structs, 'diagnostics', &
                                          diagnosticsPool)
               CALL mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
               CALL ocn_tend_freq_filtered_thickness(tendPool, statePool, &
                                                     diagnosticsPool, meshPool, &
                                                     stage1_tend_time)
               block => block % next
            END DO
            
            CALL mpas_timer_stop("se freq-filtered-thick computations")

            CALL mpas_threading_barrier()

            CALL mpas_timer_start("se freq-filtered-thick halo update")

            CALL mpas_dmpar_field_halo_exch(domain, 'tendHighFreqThickness')
            CALL mpas_dmpar_field_halo_exch(domain, 'tendLowFreqDivergence')

            CALL mpas_timer_stop("se freq-filtered-thick halo update")
            CALL mpas_threading_barrier()

            block => domain % blocklist
            
            DO WHILE (ASSOCIATED(block))
               CALL mpas_pool_get_DIMENSION(block % dimensions, 'nCells', nCellsPtr)
               CALL mpas_pool_get_DIMENSION(block % dimensions, 'nCellsArray', &
                                            nCellsArray)
               CALL mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
               CALL mpas_pool_get_subpool(block % structs, 'state', statePool)
               CALL mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
               CALL mpas_pool_get_subpool(block % structs, 'tend', tendPool)
               CALL mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
               CALL mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
               CALL mpas_pool_get_array(statePool, 'highFreqThickness', &
                                        highFreqThicknessCur, 1)
               CALL mpas_pool_get_array(statePool, 'highFreqThickness', &
                                        highFreqThicknessNew, 2)
               CALL mpas_pool_get_array(tendPool, 'highFreqThickness', &
                                        highFreqThicknessTend)
               nCells = nCellsPtr
               !$OMP DO SCHEDULE(runtime) PRIVATE(k)
               DO iCell = 1, nCells
                  DO k = 1, maxLevelCell(iCell)
                     ! This is h^{hf}_{n+1}
                     highFreqThicknessNew(k,iCell) = highFreqThicknessCur(k,iCell) &
                                                     + dt*highFreqThicknessTend(k,iCell)
                  END DO
               END DO
               !$OMP END DO
               block => block % next
            END DO

         END IF

         ! Compute velocity tendencies, T(u*,w*,p*)
         CALL mpas_timer_start("se bcl vel")
         CALL mpas_timer_start('se bcl vel tend')
         block => domain % blocklist
         
         DO WHILE (ASSOCIATED(block))
           CALL mpas_pool_get_subpool(block % structs, 'tend', tendPool)
           CALL mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
           CALL mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
           CALL mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)
           CALL mpas_pool_get_subpool(block % structs, 'state', statePool)
           CALL mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
           CALL mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
           CALL mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
           CALL mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
           CALL mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
           CALL mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, &
                                    stage1_tend_time)
           CALL mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
           CALL mpas_pool_get_array(statePool, 'highFreqThickness', &
                                    highFreqThicknessNew, 2)
           CALL mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', &
                                    layerThicknessEdge)
           CALL ocn_tend_vel(tendPool, statePool, forcingPool, diagnosticsPool, &
                             meshPool, scratchPool, stage1_tend_time, dt)
           block => block % next
         END DO
         CALL mpas_timer_stop('se bcl vel tend')
         !
         ! Begin baroclinic iterations on linear Coriolis term
         !
         DO j=1,n_bcl_iter(split_explicit_step)

            ! Use this G coefficient to avoid an if statement within the iEdge loop.
            IF (TRIM(config_time_integrator) == 'unsplit_explicit') THEN
               split = 0
            ELSE IF (TRIM(config_time_integrator) == 'split_explicit') THEN
               split = 1
            END IF

            CALL mpas_timer_start('bcl iters on linear Coriolis')            
            block => domain % blocklist
            
            DO WHILE (ASSOCIATED(block))
            
               CALL mpas_pool_get_DIMENSION(block % dimensions, 'nEdges', nEdgesPtr)
               CALL mpas_pool_get_DIMENSION(block % dimensions, 'nEdgesArray', &
                                            nEdgesArray)
               CALL mpas_pool_get_DIMENSION(block % dimensions, 'nVertLevels', &
                                            nVertLevels)

               CALL mpas_pool_get_subpool(block % structs, 'state', statePool)
               CALL mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
               CALL mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
               CALL mpas_pool_get_subpool(block % structs, 'tend', tendPool)
               CALL mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
               CALL mpas_pool_get_subpool(block % structs, 'diagnostics', &
                                          diagnosticsPool)

               CALL mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
               CALL mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
               CALL mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)

               CALL mpas_pool_get_array(statePool, 'normalVelocity', &
                                        normalVelocityNew, 2)
               CALL mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', &
                                        normalBaroclinicVelocityCur, 1)
               CALL mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', &
                                        normalBaroclinicVelocityNew, 2)
               CALL mpas_pool_get_array(statePool, 'ssh', sshNew, 2)

               CALL mpas_pool_get_array(tendPool, 'normalVelocity', normalVelocityTend)

               CALL mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', &
                                        layerThicknessEdge)
               CALL mpas_pool_get_array(diagnosticsPool, 'barotropicForcing', &
                                        barotropicForcing)

               ! Only need to loop over the 1 halo, since there is a halo exchange 
               ! immediately after this computation.
               nEdges = nEdgesArray(1)

               ! Put f*normalBaroclinicVelocity^{perp} in normalVelocityNew as a work 
               ! variable
               CALL ocn_fuperp(statePool, meshPool, 2)

               ALLOCATE(uTemp(nVertLevels))

               !$OMP DO SCHEDULE(runtime) &
               !$OMP PRIVATE(cell1, cell2, k, normalThicknessFluxSum, thicknessSum)
               DO iEdge = 1, nEdges
                  cell1 = cellsOnEdge(1,iEdge)
                  cell2 = cellsOnEdge(2,iEdge)

                  uTemp = 0.0_RKIND  
                  ! Could put this after with uTemp(maxleveledgetop+1:nvertlevels) = 0
                  DO k = 1, maxLevelEdgeTop(iEdge)
                     ! normalBaroclinicVelocityNew 
                     ! = normalBaroclinicVelocityOld 
                     !   + dt*(-f*normalBaroclinicVelocityPerp
                     !         + T(u*,w*,p*) + g*grad(SSH*))
                     ! Here uNew is a work variable containing 
                     ! -fEdge(iEdge)*normalBaroclinicVelocityPerp(k,iEdge)
                      uTemp(k) = normalBaroclinicVelocityCur(k,iEdge) &
                         + dt*(normalVelocityTend(k,iEdge) &
                         + normalVelocityNew(k,iEdge) &  
                           ! This is f*normalBaroclinicVelocity^{perp}
                         + split*gravity*(sshNew(cell2) - sshNew(cell1))/dcEdge(iEdge))
                  END DO

                  ! ThicknessSum is initialized outside the loop because on land 
                  ! boundaries.
                  ! maxLevelEdgeTop = 0, but I want to initialize thicknessSum with a
                  ! nonzero value to avoid a NaN.
                  normalThicknessFluxSum = layerThicknessEdge(1,iEdge)*uTemp(1)
                  thicknessSum  = layerThicknessEdge(1,iEdge)

                  DO k = 2, maxLevelEdgeTop(iEdge)
                     normalThicknessFluxSum &
                     = normalThicknessFluxSum + layerThicknessEdge(k,iEdge)*uTemp(k)
                     thicknessSum = thicknessSum + layerThicknessEdge(k,iEdge)
                  END DO
                  barotropicForcing(iEdge) &
                  = split*normalThicknessFluxSum/thicknessSum/dt
                  DO k = 1, maxLevelEdgeTop(iEdge)
                     ! These two steps are together here:
                     ! {\bf u}'_{k,n+1} = {\bf u}'_{k,n} - \Delta t {\overline {\bf G}}
                     ! {\bf u}'_{k,n+1/2} = \frac{1}{2}\left({\bf u}^{'}_{k,n} + {\bf 
                     ! u}'_{k,n+1}\right)
                     ! so that normalBaroclinicVelocityNew is at time n+1/2
                     normalBaroclinicVelocityNew(k,iEdge) &
                     = 0.5_RKIND*(normalBaroclinicVelocityCur(k,iEdge) + uTemp(k) &
                                  - dt*barotropicForcing(iEdge))
                  END DO
                  
               END DO ! iEdge
               !$OMP END DO

               DEALLOCATE(uTemp)

               block => block % next
               
            END DO

            CALL mpas_threading_barrier()

            CALL mpas_timer_start("se halo normalBaroclinicVelocity")
            CALL mpas_dmpar_field_halo_exch(domain, 'normalBaroclinicVelocity', &
                                            timeLevel=2)
            CALL mpas_timer_stop("se halo normalBaroclinicVelocity")

            CALL mpas_timer_stop('bcl iters on linear Coriolis')

         END DO  ! DO j=1,config_n_bcl_iter

         CALL mpas_timer_start('se halo barotropicForcing')
         CALL mpas_dmpar_field_halo_exch(domain, 'barotropicForcing')
         CALL mpas_timer_stop('se halo barotropicForcing')

         CALL mpas_timer_stop("se bcl vel")
         ! 
         ! End baroclinic iterations on linear Coriolis term

         ! Stage 2: Barotropic velocity (2D) prediction, explicitly subcycled
         !
         CALL mpas_timer_start("se btr vel")

         oldBtrSubcycleTime = 1
         newBtrSubcycleTime = 2
         
         IF (config_debug_split_explicit) THEN
            ! Open up temporary file
            CALL mpas_new_unit(debugSplitExplicitUnit)
            ! Note that it presumes we are running on a single core here ...
            OPEN (UNIT=debugSplitExplicitUnit, &
                  FILE=TRIM(config_debug_split_explicit_filename), FORM='formatted', &
                  STATUS='UNKNOWN')
         END IF 

         IF (TRIM(config_time_integrator) == 'unsplit_explicit') THEN

            CALL mpas_timer_start('btr vel ue')            
            block => domain % blocklist
            
            DO WHILE (ASSOCIATED(block))
            
               CALL mpas_pool_get_DIMENSION(block % dimensions, 'nEdges', nEdgesPtr)
               CALL mpas_pool_get_DIMENSION(block % dimensions, 'nEdgesArray', &
                                            nEdgesArray)
               CALL mpas_pool_get_DIMENSION(block % dimensions, 'nVertLevels', &
                                            nVertLevels)
               CALL mpas_pool_get_subpool(block % structs, 'state', statePool)
               CALL mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
               CALL mpas_pool_get_subpool(block % structs, 'diagnostics', &
                                          diagnosticsPool)
               CALL mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
               CALL mpas_pool_get_array(statePool, 'normalBarotropicVelocity', &
                                        normalBarotropicVelocityNew, 2)
               CALL mpas_pool_get_array(statePool, 'normalVelocity', &
                                        normalVelocityNew, 2)
               CALL mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', &
                                        normalBaroclinicVelocityNew, 2)
               CALL mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', &
                                        normalTransportVelocity)
               CALL mpas_pool_get_array(diagnosticsPool, 'normalGMBolusVelocity', &
                                        normalGMBolusVelocity)
               CALL mpas_pool_get_array(meshPool, 'edgeMask', edgeMask)

               nEdges = nEdgesPtr

               ! For Split_Explicit unsplit, simply set normalBarotropicVelocityNew = 0,
               ! normalBarotropicVelocitySubcycle = 0, and 
               ! uNew = normalBaroclinicVelocityNew

               !$OMP DO SCHEDULE(runtime) PRIVATE(k)
               DO iEdge = 1, nEdges
                  normalBarotropicVelocityNew(iEdge) = 0.0_RKIND
                  DO k = 1, nVertLevels
                     normalVelocityNew(k, iEdge) &
                     = normalBaroclinicVelocityNew(k, iEdge)
                     ! normalTransportVelocity &
                     ! = normalBaroclinicVelocity + normalGMBolusVelocity
                     ! This is u used in advective terms for layerThickness and tracers
                     ! in tendency calls in stage 3.
                     normalTransportVelocity(k,iEdge) &
                     = edgeMask(k,iEdge)*(normalBaroclinicVelocityNew(k,iEdge) &
                       + normalGMBolusVelocity(k,iEdge))
                  END DO
                  
               END DO  ! iEdge
               !$OMP END DO

               block => block % next
               
            END DO  ! block
            
            CALL mpas_timer_stop('btr vel ue')

         ELSE IF (TRIM(config_time_integrator) == 'split_explicit') THEN

            ! Initialize variables for barotropic subcycling
            CALL mpas_timer_start('btr vel se init')            
            block => domain % blocklist
            
            DO WHILE (ASSOCIATED(block))
            
               CALL mpas_pool_get_DIMENSION(block % dimensions, 'nCells', nCellsPtr)
               CALL mpas_pool_get_DIMENSION(block % dimensions, 'nEdges', nEdgesPtr)
               CALL mpas_pool_get_DIMENSION(block % dimensions, 'nCellsArray', &
                                            nCellsArray)
               CALL mpas_pool_get_DIMENSION(block % dimensions, 'nEdgesArray', &
                                            nEdgesArray)
               CALL mpas_pool_get_subpool(block % structs, 'diagnostics', &
                                          diagnosticsPool)
               CALL mpas_pool_get_subpool(block % structs, 'state', statePool)
               CALL mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

               CALL mpas_pool_get_array(diagnosticsPool, 'barotropicForcing', &
                                        barotropicForcing)
               CALL mpas_pool_get_array(diagnosticsPool, 'barotropicThicknessFlux', &
                                        barotropicThicknessFlux)
               CALL mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
               CALL mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleCur, &
                                        oldBtrSubcycleTime)
               CALL mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', &
                                        normalBarotropicVelocitySubcycleCur, &
                                        oldBtrSubcycleTime)
               CALL mpas_pool_get_array(statePool, 'normalBarotropicVelocity', &
                                        normalBarotropicVelocityCur, 1)
               CALL mpas_pool_get_array(statePool, 'normalBarotropicVelocity', &                  
                                        normalBarotropicVelocityNew, 2)

               nCells = nCellsPtr
               nEdges = nEdgesPtr

               IF (config_filter_btr_mode) THEN
                  !$OMP DO SCHEDULE(runtime)
                  DO iEdge = 1, nEdges
                     barotropicForcing(iEdge) = 0.0_RKIND
                  END DO
                  !$OMP END DO
               END IF

               !$OMP DO SCHEDULE(runtime)
               DO iCell = 1, nCells
                  ! sshSubcycleOld = sshOld
                  sshSubcycleCur(iCell) = sshCur(iCell)
               END DO
               !$OMP END DO

               !$OMP DO SCHEDULE(runtime)
               DO iEdge = 1, nEdges

                  ! normalBarotropicVelocitySubcycleOld = normalBarotropicVelocityOld
                  normalBarotropicVelocitySubcycleCur(iEdge) &
                  = normalBarotropicVelocityCur(iEdge)

                  ! normalBarotropicVelocityNew = BtrOld  
                  ! This is the first for the summation
                  normalBarotropicVelocityNew(iEdge) &
                  = normalBarotropicVelocityCur(iEdge)

                  ! barotropicThicknessFlux = 0
                  barotropicThicknessFlux(iEdge) = 0.0_RKIND
               END DO
               !$OMP END DO

               block => block % next
               
            END DO  ! block
            
            CALL mpas_timer_stop('btr vel se init')
            !
            ! Begin Barotropic subcycle loop
            !
            ! Allocate subcycled scratch fields before starting subcycle loop
            CALL mpas_pool_get_subpool(domain % blocklist % structs, 'scratch', &                  
                                       scratchPool)
            CALL mpas_pool_get_field(scratchPool, 'btrvel_temp', btrvel_tempField)
            CALL mpas_allocate_scratch_field(btrvel_tempField, .FALSE.)

            CALL mpas_threading_barrier()

            CALL mpas_timer_start('btr se subcycle loop')
            
            DO j = 1, nBtrSubcycles*config_btr_subcycle_loop_factor
                                     
               cellHaloComputeCounter = config_num_halos
               edgeHaloComputeCounter = config_num_halos + 1
               ! 
               ! Barotropic subcycle: VELOCITY PREDICTOR STEP
               ! 
               IF (config_btr_gam1_velWt1 > 1.0e-12_RKIND) THEN  
               ! Only do this part if it is needed in next SSH solve
                  uPerpTime = oldBtrSubcycleTime
                  block => domain % blocklist
                  
                  DO WHILE (ASSOCIATED(block))
                  
                     CALL mpas_pool_get_DIMENSION(block % dimensions, 'nEdges', &
                                                  nEdgesPtr)
                     CALL mpas_pool_get_DIMENSION(block % dimensions, 'nEdgesArray', &
                                                  nEdgesArray)
                     CALL mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
                     CALL mpas_pool_get_subpool(block % structs, 'state', statePool)
                     CALL mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
                     CALL mpas_pool_get_subpool(block % structs, 'diagnostics', &
                                                diagnosticsPool)
                     CALL mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
                     CALL mpas_pool_get_array(meshPool, 'nEdgesOnEdge', nEdgesOnEdge)
                     CALL mpas_pool_get_array(meshPool, 'edgesOnEdge', edgesOnEdge)
                     CALL mpas_pool_get_array(meshPool, 'weightsOnEdge', weightsOnEdge)
                     CALL mpas_pool_get_array(meshPool, 'fEdge', fEdge)
                     CALL mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
                     CALL mpas_pool_get_array(meshPool, 'edgeMask', edgeMask)
                     CALL mpas_pool_get_array(statePool, &
                                              'normalBarotropicVelocitySubcycle', &
                                              normalBarotropicVelocitySubcycleCur, &
                                              uPerpTime)
                     CALL mpas_pool_get_array(statePool, &
                                              'normalBarotropicVelocitySubcycle', &
                                              normalBarotropicVelocitySubcycleNew, &
                                              newBtrSubcycleTime)
                     CALL mpas_pool_get_array(statePool, 'sshSubcycle', &
                                              sshSubcycleCur, oldBtrSubcycleTime)
                     CALL mpas_pool_get_array(diagnosticsPool, 'barotropicForcing', &
                                              barotropicForcing)

                     nEdges = nEdgesPtr
                     nEdges = nEdgesArray(edgeHaloComputeCounter)

                     !$OMP DO SCHEDULE(runtime) &
                     !$OMP PRIVATE(cell1, cell2, CoriolisTerm, i, eoe)
                     DO iEdge = 1, nEdges

                        temp_mask = edgeMask(1, iEdge)

                        cell1 = cellsOnEdge(1,iEdge)
                        cell2 = cellsOnEdge(2,iEdge)

                        ! Compute the barotropic Coriolis term, -f*uPerp
                        CoriolisTerm = 0.0_RKIND
                        DO i = 1, nEdgesOnEdge(iEdge)
                           eoe = edgesOnEdge(i,iEdge)
                           CoriolisTerm &
                           = CoriolisTerm + weightsOnEdge(i,iEdge) &
                             *normalBarotropicVelocitySubcycleCur(eoe)*fEdge(eoe)
                        END DO

                        normalBarotropicVelocitySubcycleNew(iEdge) &
                        = temp_mask &
                          *(normalBarotropicVelocitySubcycleCur(iEdge) &
                            + dt/nBtrSubcycles*(CoriolisTerm - gravity &
                          *(sshSubcycleCur(cell2) - sshSubcycleCur(cell1)) &
                          /dcEdge(iEdge) + barotropicForcing(iEdge)))

                     END DO
                     !$OMP END DO

                     block => block % next
                     
                  END DO  ! block
                  
               END IF ! config_btr_gam1_velWt1 > 1.0e-12

               ! 1 Halo from edges is corrupted here, so reduce the edge halo layers by
               ! 1
               edgeHaloComputeCounter = edgeHaloComputeCounter - 1
               ! 
               ! Barotropic subcycle: SSH PREDICTOR STEP
               ! 
               block => domain % blocklist
               
               DO WHILE (ASSOCIATED(block))
               
                  CALL mpas_pool_get_DIMENSION(block % dimensions, 'nCells', nCellsPtr)
                  CALL mpas_pool_get_DIMENSION(block % dimensions, 'nEdges', nEdgesPtr)
                  CALL mpas_pool_get_DIMENSION(block % dimensions, 'nCellsArray', &
                                               nCellsArray)
                  CALL mpas_pool_get_DIMENSION(block % dimensions, 'nEdgesArray', &
                                               nEdgesArray)
                  CALL mpas_pool_get_subpool(block % structs, 'tend', tendPool)
                  CALL mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
                  CALL mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
                  CALL mpas_pool_get_subpool(block % structs, 'state', statePool)
                  CALL mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
                  CALL mpas_pool_get_subpool(block % structs, 'diagnostics', &
                                             diagnosticsPool)
                  CALL mpas_pool_get_array(tendPool, 'ssh', sshTend)
                  CALL mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
                  CALL mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
                  CALL mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
                  CALL mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
                  CALL mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
                  CALL mpas_pool_get_array(meshPool, 'refBottomDepthTopOfCell', &
                                           refBottomDepthTopOfCell)
                  CALL mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
                  CALL mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
                  CALL mpas_pool_get_array(meshPool, 'areaCell', areaCell)
                  CALL mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleCur, &
                                           oldBtrSubcycleTime)
                  CALL mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleNew, &
                                           newBtrSubcycleTime)
                  CALL mpas_pool_get_array(statePool, &
                                           'normalBarotropicVelocitySubcycle', &
                                           normalBarotropicVelocitySubcycleCur, &
                                           oldBtrSubcycleTime)
                  CALL mpas_pool_get_array(statePool, &
                                           'normalBarotropicVelocitySubcycle', &
                                           normalBarotropicVelocitySubcycleNew, &
                                           newBtrSubcycleTime)
                  CALL mpas_pool_get_array(diagnosticsPool, 'barotropicThicknessFlux', &
                                           barotropicThicknessFlux)
                  nCells = nCellsPtr
                  nEdges = nEdgesPtr
                  nCells = nCellsArray(cellHaloComputeCounter)
                  nEdges = nEdgesArray(edgeHaloComputeCounter)
                
                  ! Config_btr_gam1_velWt1 sets the forward weighting of velocity in the
                  ! SSH computation
                  ! config_btr_gam1_velWt1=  1     flux = normalBarotropicVelocityNew*H
                  ! config_btr_gam1_velWt1=0.5     
                  ! flux 
                  ! = 1/2*(normalBarotropicVelocityNew+normalBarotropicVelocityOld)*H
                  ! config_btr_gam1_velWt1 = 0     
                  ! flux = normalBarotropicVelocityOld*H

                  !$OMP DO SCHEDULE(runtime) &
                  !$OMP PRIVATE(i, iEdge, cell1, cell2, sshEdge, thicknessSum, flux)
                  DO iCell = 1, nCells
                     sshTend(iCell) = 0.0_RKIND
                     DO i = 1, nEdgesOnCell(iCell)
                        iEdge = edgesOnCell(i, iCell)   
                        cell1 = cellsOnEdge(1, iEdge)
                        cell2 = cellsOnEdge(2, iEdge)
                        sshEdge = 0.5_RKIND*(sshSubcycleCur(cell1) &
                                  + sshSubcycleCur(cell2))
                       ! method 0: orig, works only without pbc:
                       ! thicknessSum 
                       ! = sshEdge 
                       !   + refBottomDepthTopOfCell(maxLevelEdgeTop(iEdge)+1)
                       ! method 1, matches method 0 without pbcs, works with pbcs.
                       thicknessSum = sshEdge + MIN(bottomDepth(cell1), &
                                                    bottomDepth(cell2))
                       ! method 2: may be better than method 1.
                       ! Take average  of full thickness at two neighboring cells.
                       ! thicknessSum = sshEdge + 0.5*(bottomDepth(cell1) &
                       !                               + bottomDepth(cell2))
   
                       flux = ((1.0-config_btr_gam1_velWt1) &
                               *normalBarotropicVelocitySubcycleCur(iEdge) &
                               + config_btr_gam1_velWt1 &
                               *normalBarotropicVelocitySubcycleNew(iEdge)) &
                               *thicknessSum
   
                       sshTend(iCell) = sshTend(iCell) &
                                        + edgeSignOncell(i, iCell)*flux*dvEdge(iEdge)
                     END DO
                  ! SSHnew = SSHold + dt/J*(-div(Flux))
                  sshSubcycleNew(iCell) = sshSubcycleCur(iCell) &
                                          + dt/nBtrSubcycles*sshTend(iCell) &
                                            /areaCell(iCell)
                  END DO
                  !$OMP END DO

                  IF (config_debug_split_explicit) THEN
                     ! Add in write statement for debugging to print to a standard text 
                     ! file (get single strip out) 4.000
                     ! Add in a new line each time for the write statement
                     WRITE (debugSplitExplicitUnit,*) sshSubcycleNew(:nCells:4)  
                     ! 4 hard coded for an example for brevity of case
                     FLUSH (debugSplitExplicitUnit)
                  END IF                   
                  ! Asarje: Changed to avoid redundant computations when 
                  ! config_btr_solve_SSH2 is true
                  IF (config_btr_solve_SSH2) THEN
                     ! If config_btr_solve_SSH2=.TRUE.,
                     ! then do NOT accumulate barotropicThicknessFlux in this SSH 
                     ! predictor section, because it will be accumulated in the SSH 
                     ! corrector section.
                     barotropicThicknessFlux_coeff = 0.0_RKIND
                     ! Nothing else to do
                  ELSE
                  ! Otherwise, do accumulate barotropicThicknessFlux in this SSH 
                  ! predictor section
                  barotropicThicknessFlux_coeff = 1.0_RKIND
                  !$OMP DO SCHEDULE(runtime)
                  DO iEdge = 1, nEdges
                     cell1 = cellsOnEdge(1,iEdge)
                     cell2 = cellsOnEdge(2,iEdge)
                     sshEdge = 0.5_RKIND*(sshSubcycleCur(cell1) + sshSubcycleCur(cell2))
                     ! method 1, matches method 0 without pbcs, works with pbcs.
                     thicknessSum = sshEdge + MIN(bottomDepth(cell1), &
                                                  bottomDepth(cell2))
                     flux = ((1.0-config_btr_gam1_velWt1) &
                             *normalBarotropicVelocitySubcycleCur(iEdge) &
                             + config_btr_gam1_velWt1 &
                               *normalBarotropicVelocitySubcycleNew(iEdge)) &
                            *thicknessSum
                     barotropicThicknessFlux(iEdge) = barotropicThicknessFlux(iEdge) &
                                                      + flux
                  END DO
                  !$OMP END DO
               END IF
                
               block => block % next      
                         
            END DO  ! block

            ! 1 cell halo layer is now corrupted, so remove one from computing on.
            cellHaloComputeCounter = cellHaloComputeCounter - 1
            !
            ! Barotropic subcycle: VELOCITY CORRECTOR STEP
            !
            ! 1 edge halo is already corrupted from the predictor step.
            DO BtrCorIter = 1, config_n_btr_cor_iter
              
               uPerpTime = newBtrSubcycleTime               
               block => domain % blocklist
               
               DO WHILE (ASSOCIATED(block))
                  CALL mpas_pool_get_DIMENSION(block % dimensions, 'nEdges', nEdgesPtr)
                  CALL mpas_pool_get_DIMENSION(block % dimensions, 'nEdgesArray', &
                                               nEdgesArray)
                  CALL mpas_pool_get_subpool(block % structs, 'state', statePool)
                  CALL mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
                  CALL mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
                  CALL mpas_pool_get_subpool(block % structs, 'diagnostics', &
                                             diagnosticsPool)
                  CALL mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
                  CALL mpas_pool_get_array(statePool, &
                                           'normalBarotropicVelocitySubcycle', &
                                           normalBarotropicVelocitySubcycleCur, &
                                           oldBtrSubcycleTime)
                  CALL mpas_pool_get_array(statePool, &
                                           'normalBarotropicVelocitySubcycle', &
                                           normalBarotropicVelocitySubcycleNew, &
                                           newBtrSubcycleTime)
                  CALL mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleCur, &
                                           oldBtrSubcycleTime)
                  CALL mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleNew, &
                                           newBtrSubcycleTime)
                  CALL mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
                  CALL mpas_pool_get_array(meshPool, 'nEdgesOnEdge', nEdgesOnEdge)
                  CALL mpas_pool_get_array(meshPool, 'edgesOnEdge', edgesOnEdge)
                  CALL mpas_pool_get_array(meshPool, 'weightsOnEdge', weightsOnEdge)
                  CALL mpas_pool_get_array(meshPool, 'fEdge', fEdge)
                  CALL mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
                  CALL mpas_pool_get_array(meshPool, 'edgeMask', edgeMask)
                  CALL mpas_pool_get_array(diagnosticsPool, 'barotropicForcing', &
                                           barotropicForcing)
                  CALL mpas_pool_get_field(scratchPool, 'btrvel_temp', &
                                           btrvel_tempField)
                  btrvel_temp => btrvel_tempField % array

                  ! Need to initialize btr_vel_temp over the one more halo than we're 
                  ! computing over
                  nEdges = nEdgesPtr
                  nEdges = nEdgesArray(MIN(edgeHaloComputeCounter + 1, &
                                           config_num_halos + 1))

                  !$OMP DO SCHEDULE(runtime)
                  DO iEdge = 1, nEdges+1
                     btrvel_temp(iEdge) = normalBarotropicVelocitySubcycleNew(iEdge)
                  END DO
                  !$OMP END DO

                  nEdges = nEdgesArray(edgeHaloComputeCounter)

                  !$OMP DO SCHEDULE(runtime) &
                  !$OMP PRIVATE(cell1, cell2, eoe, CoriolisTerm, i, sshCell1, sshCell2)
                  DO iEdge = 1, nEdges

                     ! Asarje: Added to avoid redundant computations based on mask
                     temp_mask = edgeMask(1,iEdge)

                     cell1 = cellsOnEdge(1,iEdge)
                     cell2 = cellsOnEdge(2,iEdge)

                     ! Compute the barotropic Coriolis term, -f*uPerp
                     CoriolisTerm = 0.0_RKIND
                     DO i = 1, nEdgesOnEdge(iEdge)
                        eoe = edgesOnEdge(i,iEdge)
                        CoriolisTerm = CoriolisTerm &
                                       + weightsOnEdge(i,iEdge)*btrvel_temp(eoe) &
                                         *fEdge(eoe)
                     END DO

                     ! In this final solve for velocity, SSH is a linear
                     ! combination of SSHold and SSHnew.
                     sshCell1 = (1 - config_btr_gam2_SSHWt1)*sshSubcycleCur(cell1) &
                                + config_btr_gam2_SSHWt1*sshSubcycleNew(cell1)
                     sshCell2 = (1 - config_btr_gam2_SSHWt1)*sshSubcycleCur(cell2) &
                                + config_btr_gam2_SSHWt1*sshSubcycleNew(cell2)

                     ! normalBarotropicVelocityNew 
                     ! = normalBarotropicVelocityOld 
                     !   + dt/J*(-f*normalBarotropicVelocityoldPerp - g*grad(SSH) + G)
                     normalBarotropicVelocitySubcycleNew(iEdge) &
                     = temp_mask*(normalBarotropicVelocitySubcycleCur(iEdge) &
                       + dt/nBtrSubcycles &
                         *(CoriolisTerm - gravity*(sshCell2 - sshCell1)/dcEdge(iEdge) &
                           + barotropicForcing(iEdge)))

                  END DO
                  !$OMP END DO

                  block => block % next
                   
               END DO  ! block

               ! mrp 170503 This is the original line.  Go back to this once I verify 
               ! that halo regions are truly three cells wide.
               ! IF (edgeHaloComputeCounter == 1 .AND. &
               !     BtrCorIter /= config_n_btr_cor_iter) THEN
               ! replaced with this, to force a halo update on the first iteration:
               IF (edgeHaloComputeCounter == 1 .OR. BtrCorIter == 1) THEN
                  edgeHaloComputeCounter = config_num_halos + 1
                  CALL mpas_timer_start("se halo btr vel corr")
                  CALL mpas_dmpar_field_halo_exch(domain, &
                                                  'normalBarotropicVelocitySubcycle', &
                                                  timeLevel=newBtrSubcycleTime)
                  CALL mpas_timer_stop("se halo btr vel corr")
               ELSE
                  edgeHaloComputeCounter = edgeHaloComputeCounter - 1
               END IF

            END DO ! Do BtrCorIter=1,config_n_btr_cor_iter
            ! 
            ! Barotropic subcycle: SSH CORRECTOR STEP
            ! 
            IF (config_btr_solve_SSH2) THEN
            
               block => domain % blocklist
               
               DO WHILE (ASSOCIATED(block))
               
                  CALL mpas_pool_get_DIMENSION(block % dimensions, 'nCells', nCellsPtr)
                  CALL mpas_pool_get_DIMENSION(block % dimensions, 'nEdges', nEdgesPtr)
                  CALL mpas_pool_get_DIMENSION(block % dimensions, 'nCellsArray', &
                                               nCellsArray)
                  CALL mpas_pool_get_DIMENSION(block % dimensions, 'nEdgesArray', &
                                               nEdgesArray)
                  CALL mpas_pool_get_subpool(block % structs, 'tend', tendPool)
                  CALL mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
                  CALL mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
                  CALL mpas_pool_get_subpool(block % structs, 'state', statePool)
                  CALL mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
                  CALL mpas_pool_get_subpool(block % structs, 'diagnostics', &
                                             diagnosticsPool)
                  CALL mpas_pool_get_array(tendPool, 'ssh', sshTend)
                  CALL mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
                  CALL mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
                  CALL mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
                  CALL mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
                  CALL mpas_pool_get_array(meshPool, 'refBottomDepthTopOfCell', &
                                           refBottomDepthTopOfCell)
                  CALL mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
                  CALL mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
                  CALL mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
                  CALL mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleCur, &
                                           oldBtrSubcycleTime)
                  CALL mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleNew, &
                                           newBtrSubcycleTime)
                  CALL mpas_pool_get_array(statePool, &
                                           'normalBarotropicVelocitySubcycle', &
                                           normalBarotropicVelocitySubcycleCur, &
                                           oldBtrSubcycleTime)
                  CALL mpas_pool_get_array(statePool, &
                                           'normalBarotropicVelocitySubcycle', &
                                           normalBarotropicVelocitySubcycleNew, &
                                           newBtrSubcycleTime)
                  CALL mpas_pool_get_array(diagnosticsPool, 'barotropicThicknessFlux', &
                                           barotropicThicknessFlux)
                  nCells = nCellsPtr
                  nEdges = nEdgesPtr
                  nCells = nCellsArray(cellHaloComputeCounter)
                  nEdges = nEdgesArray(edgeHaloComputeCounter)
                  ! config_btr_gam3_velWt2 sets the forward weighting of velocity in the 
                  ! SSH computation
                  ! config_btr_gam3_velWt2 = 1     
                  ! flux = normalBarotropicVelocityNew*H
                  ! config_btr_gam3_velWt2 = 0.5     
                  ! flux = 1/2*(normalBarotropicVelocityNew 
                  !             + normalBarotropicVelocityOld)*H
                  ! config_btr_gam3_velWt2=  0     
                  ! flux = normalBarotropicVelocityOld*H

                  !$OMP DO SCHEDULE(runtime) &
                  !$OMP PRIVATE(i, iEdge, cell1, cell2, sshCell1, sshCell2, sshEdge, &
                  !$OMP&        thicknessSum, flux)
                  DO iCell = 1, nCells
                     sshTend(iCell) = 0.0_RKIND
                     DO i = 1, nEdgesOnCell(iCell)
                        iEdge = edgesOnCell(i, iCell)

                        cell1 = cellsOnEdge(1,iEdge)
                        cell2 = cellsOnEdge(2,iEdge)

                        ! SSH is a linear combination of SSHold and SSHnew.
                        sshCell1 = (1-config_btr_gam2_SSHWt1)* sshSubcycleCur(cell1) &
                                   +   config_btr_gam2_SSHWt1*sshSubcycleNew(cell1)
                        sshCell2 = (1-config_btr_gam2_SSHWt1)* sshSubcycleCur(cell2) &
                                 +   config_btr_gam2_SSHWt1*sshSubcycleNew(cell2)

                        sshEdge = 0.5_RKIND*(sshCell1 + sshCell2)

                        ! method 0: orig, works only without pbc:
                        ! thicknessSum 
                        ! = sshEdge + refBottomDepthTopOfCell(maxLevelEdgeTop(iEdge)+1)

                        ! method 1, matches method 0 without pbcs, works with pbcs.
                        thicknessSum = sshEdge + MIN(bottomDepth(cell1), &
                                                     bottomDepth(cell2))

                        ! method 2: may be better than method 1.
                        ! take average  of full thickness at two neighboring cells
                        ! thicknessSum = sshEdge + 0.5 *(bottomDepth(cell1) 
                        !                                + bottomDepth (cell2))
                        flux = ((1.0-config_btr_gam3_velWt2) &
                                *normalBarotropicVelocitySubcycleCur(iEdge) &
                                + config_btr_gam3_velWt2 &
                                  *normalBarotropicVelocitySubcycleNew(iEdge)) &
                               *thicknessSum

                        sshTend(iCell) = sshTend(iCell) + edgeSignOnCell(i, iCell) &
                                                          *flux*dvEdge(iEdge)
                     END DO
                     ! SSHnew = SSHold + dt/J*(-div(Flux))
                     sshSubcycleNew(iCell) = sshSubcycleCur(iCell) &
                                             + dt/nBtrSubcycles*sshTend(iCell) &
                                               /areaCell(iCell)
                  END DO
                  !$OMP END DO

                  !$OMP DO SCHEDULE(runtime) &
                  !$OMP PRIVATE(cell1, cell2, sshCell1, sshCell2, sshEdge, &
                  !$OMP&        thicknessSum, flux)
                  DO iEdge = 1, nEdges
                     cell1 = cellsOnEdge(1,iEdge)
                     cell2 = cellsOnEdge(2,iEdge)

                     ! SSH is a linear combination of SSHold and SSHnew.
                     sshCell1 = (1 - config_btr_gam2_SSHWt1)*sshSubcycleCur(cell1) &
                                 + config_btr_gam2_SSHWt1*sshSubcycleNew(cell1)
                     sshCell2 = (1-config_btr_gam2_SSHWt1)* sshSubcycleCur(cell2) &
                                 + config_btr_gam2_SSHWt1*sshSubcycleNew(cell2)
                     sshEdge = 0.5_RKIND*(sshCell1 + sshCell2)

                     ! method 0: orig, works only without pbc:
                     ! thicknessSum 
                     ! = sshEdge + refBottomDepthTopOfCell(maxLevelEdgeTop(iEdge)+1)

                     ! method 1, matches method 0 without pbcs, works with pbcs.
                     thicknessSum = sshEdge + MIN(bottomDepth(cell1), &
                                                   bottomDepth(cell2))

                     ! method 2, better, I think.
                     ! take average  of full thickness at two neighboring cells
                     ! thicknessSum = sshEdge + 0.5*(bottomDepth(cell1) 
                     !                               + bottomDepth(cell2))

                     flux = ((1.0-config_btr_gam3_velWt2) &
                             *normalBarotropicVelocitySubcycleCur(iEdge) &
                             + config_btr_gam3_velWt2 &
                               *normalBarotropicVelocitySubcycleNew(iEdge))*thicknessSum

                     barotropicThicknessFlux(iEdge) = barotropicThicknessFlux(iEdge) &
                                                      + flux
                  END DO
                  !$OMP END DO

                  block => block % next
                  
               END DO ! block
               
            END IF ! config_btr_solve_SSH2
            
            ! Boundary update on SSHnew
            CALL mpas_timer_start("se halo subcycle")
            CALL mpas_dmpar_exch_group_create(domain, subcycleGroupName)
            CALL mpas_dmpar_exch_group_add_field(domain, subcycleGroupName, &
                                                 'sshSubcycle', &
                                                 timeLevel=newBtrSubcycleTime)
            CALL mpas_dmpar_exch_group_add_field(domain, subcycleGroupName, &
                                                 'normalBarotropicVelocitySubcycle', &
                                                 timeLevel=newBtrSubcycleTime)

            CALL mpas_threading_barrier()
            CALL mpas_dmpar_exch_group_full_halo_exch(domain, subcycleGroupName)

            CALL mpas_dmpar_exch_group_destroy(domain, subcycleGroupName)
            CALL mpas_timer_stop("se halo subcycle")

            ! Reset the halo counters
            edgeHaloComputeCounter = config_num_halos + 1
            cellHaloComputeCounter = config_num_halos            
            ! 
            ! Barotropic subcycle: Accumulate running sums, advance timestep pointers
            ! 
            block => domain % blocklist
            
            DO WHILE (ASSOCIATED(block))
            
               CALL mpas_pool_get_DIMENSION(block % dimensions, 'nEdges', nEdgesPtr)
               CALL mpas_pool_get_DIMENSION(block % dimensions, 'nEdgesArray', &
                                            nEdgesArray)

               CALL mpas_pool_get_subpool(block % structs, 'state', statePool)
               CALL mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

               CALL mpas_pool_get_array(statePool, 'normalBarotropicVelocity', &
                                        normalBarotropicVelocityNew, 2)
               CALL mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', &
                                        normalBarotropicVelocitySubcycleNew, &
                                        newBtrSubcycleTime)

               ! normalBarotropicVelocityNew 
               ! = normalBarotropicVelocityNew + normalBarotropicVelocitySubcycleNEW
               ! This accumulates the sum.
               ! If the Barotropic Coriolis iteration is limited to one, this could
               ! be merged with the above code.

               nEdges = nEdgesPtr

               !$OMP DO SCHEDULE(runtime)
               DO iEdge = 1, nEdges
                  normalBarotropicVelocityNew(iEdge) &
                  = normalBarotropicVelocityNew(iEdge) &
                    + normalBarotropicVelocitySubcycleNew(iEdge)
               END DO ! iEdge
               !$OMP END DO

               block => block % next
               
            END DO  ! block

            ! Advance time pointers
            oldBtrSubcycleTime = MOD(oldBtrSubcycleTime,2)+1
            newBtrSubcycleTime = MOD(newBtrSubcycleTime,2)+1

            END DO ! j = 1, nBtrSubcycles
            CALL mpas_timer_stop('btr se subcycle loop')

            CALL mpas_threading_barrier()
            CALL mpas_pool_get_subpool(domain % blocklist % structs, 'scratch', &
                                       scratchPool)
            CALL mpas_pool_get_field(scratchPool, 'btrvel_temp', btrvel_tempField)
            CALL mpas_deallocate_scratch_field(btrvel_tempField, .FALSE.)
            ! 
            ! END Barotropic subcycle loop
            ! 
            ! Normalize Barotropic subcycle sums: ssh, normalBarotropicVelocity, and F
            CALL mpas_timer_start('btr se norm')
            block => domain % blocklist
            
            DO WHILE (ASSOCIATED(block))
            
               CALL mpas_pool_get_DIMENSION(block % dimensions, 'nEdges', nEdgesPtr)
               CALL mpas_pool_get_DIMENSION(block % dimensions, 'nEdgesArray', &
                                            nEdgesArray)

               CALL mpas_pool_get_subpool(block % structs, 'state', statePool)
               CALL mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
               CALL mpas_pool_get_subpool(block % structs, 'diagnostics', &
                                          diagnosticsPool)

               CALL mpas_pool_get_array(statePool, 'normalBarotropicVelocity', &
                                        normalBarotropicVelocityNew, 2)

               CALL mpas_pool_get_array(diagnosticsPool, 'barotropicThicknessFlux', &
                                        barotropicThicknessFlux)

               nEdges = nEdgesPtr

               nEdges = nEdgesArray(1)

               !$OMP DO SCHEDULE(runtime)
               DO iEdge = 1, nEdges
                  barotropicThicknessFlux(iEdge) &
                  = barotropicThicknessFlux(iEdge) &
                    /(nBtrSubcycles*config_btr_subcycle_loop_factor)

                  normalBarotropicVelocityNew(iEdge) &
                  = normalBarotropicVelocityNew(iEdge) &
                    /(nBtrSubcycles*config_btr_subcycle_loop_factor + 1)
               END DO
               !$OMP END DO

               block => block % next
               
            END DO  ! block
            
            CALL mpas_timer_stop('btr se norm')
            CALL mpas_threading_barrier()

            ! Boundary update on F
            CALL mpas_timer_start("se halo F and btr vel")
            CALL mpas_dmpar_exch_group_create(domain, finalBtrGroupName)
            CALL mpas_dmpar_exch_group_add_field(domain, finalBtrGroupName, &
                                                 'barotropicThicknessFlux')
            CALL mpas_dmpar_exch_group_add_field(domain, finalBtrGroupName, &
                                                 'normalBarotropicVelocity', &
                                                 timeLevel=2)
            CALL mpas_threading_barrier()
            CALL mpas_dmpar_exch_group_full_halo_exch(domain, finalBtrGroupName)
            CALL mpas_dmpar_exch_group_destroy(domain, finalBtrGroupName)
            CALL mpas_timer_stop("se halo F and btr vel")
            CALL mpas_threading_barrier()

            ! Check that you can compute SSH using the total sum or the individual 
            ! increments over the barotropic subcycles.
            ! Efficiency: This next block of code is really a check for debugging, and 
            ! can be removed later.
            CALL mpas_timer_start('btr se ssh verif')
            block => domain % blocklist
            
            DO WHILE (ASSOCIATED(block))
            
               CALL mpas_pool_get_DIMENSION(block % dimensions, 'nEdges', nEdgesPtr)
               CALL mpas_pool_get_DIMENSION(block % dimensions, 'nEdgesArray', &
                                            nEdgesArray)
               CALL mpas_pool_get_DIMENSION(block % dimensions, 'nVertLevels', &
                                            nVertLevels)
               CALL mpas_pool_get_subpool(block % structs, 'state', statePool)
               CALL mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
               CALL mpas_pool_get_subpool(block % structs, 'diagnostics', &
                                          diagnosticsPool)
               CALL mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
               CALL mpas_pool_get_array(statePool, 'normalBarotropicVelocity', &
                                        normalBarotropicVelocityNew, 2)
               CALL mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', &
                                        normalBaroclinicVelocityNew, 2)
               CALL mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', &
                                        normalTransportVelocity)
               CALL mpas_pool_get_array(diagnosticsPool, 'normalGMBolusVelocity', &
                                        normalGMBolusVelocity)
               CALL mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', &
                                        layerThicknessEdge)
               CALL mpas_pool_get_array(diagnosticsPool, 'barotropicThicknessFlux', &
                                        barotropicThicknessFlux)
               CALL mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
               CALL mpas_pool_get_array(meshPool, 'edgeMask', edgeMask)

               nEdges = nEdgesPtr
               nEdges = nEdgesArray(config_num_halos)

               ALLOCATE(uTemp(nVertLevels))

               ! Correction velocity normalVelocityCorrection = (Flux - Sum(h u*))/H
               ! or, for the full latex version:
               ! {\bf u}^{corr} = \left({\overline{\bf F}}
               ! - \sum_{k=1}^{N^{edge}} h_{k,*}^{edge} {\bf u}_k^{avg} \right)
               ! \left/ \sum_{k=1}^{N^{edge}} h_{k,*}^{edge} \right.

               IF (config_vel_correction) THEN
                  useVelocityCorrection = 1
               ELSE
                  useVelocityCorrection = 0
               END IF

               !$OMP DO SCHEDULE(runtime) &
               !$OMP PRIVATE(k, normalThicknessFluxSum, thicknessSum, &
               !$OMP&        normalVelocityCorrection)
               DO iEdge = 1, nEdges

                  ! velocity for normalVelocityCorrectionection is 
                  ! normalBarotropicVelocity + normalBaroclinicVelocity + uBolus
                  uTemp(:) = normalBarotropicVelocityNew(iEdge) &
                             + normalBaroclinicVelocityNew(:,iEdge) &
                             + normalGMBolusVelocity(:,iEdge)

                  ! thicknessSum is initialized outside the loop because on land 
                  ! boundaries maxLevelEdgeTop=0, but I want to initialize thicknessSum
                  ! with a nonzero value to avoid a NaN.
                  normalThicknessFluxSum = layerThicknessEdge(1,iEdge)*uTemp(1)
                  thicknessSum  = layerThicknessEdge(1,iEdge)

                  DO k = 2, maxLevelEdgeTop(iEdge)
                     normalThicknessFluxSum &
                     = normalThicknessFluxSum + layerThicknessEdge(k,iEdge)*uTemp(k)
                     thicknessSum  =  thicknessSum + layerThicknessEdge(k,iEdge)
                  END DO

                  normalVelocityCorrection &
                  = useVelocityCorrection*((barotropicThicknessFlux(iEdge) &
                                           - normalThicknessFluxSum)/thicknessSum)

                  DO k = 1, nVertLevels

                     ! normalTransportVelocity = normalBarotropicVelocity 
                     !                           + normalBaroclinicVelocity 
                     !                           + normalGMBolusVelocity
                     !                           + normalVelocityCorrection
                     ! This is u used in advective terms for layerThickness and tracers
                     ! in tendency calls in stage 3.
                     ! mrp note: in QC version, there is an IF (config_use_standardGM) 
                     ! on adding normalGMBolusVelocity. I think it is not needed because 
                     ! normalGMBolusVelocity = 0 when GM not on.
                     
                     normalTransportVelocity(k,iEdge) &
                     = edgeMask(k,iEdge)*(normalBarotropicVelocityNew(iEdge) &
                                          + normalBaroclinicVelocityNew(k,iEdge) &
                                          + normalGMBolusVelocity(k,iEdge) &
                                          + normalVelocityCorrection)
                  END DO

               END DO ! iEdge
               !$OMP END DO

               DEALLOCATE(uTemp)

               block => block % next
               
            END DO ! block
            CALL mpas_timer_stop('btr se ssh verif')

         END IF ! split_explicit
         
         IF (config_debug_split_explicit) THEN
            ! Close file
            FLUSH (debugSplitExplicitUnit)
            CLOSE (debugSplitExplicitUnit)
            CALL mpas_release_unit(debugSplitExplicitUnit)
         END IF

         CALL mpas_timer_stop("se btr vel")
         !
         ! Stage 3: Tracer, density, pressure, vertical velocity prediction
         !
         ! Only compute tendencies for active tracers on last large iteration
         IF (split_explicit_step < config_n_ts_iter) THEN
            activeTracersOnly = .TRUE.
         ELSE
            activeTracersOnly = .FALSE.
         END IF

         ! Thickness tendency computations and thickness halo updates are completed 
         ! before tracer tendency computations to allow monotonic advection.
         CALL mpas_timer_start('se thick tend')
         block => domain % blocklist
         
         DO WHILE (ASSOCIATED(block))       
            CALL mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            CALL mpas_pool_get_subpool(block % structs, 'verticalMesh', &
                                       verticalMeshPool)
            CALL mpas_pool_get_subpool(block % structs, 'state', statePool)
            CALL mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
            CALL mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
            CALL mpas_pool_get_subpool(block % structs, 'tend', tendPool)
            CALL mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
            CALL mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
            CALL mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
            CALL mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
            CALL mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
            CALL mpas_pool_get_array(statePool, 'highFreqThickness', &
                                     highFreqThicknessNew, 2)
            CALL mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', &
                                     layerThicknessEdge)
            CALL mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', &
                                     normalTransportVelocity)
            CALL mpas_pool_get_array(diagnosticsPool, 'vertAleTransportTop', &
                                     vertAleTransportTop)
            ! Compute vertAleTransportTop. Use normalTransportVelocity for advection of 
            ! layerThickness and tracers. Use time level 1 values of layerThickness and 
            ! layerThicknessEdge because layerThickness has not yet been computed for 
            ! time level 2.
            CALL mpas_timer_start('thick vert trans vel top')
            IF (ASSOCIATED(highFreqThicknessNew)) THEN
               CALL ocn_vert_transport_velocity_top &
                    (meshPool, verticalMeshPool, scratchPool, layerThicknessCur, &
                     layerThicknessEdge, normalTransportVelocity, sshCur, dt, &
                     vertAleTransportTop, err, highFreqThicknessNew)
            ELSE
               CALL ocn_vert_transport_velocity_top &
                    (meshPool, verticalMeshPool, scratchPool, layerThicknessCur, &
                     layerThicknessEdge, normalTransportVelocity, sshCur, dt, &
                     vertAleTransportTop, err)
            END IF
            CALL mpas_timer_stop('thick vert trans vel top')

            CALL ocn_tend_thick(tendPool, forcingPool, diagnosticsPool, meshPool)

            block => block % next
            
         END DO
         
         CALL mpas_timer_stop('se thick tend')
         CALL mpas_threading_barrier()
         ! update halo for thickness tendencies
         CALL mpas_timer_start("se halo thickness")
         CALL mpas_dmpar_field_halo_exch(domain, 'tendLayerThickness')
         CALL mpas_timer_stop("se halo thickness")
         CALL mpas_threading_barrier()
         CALL mpas_timer_start('se tracer tend', .FALSE.)
         block => domain % blocklist
         
         DO WHILE (ASSOCIATED(block))
            CALL mpas_pool_get_subpool(block % structs, 'tend', tendPool)
            CALL mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
            CALL mpas_pool_get_subpool(block % structs, 'state', statePool)
            CALL mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
            CALL mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
            CALL mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
            CALL mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            CALL mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
            CALL mpas_pool_get_subpool(block % structs, 'shortwave', swForcingPool)
            CALL ocn_tend_tracer(tendPool, statePool, forcingPool, diagnosticsPool, &
                                 meshPool, swForcingPool, scratchPool, dt, &
                                 activeTracersOnly, 2)
            block => block % next
         END DO
         
         CALL mpas_timer_stop('se tracer tend')
         CALL mpas_threading_barrier()

         ! Update halo for tracer tendencies
         CALL mpas_timer_start("se halo tracers")
         CALL mpas_pool_get_subpool(domain % blocklist % structs, 'tend', tendPool)
         CALL mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)

         CALL mpas_pool_begin_iteration(tracersTendPool)
         DO WHILE (mpas_pool_get_next_member(tracersTendPool, groupItr))
            IF (groupItr % memberType == MPAS_POOL_FIELD) THEN
               ! Only compute tendencies for active tracers if activeTracersOnly flag is 
               ! true.
               IF (.NOT. activeTracersOnly .OR. &
                   TRIM(groupItr % memberName)=='activeTracersTend') THEN
                  CALL mpas_dmpar_field_halo_exch(domain, groupItr % memberName)
               END IF
            END IF
         END DO
         
         CALL mpas_timer_stop("se halo tracers")
         CALL mpas_threading_barrier()
         CALL mpas_timer_start('se loop fini')
         block => domain % blocklist
         
         DO WHILE (ASSOCIATED(block))
         
            CALL mpas_pool_get_DIMENSION(block % dimensions, 'nCells', nCellsPtr)
            CALL mpas_pool_get_DIMENSION(block % dimensions, 'nEdges', nEdgesPtr)
            CALL mpas_pool_get_DIMENSION(block % dimensions, 'nCellsArray', nCellsArray)
            CALL mpas_pool_get_DIMENSION(block % dimensions, 'nEdgesArray', nEdgesArray)
            CALL mpas_pool_get_DIMENSION(block % dimensions, 'nVertLevels', nVertLevels)
            CALL mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            CALL mpas_pool_get_subpool(block % structs, 'state', statePool)
            CALL mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
            CALL mpas_pool_get_subpool(block % structs, 'tend', tendPool)
            CALL mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
            CALL mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
            CALL mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
            CALL mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
            CALL mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
            CALL mpas_pool_get_array(meshPool, 'edgeMask', edgeMask)
            CALL mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
            CALL mpas_pool_get_array(tracersPool, 'activeTracers', tracersGroupCur, 1)
            CALL mpas_pool_get_array(tracersPool, 'activeTracers', tracersGroupNew, 2)
            CALL mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
            CALL mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)
            CALL mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
            CALL mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
            CALL mpas_pool_get_array(statePool, 'highFreqThickness', &
                                     highFreqThicknessCur, 1)
            CALL mpas_pool_get_array(statePool, 'highFreqThickness', &
                                     highFreqThicknessNew, 2)
            CALL mpas_pool_get_array(statePool, 'lowFreqDivergence', &
                                     lowFreqDivergenceCur, 1)
            CALL mpas_pool_get_array(statePool, 'lowFreqDivergence', &
                                     lowFreqDivergenceNew, 2)
            CALL mpas_pool_get_array(statePool, 'normalBarotropicVelocity', &
                                     normalBarotropicVelocityCur, 1)
            CALL mpas_pool_get_array(statePool, 'normalBarotropicVelocity', &
                                     normalBarotropicVelocityNew, 2)
            CALL mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', &
                                     normalBaroclinicVelocityCur, 1)
            CALL mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', &
                                     normalBaroclinicVelocityNew, 2)
            CALL mpas_pool_get_array(tendPool, 'layerThickness', layerThicknessTend)
            CALL mpas_pool_get_array(tendPool, 'normalVelocity', normalVelocityTend)
            CALL mpas_pool_get_array(tendPool, 'highFreqThickness', &
                                     highFreqThicknessTend)
            CALL mpas_pool_get_array(tendPool, 'lowFreqDivergence', &
                                     lowFreqDivergenceTend)
            CALL mpas_pool_get_array(tracersTendPool, 'activeTracersTend', &
                                     activeTracersTend)

            nCells = nCellsPtr
            nEdges = nEdgesPtr
            !
            ! If iterating, reset variables for next iteration
            !
            IF (split_explicit_step < config_n_ts_iter) THEN

               ! Get indices for dynamic tracers (Includes T&S).
               CALL mpas_pool_get_DIMENSION(tracersPool, 'activeGRP_start', startIndex)
               CALL mpas_pool_get_DIMENSION(tracersPool, 'activeGRP_end', endIndex)

               ! Only need T & S for earlier iterations,
               ! then all the tracers needed the last time through.

               !$OMP DO SCHEDULE(runtime) PRIVATE(k, temp_h, temp, i)
               DO iCell = 1, nCells
                  ! sshNew is a pointer, defined above.
                  DO k = 1, maxLevelCell(iCell)
                     ! This is h_{n+1}
                     temp_h = layerThicknessCur(k,iCell) &
                              + dt*layerThicknessTend(k,iCell)
                     ! This is h_{n+1/2}
                     layerThicknessNew(k,iCell) = 0.5*(layerThicknessCur(k,iCell) &
                                                       + temp_h)
                     DO i = startIndex, endIndex
                        ! This is Phi at n+1
                        temp = (tracersGroupCur(i,k,iCell)*layerThicknessCur(k,iCell) &
                               + dt*activeTracersTend(i,k,iCell))/temp_h
                        ! This is Phi at n+1/2
                        tracersGroupNew(i,k,iCell) &
                        = 0.5_RKIND*(tracersGroupCur(i,k,iCell) + temp)
                     END DO
                  END DO
               END DO ! iCell
               !$OMP END DO

               IF (config_use_freq_filtered_thickness) THEN
                  !$OMP DO SCHEDULE(runtime) PRIVATE(k, temp)
                  DO iCell = 1, nCells
                     DO k = 1, maxLevelCell(iCell)
                        ! h^{hf}_{n+1} was computed in Stage 1
                        ! This is h^{hf}_{n+1/2}
                        highFreqThicknessnew(k,iCell) &
                        = 0.5_RKIND*(highFreqThicknessCur(k,iCell) &
                          + highFreqThicknessNew(k,iCell))
                        ! This is D^{lf}_{n+1}
                        temp = lowFreqDivergenceCur(k,iCell) &
                               + dt*lowFreqDivergenceTend(k,iCell)
                        ! This is D^{lf}_{n+1/2}
                        lowFreqDivergenceNew(k,iCell) &
                        = 0.5_RKIND*(lowFreqDivergenceCur(k,iCell) + temp)
                     END DO
                  END DO
                  !$OMP END DO
               END IF

               !$OMP DO SCHEDULE(runtime) PRIVATE(k)
               DO iEdge = 1, nEdges
                  DO k = 1, nVertLevels
                     ! u = normalBarotropicVelocity + normalBaroclinicVelocity
                     ! here normalBaroclinicVelocity is at time n+1/2
                     ! This is u used in next iteration or step
                     normalVelocityNew(k,iEdge) &
                     = edgeMask(k,iEdge)*(normalBarotropicVelocityNew(iEdge) &
                                          + normalBaroclinicVelocityNew(k,iEdge))
                  END DO
               END DO ! iEdge
               !$OMP END DO

               ! Efficiency note: We really only need this to compute 
               ! layerThicknessEdge, density, pressure, and SSH in this diagnostics 
               ! solve.
               CALL ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool, &
                                         diagnosticsPool, scratchPool, tracersPool, 2)
            !
            !  If large iteration complete, compute all variables at time n+1
            !
            ELSE IF (split_explicit_step == config_n_ts_iter) THEN

               !$OMP DO SCHEDULE(runtime) PRIVATE(k)
               DO iCell = 1, nCells
                  DO k = 1, maxLevelCell(iCell)
                     ! This is h_{n+1}
                     layerThicknessNew(k,iCell) = layerThicknessCur(k,iCell) &
                                                  + dt*layerThicknessTend(k,iCell)
                  END DO
               END DO
               !$OMP END DO

               IF (config_compute_active_tracer_budgets) THEN
                  CALL mpas_pool_get_array(diagnosticsPool, &
                                           'activeTracerHorizontalAdvectionTendency', &
                                           activeTracerHorizontalAdvectionTendency)
                  CALL mpas_pool_get_array(diagnosticspool, &
                                           'activeTracerVerticalAdvectionTendency', &
                                           activeTracerVerticalAdvectionTendency)
                  CALL mpas_pool_get_array(diagnosticsPool, &
                                           'activeTracerSurfaceFluxTendency', &
                                           activeTracerSurfaceFluxTendency)
                  CALL mpas_pool_get_array(diagnosticsPool, &
                                           'temperatureShortWaveTendency', &
                                           temperatureShortWaveTendency)
                  CALL mpas_pool_get_array(diagnosticsPool, &
                                           'activeTracerNonLocalTendency', &
                                           activeTracerNonLocalTendency)
                  !$OMP DO SCHEDULE(runtime) PRIVATE(k)
                  DO iCell = 1, nCells
                     DO k= 1, maxLevelCell(iCell)
                        activeTracerHorizontalAdvectionTendency(:,k,iCell) &
                        = activeTracerHorizontalAdvectionTendency(:,k,iCell) &
                          /layerThicknessNew(k,iCell)
        
                        activeTracerVerticalAdvectionTendency(:,k,iCell) &
                        = activeTracerVerticalAdvectionTendency(:,k,iCell) &
                          /layerThicknessNew(k,iCell)

                        activeTracerSurfaceFluxTendency(:,k,iCell) &
                        = activeTracerSurfaceFluxTendency(:,k,iCell) &
                          /layerThicknessNew(k,iCell)

                        temperatureShortWaveTendency(k,iCell) &
                        = temperatureShortWaveTendency(k,iCell) &
                          /layerThicknessNew(k,iCell)

                        activeTracerNonLocalTendency(:,k,iCell) &
                        = activeTracerNonLocalTendency(:,k,iCell) &
                          /layerThicknessNew(k,iCell)
                     END DO
                  END DO 
                  !$OMP END DO
               END IF

               CALL mpas_pool_begin_iteration(tracersPool)
               DO WHILE (mpas_pool_get_next_member(tracersPool, groupItr))
                  IF (groupItr % memberType == MPAS_POOL_FIELD) THEN
                     configName = 'config_use_' // TRIM(groupItr % memberName)
                     CALL mpas_pool_get_config(domain % configs, configName, &
                                               config_use_tracerGroup)
                     IF (config_use_tracerGroup) THEN
                        CALL mpas_pool_get_array(tracersPool, groupItr % memberName, &
                                                 tracersGroupCur, 1)
                        CALL mpas_pool_get_array(tracersPool, groupItr % memberName, &
                                                 tracersGroupNew, 2)
                        modifiedGroupName = TRIM(groupItr % memberName) // 'Tend'
                        CALL mpas_pool_get_array(tracersTendPool, modifiedGroupName, &
                                                 tracersGroupTend)
                        !$OMP DO SCHEDULE(runtime) PRIVATE(k)
                        DO iCell = 1, nCells
                           DO k = 1, maxLevelCell(iCell)
                              tracersGroupNew(:,k,iCell) &
                              = (tracersGroupCur(:,k,iCell)*layerThicknessCur(k,iCell) &
                                 + dt*tracersGroupTend(:,k,iCell)) &
                                /layerThicknessNew(k,iCell)
                           END DO
                        END DO
                        !$OMP END DO
                        
                        ! Limit salinity in separate loop
                        IF (TRIM(groupItr % memberName) == 'activeTracers') THEN
                           !$OMP DO SCHEDULE(runtime) PRIVATE(k)
                           DO iCell = 1, nCells
                              DO k = 1, maxLevelCell(iCell)
                                 tracersGroupNew(indexSalinity,k,iCell) &
                                 = MAX(0.001_RKIND, &
                                       tracersGroupNew(indexSalinity,k,iCell))
                              END DO
                           END DO
                           !$OMP END DO
                        END IF

                     END IF
                  END IF
               END DO

               IF (config_use_freq_filtered_thickness) THEN
                  !$OMP DO SCHEDULE(runtime) PRIVATE(k)
                  DO iCell = 1, nCells
                     DO k = 1, maxLevelCell(iCell)
                        ! h^{hf}_{n+1} was computed in Stage 1
                        ! This is D^{lf}_{n+1}
                        lowFreqDivergenceNew(k,iCell) &
                        = lowFreqDivergenceCur(k,iCell) &
                          + dt*lowFreqDivergenceTend(k,iCell)
                     END DO
                  END DO
                  !$OMP END DO
               END IF

               ! Recompute final u to go on to next step.
               ! u_{n+1} = normalBarotropicVelocity_{n+1} 
               !           + normalBaroclinicVelocity_{n+1}
               ! Right now normalBaroclinicVelocityNew is at time n+1/2, so back compute 
               ! to get normalBaroclinicVelocity at time n+1 using 
               ! normalBaroclinicVelocity_{n+1/2} = 1/2*(normalBaroclinicVelocity_n 
               !                                         + u_Bcl_{n+1})
               ! So the following lines are
               ! u_{n+1} = normalBarotropicVelocity_{n+1} 
               !           + 2*normalBaroclinicVelocity_{n+1/2} 
               !           - normalBaroclinicVelocity_n
               ! Note that normalBaroclinicVelocity is recomputed at the beginning of 
               ! the next timestep due to Imp Vert mixing, so normalBaroclinicVelocity 
               ! does not have to be recomputed here.

               !$OMP DO SCHEDULE(runtime) PRIVATE(k)
               DO iEdge = 1, nEdges
                  DO k = 1, maxLevelEdgeTop(iEdge)
                     normalVelocityNew(k,iEdge) &
                     = normalBarotropicVelocityNew(iEdge) &
                       + 2*normalBaroclinicVelocityNew(k,iEdge) &
                       - normalBaroclinicVelocityCur(k,iEdge)
                  END DO
               END DO ! iEdges
               !$OMP END DO

            END IF ! split_explicit_step

            block => block % next
            
         END DO

         CALL mpas_timer_stop('se loop fini')
         CALL mpas_timer_stop('se loop')

      END DO ! split_explicit_step = 1, config_n_ts_iter
      !
      ! END large iteration loop
      !
      CALL mpas_timer_start("se implicit vert mix")

      block => domain % blocklist
      
      DO WHILE(ASSOCIATED(block))
      
        CALL mpas_pool_get_subpool(block % structs, 'state', statePool)
        CALL mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
        CALL mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
        CALL mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
        CALL mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
        CALL mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)

        ! Call ocean diagnostic solve in preparation for vertical mixing.  Note
        ! it is called again after vertical mixing, because u and tracers change.
        ! For Richardson vertical mixing, only density, layerThicknessEdge, and 
        ! kineticEnergyCell need to be computed.  For kpp, more variables may be needed.  
        ! Either way, this could be made more efficient by only computing what is needed 
        ! for the implicit vmix routine that follows.
        CALL ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool, &
                                  diagnosticsPool, scratchPool, tracersPool, 2)

        ! Compute normalGMBolusVelocity; it will be added to the baroclinic modes in 
        ! Stage 2 above.
        IF (config_use_standardGM) THEN
           CALL ocn_gm_compute_Bolus_velocity(diagnosticsPool, meshPool, scratchPool)
        END IF
        CALL ocn_vmix_implicit(dt, meshPool, diagnosticsPool, statePool, forcingPool, &
                               scratchPool, err, 2)

        block => block % next
        
      END DO

      CALL mpas_threading_barrier()

      ! Update halo on u and tracers, which were just updated for implicit vertical 
      ! mixing. If not done, this leads to lack of volume conservation. It is required
      ! because halo updates in stage 3 are only conducted on tendencies, not on the 
      ! velocity and tracer fields.  So this update is required to communicate the 
      ! change due to implicit vertical mixing across the boundary.
      CALL mpas_timer_start('se vmix halos')
      CALL mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)
      CALL mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

      CALL mpas_timer_start('se vmix halos normalVelFld')
      CALL mpas_dmpar_field_halo_exch(domain, 'normalVelocity', timeLevel=2)
      CALL mpas_timer_stop('se vmix halos normalVelFld')

      CALL mpas_pool_begin_iteration(tracersPool)
      DO WHILE (mpas_pool_get_next_member(tracersPool, groupItr))
         IF (groupItr % memberType == MPAS_POOL_FIELD) THEN
            CALL mpas_dmpar_field_halo_exch(domain, groupItr % memberName, timeLevel=2)
         END IF
      END DO
      CALL mpas_timer_stop('se vmix halos')

      CALL mpas_timer_stop("se implicit vert mix")

      CALL mpas_threading_barrier()

      CALL mpas_timer_start('se fini')
      block => domain % blocklist
      
      DO WHILE (ASSOCIATED(block))
      
         CALL mpas_pool_get_subpool(block % structs, 'state', statePool)
         CALL mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         CALL mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
         CALL mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         CALL mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         CALL mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)

         CALL mpas_pool_get_DIMENSION(block % dimensions, 'nCells', nCellsPtr)
         CALL mpas_pool_get_DIMENSION(block % dimensions, 'nEdges', nEdgesPtr)
         CALL mpas_pool_get_DIMENSION(block % dimensions, 'nCellsArray', nCellsArray)
         CALL mpas_pool_get_DIMENSION(block % dimensions, 'nEdgesArray', nEdgesArray)

         CALL mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
         CALL mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
         CALL mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
         CALL mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)
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

         CALL mpas_pool_get_DIMENSION(diagnosticsPool, 'index_surfaceVelocityZonal', &                     
                                      indexSurfaceVelocityZonal)
         CALL mpas_pool_get_DIMENSION(diagnosticsPool, &
                                      'index_surfaceVelocityMeridional', &
                                      indexSurfaceVelocityMeridional)
         CALL mpas_pool_get_DIMENSION(diagnosticsPool, 'index_SSHGradientZonal', &
                                      indexSSHGradientZonal)
         CALL mpas_pool_get_DIMENSION(diagnosticsPool, 'index_SSHGradientMeridional', &
                                      indexSSHGradientMeridional)

         nCells = nCellsPtr
         nEdges = nEdgesPtr

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

         ! Update the effective desnity in land ice if we're coupling to land ice
         CALL ocn_effective_density_in_land_ice_update(meshPool, forcingPool, &
                                                       statePool, scratchPool, err)

         ! Compute normalGMBolusVelocity; it will be added to normalVelocity in Stage 2 
         ! of the next cycle.
         IF (config_use_standardGM) THEN
            CALL ocn_gm_compute_Bolus_velocity(diagnosticsPool, meshPool, scratchPool)
         END IF

         CALL mpas_timer_start('se final mpas reconstruct', .FALSE.)

         CALL mpas_reconstruct(meshPool, normalVelocityNew, velocityX, velocityY, &
                               velocityZ, velocityZonal, velocityMeridional, &
                               includeHalos = .TRUE.)

         CALL mpas_reconstruct(meshPool, gradSSH, gradSSHX, gradSSHY, gradSSHZ, &
                               gradSSHZonal, gradSSHMeridional, includeHalos = .TRUE.)

         CALL mpas_timer_stop('se final mpas reconstruct')

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
         CALL mpas_threading_barrier()

         IF (config_use_standardGM) THEN
            CALL ocn_reconstruct_gm_vectors(diagnosticsPool, meshPool)
         END IF

         block => block % next
         
      END DO

      IF (TRIM(config_land_ice_flux_mode) == 'coupled') THEN
         CALL mpas_timer_start("se effective density halo")
         CALL mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)
         CALL mpas_pool_get_field(statePool, 'effectiveDensityInLandIce', &
                                  effectiveDensityField, 2)
         CALL mpas_dmpar_exch_halo_field(effectiveDensityField)
         CALL mpas_timer_stop("se effective density halo")
      END IF

      CALL mpas_timer_stop('se fini')
      CALL mpas_timer_stop("se timestep")

      DEALLOCATE(n_bcl_iter)

   END SUBROUTINE ocn_time_integrator_split
   !-----------------------------------------------------------------------------------!
   !                                                                                   !
   !-----------------------------------------------------------------------------------!
   !
   !  Routine ocn_time_integration_split_init
   !
   !> \brief   Initialize split-explicit time stepping within MPAS-Ocean core
   !> \author  Mark Petersen
   !> \date    September 2011
   !> \details
   !> This routine initializes variables required for the split-explicit time stepper.
   !
   SUBROUTINE ocn_time_integration_split_init(domain)
   
      IMPLICIT NONE
   
      ! Initialize splitting variables

      TYPE(domain_type), INTENT(INOUT) :: domain
      
      ! Local variables
      CHARACTER(LEN=StrKIND), POINTER :: config_time_integrator, config_btr_dt, &
                                         config_dt
      INTEGER :: i, iCell, iEdge, iVertex, k, iTracer, cell, cell1, cell2, iErr
      INTEGER(KIND=I8KIND) :: nBtrSubcyclesI8
      INTEGER, DIMENSION(:), POINTER :: maxLevelEdgeTop
      INTEGER, DIMENSION(:,:), POINTER :: cellsOnEdge      
      INTEGER, POINTER :: nVertLevels, nCells, nEdges
      LOGICAL, POINTER :: config_filter_btr_mode, config_do_restart
      REAL(KIND=RKIND) :: normalThicknessFluxSum, layerThicknessSum, layerThicknessEdge1
      REAL(KIND=RKIND), DIMENSION(:), POINTER :: refBottomDepth, &
                                                 normalBarotropicVelocity
      REAL(KIND=RKIND), DIMENSION(:,:), POINTER &
      :: layerThickness, normalBaroclinicVelocity, normalVelocity      

      TYPE(mpas_time_type) :: nowTime
      TYPE(mpas_timeInterval_type) :: fullTimeStep, barotropicTimeStep, remainder, &
                                      zeroInterval
      TYPE(block_type), POINTER :: block
      TYPE(mpas_pool_type), POINTER :: statePool, meshPool, tracersPool
      
      CALL mpas_pool_get_config(domain % configs, 'config_do_restart', &
                                config_do_restart)

      ! Determine the number of barotropic subcycles based on the ratio of time steps
      CALL mpas_pool_get_config(domain % configs, 'config_time_integrator', &
                                config_time_integrator)
      CALL mpas_pool_get_config(domain % configs, 'config_btr_dt', config_btr_dt)
      CALL mpas_pool_get_config(domain % configs, 'config_dt', config_dt)

      nowTime = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr)
      CALL mpas_set_timeInterval(zeroInterval, S=0)

      CALL mpas_set_timeInterval(fullTimeStep , timeString=config_dt)
      CALL mpas_set_timeInterval(barotropicTimeStep, timeString=config_btr_dt)

      ! transfer to I8 for division step
      nBtrSubcyclesI8 = nBtrSubcycles
      CALL mpas_interval_division(nowTime, fullTimeStep, barotropicTimeStep, &
                                  nBtrSubcyclesI8, remainder)
      nBtrSubcycles = nBtrSubcyclesI8

      IF (remainder > zeroInterval) THEN
         nBtrSubcycles = nBtrSubcycles + 1
      END IF

      IF (TRIM(config_time_integrator) == 'split_explicit') THEN
         CALL mpas_log_write('********************************************************&
                              &***********************')
         CALL mpas_log_write('The split explicit time integration is configured to &
                              &use: $i barotropic subcycles', &
                             intArgs=(/nBtrSubcycles/))
         CALL mpas_log_write('********************************************************&
                              &***********************')
      END IF

      IF (.NOT. config_do_restart) THEN
      
         ! Initialize z-level mesh variables from h, read in from input file.
         block => domain % blocklist
         
         DO WHILE (ASSOCIATED(block))
         
            CALL mpas_pool_get_config(block % configs, 'config_time_integrator', &
                                      config_time_integrator)
            CALL mpas_pool_get_config(block % configs, 'config_filter_btr_mode', &
                                      config_filter_btr_mode)
            CALL mpas_pool_get_subpool(block % structs, 'state', statePool)
            CALL mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
            CALL mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            CALL mpas_pool_get_DIMENSION(block % dimensions, 'nVertLevels', nVertLevels)
            CALL mpas_pool_get_DIMENSION(block % dimensions, 'nCells', nCells)
            CALL mpas_pool_get_DIMENSION(block % dimensions, 'nEdges', nEdges)
            CALL mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)
            CALL mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, 1)
            CALL mpas_pool_get_array(statePool, 'normalBarotropicVelocity', &
                                     normalBarotropicVelocity, 1)
            CALL mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', &
                                     normalBaroclinicVelocity, 1)
            CALL mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
            CALL mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
            CALL mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)

            ! Compute barotropic velocity at first timestep
            ! This is only done upon start-up.
            IF (TRIM(config_time_integrator) == 'unsplit_explicit') THEN
               CALL mpas_pool_get_array(statePool, 'normalBarotropicVelocity', &
                                        normalBarotropicVelocity)

               DO iEdge = 1, nEdges
                  normalBarotropicVelocity(iEdge) = 0.0_RKIND
                  normalBaroclinicVelocity(:, iEdge) = normalVelocity(:, iEdge)
               END DO

            ELSE IF (TRIM(config_time_integrator) == 'split_explicit') THEN

               CALL mpas_log_write('**************************************************&
                                    &*****************************')
               CALL mpas_log_write('The split explicit time integration is configured &
                                    &to use: $i barotropic subcycles', &
                                   intArgs=(/nBtrSubcycles/))
               CALL mpas_log_write('**************************************************&
                                    &*****************************')

               IF (config_filter_btr_mode) THEN
                  DO iCell = 1, nCells
                     layerThickness(1,iCell) = refBottomDepth(1)
                  END DO
               END IF

               DO iEdge = 1, nEdges
                  cell1 = cellsOnEdge(1,iEdge)
                  cell2 = cellsOnEdge(2,iEdge)

                  ! normalBarotropicVelocity = sum(h*u)/sum(h) on each edge
                  ! ocn_diagnostic_solve has not yet been called, so compute hEdge
                  ! just for this edge.

                  ! thicknessSum is initialized outside the loop because on land 
                  ! boundaries
                  ! maxLevelEdgeTop = 0, but I want to initialize thicknessSum with a
                  ! nonzero value to avoid a NaN.
                  layerThicknessEdge1 = 0.5_RKIND*(layerThickness(1,cell1) &
                                                   + layerThickness(1,cell2))
                  normalThicknessFluxSum = layerThicknessEdge1*normalVelocity(1,iEdge)
                  layerThicknessSum = layerThicknessEdge1

                  DO k=2, maxLevelEdgeTop(iEdge)
                     ! ocn_diagnostic_solve has not yet been called, so compute hEdge
                     ! just for this edge.
                     layerThicknessEdge1 = 0.5_RKIND*(layerThickness(k,cell1) &
                                                      + layerThickness(k,cell2))

                     normalThicknessFluxSum = normalThicknessFluxSum &
                        + layerThicknessEdge1*normalVelocity(k,iEdge)
                     layerThicknessSum = layerThicknessSum + layerThicknessEdge1

                  END DO
                  normalBarotropicVelocity(iEdge) &
                  = normalThicknessFluxSum/layerThicknessSum

                  ! normalBaroclinicVelocity(k,iEdge) 
                  ! = normalVelocity(k,iEdge) - normalBarotropicVelocity(iEdge)
                  DO k = 1, maxLevelEdgeTop(iEdge)
                     normalBaroclinicVelocity(k,iEdge) &
                     = normalVelocity(k,iEdge) - normalBarotropicVelocity(iEdge)
                  END DO

                  ! normalBaroclinicVelocity=0, normalVelocity=0 on land cells
                  DO k = maxLevelEdgeTop(iEdge)+1, nVertLevels
                     normalBaroclinicVelocity(k,iEdge) = 0.0_RKIND
                     normalVelocity(k,iEdge) = 0.0_RKIND
                  END DO
               END DO

               IF (config_filter_btr_mode) THEN
                  ! filter normalBarotropicVelocity out of initial condition

                   normalVelocity(:,:) = normalBaroclinicVelocity(:,:)
                   normalBarotropicVelocity(:) = 0.0_RKIND

               END IF

            END IF

         block => block % next
         
         END DO
         
      END IF

   END SUBROUTINE ocn_time_integration_split_init

END MODULE ocn_time_integration_split
