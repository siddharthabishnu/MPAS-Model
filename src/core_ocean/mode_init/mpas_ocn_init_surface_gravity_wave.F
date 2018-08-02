!
! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!  ocn_init_surface_gravity_wave
!
!> \brief   MPAS ocean initialize case -- surface gravity wave
!> \author  Siddhartha Bishnu
!> \date    07/23/2018
!> \details
!> This module contains the routines for initializing the surface gravity wave test case.
!
MODULE ocn_init_surface_gravity_wave

   USE mpas_kind_types
   USE mpas_io_units
   USE mpas_derived_types
   USE mpas_pool_routines
   USE mpas_constants
   USE mpas_dmpar

   USE ocn_constants
   USE ocn_init_vertical_grids
   USE ocn_init_cell_markers

   IMPLICIT NONE
   PRIVATE
   SAVE

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

   PUBLIC :: ocn_init_setup_surface_gravity_wave, ocn_init_validate_surface_gravity_wave

   !------------------------------------------------------------------------------------!
   !
   ! Private module variables
   !
   !------------------------------------------------------------------------------------!
    
   CONTAINS 
   !
   !  routine ocn_init_setup_surface_gravity_wave
   !
   !> \brief   Setup for surface_gravity_wave test case
   !> \author  Siddhartha Bishnu
   !> \date    07/23/2018
   !> \details
   !> This routine sets up the initial conditions for the surface gravity wave test 
   !> case. It should also ensure the mesh that was input is valid for the 
   !> configuration.
   !
   SUBROUTINE ocn_init_setup_surface_gravity_wave(domain, iErr)
   
      IMPLICIT NONE

      TYPE(domain_type), INTENT(INOUT) :: domain
      INTEGER, INTENT(OUT) :: iErr
      
      ! Local variables
      REAL(KIND=RKIND) :: yMin, yMax, dcEdgeMin, dcEdgeMinGlobal, yMinGlobal, &
                          yMaxGlobal, yMidGlobal, yOffset, xPerturbationMin, &
                          temperature

      TYPE(block_type), POINTER :: block_ptr
      TYPE(mpas_pool_type), POINTER :: meshPool, statePool, tracersPool, verticalMeshPool

      INTEGER :: iCell, k, idx

      ! Define config variable pointers
      CHARACTER(LEN=StrKIND), POINTER :: config_init_configuration, config_vertical_grid
      LOGICAL, POINTER :: config_surface_gravity_wave_use_distances
      REAL(KIND=RKIND), POINTER :: config_surface_gravity_wave_bottom_depth, &
                                   config_surface_gravity_wave_surface_temperature, &
                                   config_surface_gravity_wave_bottom_temperature, &
                                   config_surface_gravity_wave_salinity, &
                                   config_surface_gravity_wave_coriolis_parameter, &
                                   config_surface_elevation_to_mean_depth_ratio, &
                                   config_Gaussian_surface_elevation_decay_scale

      ! Define dimension pointers
      INTEGER, POINTER :: nCellsSolve, nEdgesSolve, nVertLevels, nVertLevelsP1, &
                          index_temperature, index_salinity, index_tracer1

      ! Define variable pointers
      INTEGER, DIMENSION(:), POINTER :: maxLevelCell
      REAL(KIND=RKIND), DIMENSION(:), POINTER :: yCell, refBottomDepth, refZMid, &
                                                 vertCoordMovementWeights, bottomDepth, &
                                                 fCell, fEdge, fVertex, dcEdge
      REAL(KIND=RKIND), DIMENSION(:,:), POINTER :: layerThickness, restingThickness
      REAL(KIND=RKIND), DIMENSION(:,:,:), POINTER :: activeTracers, debugTracers

      ! Define local interfaceLocations variable
      REAL(KIND=RKIND), DIMENSION(:), POINTER :: interfaceLocations

      iErr = 0

      CALL mpas_pool_get_config(ocnConfigs, 'config_init_configuration', &
                                config_init_configuration)

      IF (config_init_configuration /= TRIM('surface_gravity_wave')) RETURN

      CALL mpas_pool_get_config(ocnConfigs, 'config_vertical_grid', config_vertical_grid)

      CALL mpas_pool_get_config(ocnConfigs, &
                                'config_surface_gravity_wave_use_distances', &
                                config_surface_gravity_wave_use_distances)
      CALL mpas_pool_get_config(ocnConfigs, 'config_surface_gravity_wave_bottom_depth', &
                                config_surface_gravity_wave_bottom_depth)
      CALL mpas_pool_get_config(ocnConfigs, &
                                'config_surface_gravity_wave_surface_temperature', &
                                config_surface_gravity_wave_surface_temperature)
      CALL mpas_pool_get_config(ocnConfigs, &
                                'config_surface_gravity_wave_bottom_temperature', &
                                config_surface_gravity_wave_bottom_temperature)
      CALL mpas_pool_get_config(ocnConfigs, 'config_surface_gravity_wave_salinity', &
                                config_surface_gravity_wave_salinity)
      CALL mpas_pool_get_config(ocnConfigs, &
                                'config_surface_gravity_wave_coriolis_parameter', &
                                config_surface_gravity_wave_coriolis_parameter)
      CALL mpas_pool_get_config(ocnConfigs, &
                                'config_surface_elevation_to_mean_depth_ratio', &
                                config_surface_elevation_to_mean_depth_ratio)
      CALL mpas_pool_get_config(ocnConfigs, &
                                'config_Gaussian_surface_elevation_decay_scale', &
                                config_Gaussian_surface_elevation_decay_scale)
         
      PRINT *, 'Stage 1: config_surface_elevation_to_mean_depth_ratio', &
               config_surface_elevation_to_mean_depth_ratio
      PRINT *, 'Stage 1: config_Gaussian_surface_elevation_decay_scale', &
               config_Gaussian_surface_elevation_decay_scale
                                
      ! Determine vertical grid for configuration
      CALL mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
      CALL mpas_pool_get_dimension(meshPool, 'nVertLevelsP1', nVertLevelsP1)

      ALLOCATE (interfaceLocations(nVertLevelsP1))
      CALL ocn_generate_vertical_grid(config_vertical_grid, interfaceLocations)

      ! Initalize min/max values to large positive and negative values
      yMin = 1.0E10_RKIND
      yMax = -1.0E10_RKIND
      dcEdgeMin = 1.0E10_RKIND

      ! Determine local min and max values.
      block_ptr => domain % blocklist
      
      DO WHILE (ASSOCIATED(block_ptr))
         CALL mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
         CALL mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
         CALL mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
         CALL mpas_pool_get_array(meshPool, 'yCell', yCell)
         CALL mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
         yMin = MIN(yMin, MINVAL(yCell(1:nCellsSolve)))
         yMax = MAX(yMax, MAXVAL(yCell(1:nCellsSolve)))
         dcEdgeMin = MIN(dcEdgeMin, MINVAL(dcEdge(1:nEdgesSolve)))
         block_ptr => block_ptr % next
      END DO

      ! Determine global min and max values.
      CALL mpas_dmpar_min_real(domain % dminfo, yMin, yMinGlobal)
      CALL mpas_dmpar_max_real(domain % dminfo, yMax, yMaxGlobal)
      CALL mpas_dmpar_min_real(domain % dminfo, dcEdgeMin, dcEdgeMinGlobal)

      yMidGlobal = (yMinGlobal + yMaxGlobal)*0.5_RKIND

      block_ptr => domain % blocklist
      
      DO WHILE (ASSOCIATED(block_ptr))
      
         CALL mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
         CALL mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
         CALL mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', &
                                    verticalMeshPool)
         CALL mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

         CALL mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
         CALL mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

         CALL mpas_pool_get_dimension(tracersPool, 'index_temperature', &
                                      index_temperature)
         CALL mpas_pool_get_dimension(tracersPool, 'index_salinity', index_salinity)
         CALL mpas_pool_get_dimension(tracersPool, 'index_tracer1', index_tracer1)

         CALL mpas_pool_get_array(meshPool, 'yCell', yCell)
         CALL mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
         CALL mpas_pool_get_array(meshPool, 'vertCoordMovementWeights', &
                                  vertCoordMovementWeights)
         CALL mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
         CALL mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
         CALL mpas_pool_get_array(meshPool, 'fCell', fCell)
         CALL mpas_pool_get_array(meshPool, 'fEdge', fEdge)
         CALL mpas_pool_get_array(meshPool, 'fVertex', fVertex)

         CALL mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)
         CALL mpas_pool_get_array(tracersPool, 'debugTracers', debugTracers, 1)
         CALL mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)

         CALL mpas_pool_get_array(verticalMeshPool, 'refZMid', refZMid)
         CALL mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness)

         CALL ocn_mark_north_boundary(meshPool, yMaxGlobal, dcEdgeMinGlobal, iErr)
         CALL ocn_mark_south_boundary(meshPool, yMinGlobal, dcEdgeMinGlobal, iErr)

         PRINT *, 'Stage 2: config_surface_elevation_to_mean_depth_ratio', &
                  config_surface_elevation_to_mean_depth_ratio
         PRINT *, 'Stage 2: config_Gaussian_surface_elevation_decay_scale', &
                  config_Gaussian_surface_elevation_decay_scale
         
         ! Set refBottomDepth and refZMid
         DO k = 1, nVertLevels
            refBottomDepth(k) = config_surface_gravity_wave_bottom_depth &
                                *interfaceLocations(k+1)
            refZMid(k) = -0.5_RKIND*(interfaceLocations(k+1) + interfaceLocations(k)) &
                          *config_surface_gravity_wave_bottom_depth            
         END DO
         
         ! Set vertCoordMovementWeights
         vertCoordMovementWeights(:) = 1.0_RKIND        

         DO iCell = 1, nCellsSolve
         
            ! Set layerThickness and restingThickness
            DO k = 1, nVertLevels
               restingThickness(k, iCell) &
               = config_surface_gravity_wave_bottom_depth*(interfaceLocations(k+1) &
                                                           - interfaceLocations(k))          
!              layerThickness(k, iCell) = restingThickness(k, iCell)  
               layerThickness(k, iCell) &
               = restingThickness(k, iCell) &
                 *(1._RKIND &
                   + config_surface_elevation_to_mean_depth_ratio &
                     *EXP(-((yCell(iCell) - yMidGlobal) &
                            /config_Gaussian_surface_elevation_decay_scale)**2._RKIND))
!              layerThickness(k, iCell) &
!              = restingThickness(k, iCell) &
!                *(1._RKIND &
!                  + 0.01_RKIND &
!                    *EXP(-((yCell(iCell) - yMidGlobal) &
!                           /(5._RKIND*10._RKIND**4._RKIND))**2._RKIND))
!              layerThickness(k, iCell) &
!              = restingThickness(k, iCell) &
!                *(1._RKIND &
!                  + config_surface_elevation_to_mean_depth_ratio &
!                    *EXP(-((yCell(iCell) - yMidGlobal) &
!                           /(0.1_RKIND*(yMaxGlobal - yMinGlobal)))**2._RKIND))
            END DO

            ! Set debug tracer
            IF (ASSOCIATED(debugTracers)) THEN
               idx = index_tracer1
               DO k = 1, nVertLevels
                  debugTracers(idx, k, iCell) = 1.0_RKIND
               END DO
            END IF

            ! Set stratification based on northern half of domain temperature
            IF (ASSOCIATED(activeTracers)) THEN
               
               ! Set temperature
               idx = index_temperature
               DO k = nVertLevels, 1, -1
                  temperature = config_surface_gravity_wave_bottom_temperature &
                                + (config_surface_gravity_wave_surface_temperature &
                                   - config_surface_gravity_wave_bottom_temperature) &
                                *((refZMid(k) + refBottomDepth(nVertLevels)) &
                                /refBottomDepth(nVertLevels))
                  activeTracers(idx, k, iCell) = temperature
               END DO               

               ! Set salinity
               idx = index_salinity
               activeTracers(idx, :, iCell) = config_surface_gravity_wave_salinity
               
            END IF
            
            ! Set bottomDepth
            bottomDepth(iCell) = config_surface_gravity_wave_bottom_depth

            ! Set maxLevelCell
            maxLevelCell(iCell) = nVertLevels
            
         END DO
   
         ! Set Coriolis parameters
         fCell(:) = config_surface_gravity_wave_coriolis_parameter
         fEdge(:) = config_surface_gravity_wave_coriolis_parameter
         fVertex(:) = config_surface_gravity_wave_coriolis_parameter
   
         block_ptr => block_ptr % next
         
      END DO

      DEALLOCATE (interfaceLocations)

   END SUBROUTINE ocn_init_setup_surface_gravity_wave
   !------------------------------------------------------------------------------------!
   !                                                                                    !
   !------------------------------------------------------------------------------------!
   !
   !  routine ocn_init_validate_surface_gravity_wave
   !
   !> \brief   Validation for surface gravity wave test case
   !> \author  Siddhartha Bishnu
   !> \date    07/23/2018
   !> \details
   !> This routine validates the configuration options for the surface gravity wave test          
   !  case.
   !
   SUBROUTINE ocn_init_validate_surface_gravity_wave(configPool, packagePool, &
                                                     iocontext, iErr)

      IMPLICIT NONE

      TYPE(mpas_pool_type), INTENT(INOUT) :: configPool, packagePool
      TYPE(mpas_io_context_type), INTENT(INOUT) :: iocontext
      INTEGER, INTENT(OUT) :: iErr
      
      ! Local variables
      CHARACTER(LEN=StrKIND), POINTER :: config_init_configuration
      INTEGER, POINTER :: config_vert_levels, config_surface_gravity_wave_vert_levels

      iErr = 0

      CALL mpas_pool_get_config(configPool, 'config_init_configuration', &
                                config_init_configuration)

      IF (config_init_configuration /= TRIM('surface_gravity_wave')) RETURN

      CALL mpas_pool_get_config(configPool, 'config_vert_levels', config_vert_levels)
      CALL mpas_pool_get_config(configPool, 'config_surface_gravity_wave_vert_levels', &
                                config_surface_gravity_wave_vert_levels)

      IF (config_vert_levels <= 0 .AND. config_surface_gravity_wave_vert_levels > 0) THEN
         config_vert_levels = config_surface_gravity_wave_vert_levels
      ELSE IF (config_vert_levels <= 0) THEN
         CALL mpas_log_write('Validation failed for surface gravity wave. Not given a &
                              &usable value for vertical levels.', MPAS_LOG_CRIT)
         iErr = 1
      END IF

   END SUBROUTINE ocn_init_validate_surface_gravity_wave

END MODULE ocn_init_surface_gravity_wave
