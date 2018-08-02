!
! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!  ocn_init_baroclinic_channel
!
!> \brief   MPAS ocean initialize case -- Baroclinic Channel
!> \author  Doug Jacobsen
!> \date    02/18/2014
!> \details
!> This module contains the routines for initializing the baroclinic channel test case.
!
MODULE ocn_init_baroclinic_channel

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

   PUBLIC :: ocn_init_setup_baroclinic_channel, ocn_init_validate_baroclinic_channel

   !------------------------------------------------------------------------------------!
   !
   ! Private module variables
   !
   !------------------------------------------------------------------------------------!
    
   CONTAINS 
   !
   !  routine ocn_init_setup_baroclinic_channel
   !
   !> \brief   Setup for baroclinic channel test case
   !> \author  Doug Jacobsen
   !> \date    02/19/2014
   !> \details
   !>  This routine sets up the initial conditions for the baroclinic channel test case.
   !>  It should also ensure the mesh that was input is valid for the configuration.
   !
   SUBROUTINE ocn_init_setup_baroclinic_channel(domain, iErr)
   
      IMPLICIT NONE

      TYPE(domain_type), INTENT(INOUT) :: domain
      INTEGER, INTENT(OUT) :: iErr
      
      ! Local variables
      REAL(KIND=RKIND) :: yMin, yMax, xMin, xMax, dcEdgeMin, dcEdgeMinGlobal, &
                          yMinGlobal, yMaxGlobal, yMidGlobal, xMinGlobal, xMaxGlobal, &
                          temperature, yOffset, xPerturbationMin, xPerturbationMax, &
                          perturbationWidth

      TYPE(block_type), POINTER :: block_ptr
      TYPE(mpas_pool_type), POINTER :: meshPool, statePool, tracersPool, verticalMeshPool

      INTEGER :: iCell, k, idx

      ! Define config variable pointers
      CHARACTER(LEN=StrKIND), POINTER :: config_init_configuration, config_vertical_grid
      LOGICAL, POINTER :: config_baroclinic_channel_use_distances
      REAL(KIND=RKIND), POINTER :: config_baroclinic_channel_gradient_width_dist, &
                                   config_baroclinic_channel_gradient_width_frac, &
                                   config_baroclinic_channel_bottom_depth, &
                                   config_baroclinic_channel_surface_temperature, &
                                   config_baroclinic_channel_bottom_temperature, &
                                   config_baroclinic_channel_temperature_difference, &
                                   config_baroclinic_channel_salinity, &
                                   config_baroclinic_channel_coriolis_parameter

      ! Define dimension pointers
      INTEGER, POINTER :: nCellsSolve, nEdgesSolve, nVertLevels, nVertLevelsP1, &
                          index_temperature, index_salinity, index_tracer1

      ! Define variable pointers
      INTEGER, DIMENSION(:), POINTER :: maxLevelCell
      REAL(KIND=RKIND), DIMENSION(:), POINTER :: xCell, yCell, refBottomDepth, refZMid, &
                                                 vertCoordMovementWeights, bottomDepth, &
                                                 fCell, fEdge, fVertex, dcEdge
      REAL(KIND=RKIND), DIMENSION(:,:), POINTER :: layerThickness, restingThickness
      REAL(KIND=RKIND), DIMENSION(:,:,:), POINTER :: activeTracers, debugTracers

      ! Define local interfaceLocations variable
      REAL(KIND=RKIND), DIMENSION(:), POINTER :: interfaceLocations

      LOGICAL, POINTER :: on_a_sphere

      iErr = 0

      CALL mpas_pool_get_config(ocnConfigs, 'config_init_configuration', &
                                config_init_configuration)

      IF (config_init_configuration /= TRIM('baroclinic_channel')) RETURN

      CALL mpas_pool_get_config(ocnConfigs, 'config_vertical_grid', config_vertical_grid)

      CALL mpas_pool_get_config(ocnConfigs, 'config_baroclinic_channel_use_distances', &
                                config_baroclinic_channel_use_distances)
      CALL mpas_pool_get_config(ocnConfigs, &
                                'config_baroclinic_channel_gradient_width_dist', &
                                config_baroclinic_channel_gradient_width_dist)
      CALL mpas_pool_get_config(ocnConfigs, &
                                'config_baroclinic_channel_gradient_width_frac', &
                                config_baroclinic_channel_gradient_width_frac)
      CALL mpas_pool_get_config(ocnConfigs, 'config_baroclinic_channel_bottom_depth', &
                                config_baroclinic_channel_bottom_depth)
      CALL mpas_pool_get_config(ocnConfigs, &
                                'config_baroclinic_channel_surface_temperature', &
                                config_baroclinic_channel_surface_temperature)
      CALL mpas_pool_get_config(ocnConfigs, &
                                'config_baroclinic_channel_bottom_temperature', &
                                config_baroclinic_channel_bottom_temperature)
      CALL mpas_pool_get_config(ocnConfigs, &
                                'config_baroclinic_channel_temperature_difference', &
                                config_baroclinic_channel_temperature_difference)
      CALL mpas_pool_get_config(ocnConfigs, 'config_baroclinic_channel_salinity', &
                                config_baroclinic_channel_salinity)
      CALL mpas_pool_get_config(ocnConfigs, &
                                'config_baroclinic_channel_coriolis_parameter', &
                                config_baroclinic_channel_coriolis_parameter)

      ! Determine vertical grid for configuration
      CALL mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
      CALL mpas_pool_get_dimension(meshPool, 'nVertLevelsP1', nVertLevelsP1)
      CALL mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)

      IF (on_a_sphere) &
      CALL mpas_log_write('The baroclinic channel configuration ' &
                          // 'can only be applied to a planar mesh. Exiting...', &
                          MPAS_LOG_CRIT)

      ALLOCATE (interfaceLocations(nVertLevelsP1))
      CALL ocn_generate_vertical_grid(config_vertical_grid, interfaceLocations)

      ! Initalize min/max values to large positive and negative values
      yMin = 1.0E10_RKIND
      yMax = -1.0E10_RKIND
      xMin = 1.0E10_RKIND
      xMax = -1.0E10_RKIND
      dcEdgeMin = 1.0E10_RKIND

      ! Determine local min and max values.
      block_ptr => domain % blocklist
      
      DO WHILE (ASSOCIATED(block_ptr))
         CALL mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
         CALL mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
         CALL mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
         CALL mpas_pool_get_array(meshPool, 'xCell', xCell)
         CALL mpas_pool_get_array(meshPool, 'yCell', yCell)
         CALL mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
         yMin = MIN(yMin, MINVAL(yCell(1:nCellsSolve)))
         yMax = MAX(yMax, MAXVAL(yCell(1:nCellsSolve)))
         xMin = MIN(xMin, MINVAL(xCell(1:nCellsSolve)))
         xMax = MAX(xMax, MAXVAL(xCell(1:nCellsSolve)))
         dcEdgeMin = MIN(dcEdgeMin, MINVAL(dcEdge(1:nEdgesSolve)))
         block_ptr => block_ptr % next
      END DO

      ! Determine global min and max values.
      CALL mpas_dmpar_min_real(domain % dminfo, yMin, yMinGlobal)
      CALL mpas_dmpar_max_real(domain % dminfo, yMax, yMaxGlobal)
      CALL mpas_dmpar_min_real(domain % dminfo, xMin, xMinGlobal)
      CALL mpas_dmpar_max_real(domain % dminfo, xMax, xMaxGlobal)
      CALL mpas_dmpar_min_real(domain % dminfo, dcEdgeMin, dcEdgeMinGlobal)

      yMidGlobal = (yMinGlobal + yMaxGlobal)*0.5_RKIND
      xPerturbationMin = xMinGlobal + 4.0_RKIND*(xMaxGlobal - xMinGlobal)/6.0_RKIND
      xPerturbationMax = xMinGlobal + 5.0_RKIND*(xMaxGlobal - xMinGlobal)/6.0_RKIND
      IF (config_baroclinic_channel_use_distances) THEN
         perturbationWidth = config_baroclinic_channel_gradient_width_dist
      ELSE
         perturbationWidth = (yMaxGlobal - yMinGlobal) &
                             *config_baroclinic_channel_gradient_width_frac
      END IF

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

         CALL mpas_pool_get_array(meshPool, 'xCell', xCell)
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

         ! Set refBottomDepth and refZMid
         DO k = 1, nVertLevels
            refBottomDepth(k) = config_baroclinic_channel_bottom_depth &
                                *interfaceLocations(k+1)
            refZMid(k) = -0.5_RKIND*(interfaceLocations(k+1) + interfaceLocations(k)) &
                          *config_baroclinic_channel_bottom_depth            
         END DO

         ! Set vertCoordMovementWeights
         vertCoordMovementWeights(:) = 1.0_RKIND

         DO iCell = 1, nCellsSolve
         
            ! Determine cutoff location for large sin wave
            yOffset = perturbationWidth*SIN(6.0_RKIND*pii*(xCell(iCell) - xMinGlobal) &
                      /(xMaxGlobal - xMinGlobal))

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
                  temperature = config_baroclinic_channel_bottom_temperature &
                                + (config_baroclinic_channel_surface_temperature &
                                   - config_baroclinic_channel_bottom_temperature) &
                                *((refZMid(k) + refBottomDepth(nVertLevels)) &
                                /refBottomDepth(nVertLevels))
                  activeTracers(idx, k, iCell) = temperature
               END DO

               IF (yCell(iCell) < yMidGlobal - yOffset) THEN
               ! If cell is in the southern half, outside the sin width, subtract 
               ! temperature difference
                  activeTracers(idx, :, iCell) &
                  = activeTracers(idx, :, iCell) &
                    - config_baroclinic_channel_temperature_difference
               ELSE IF (yCell(iCell) >= yMidGlobal - yOffset .AND. &
                        yCell(iCell) < yMidGlobal - yOffset + perturbationWidth) THEN
                  activeTracers(idx, :, iCell) &
                  = activeTracers(idx, :, iCell) &
                    - config_baroclinic_channel_temperature_difference &
                      *(1.0_RKIND - (yCell(iCell) - ((yMaxGlobal + yMinGlobal) &
                      *0.5_RKIND - yOffset))/perturbationWidth)
               END IF
   
               ! Determine yOffset for 3rd crest in sin wave.
               yOffset = 0.5_RKIND*perturbationWidth*SIN(pii*(xCell(iCell) &
                         - xPerturbationMin)/(xPerturbationMax - xPerturbationMin))
   
               IF (yCell(iCell) >= yMidGlobal - yOffset - 0.5_RKIND*perturbationWidth &
                   .AND. yCell(iCell) <= yMidGlobal - yOffset &
                                         + 0.5_RKIND*perturbationWidth .AND. &
                  xCell(iCell) >= xPerturbationMin .AND. &
                  xCell(iCell) <= xPerturbationMax) THEN
                  DO k = 1, nVertLevels
                     activeTracers(idx, k, iCell) &
                     = activeTracers(idx, k, iCell) &
                       + 0.3_RKIND*(1.0_RKIND - ((yCell(iCell) - (yMidGlobal &
                                    - yOffset))/(0.5_RKIND*perturbationWidth)))
                  END DO
               END IF
   
               ! Set salinity
               idx = index_salinity
               activeTracers(idx, :, iCell) = config_baroclinic_channel_salinity
               
            END IF

            ! Set layerThickness and restingThickness
            DO k = 1, nVertLevels
               layerThickness(k, iCell) &
               = config_baroclinic_channel_bottom_depth*(interfaceLocations(k+1) &
                                                         - interfaceLocations(k))
               restingThickness(k, iCell) &
               = config_baroclinic_channel_bottom_depth*(interfaceLocations(k+1) &
                                                         - interfaceLocations(k))
            END DO
   
            ! Set bottomDepth
            bottomDepth(iCell) = config_baroclinic_channel_bottom_depth
   
            ! Set maxLevelCell
            maxLevelCell(iCell) = nVertLevels
            
         END DO
   
         ! Set Coriolis parameters
         fCell(:) = config_baroclinic_channel_coriolis_parameter
         fEdge(:) = config_baroclinic_channel_coriolis_parameter
         fVertex(:) = config_baroclinic_channel_coriolis_parameter
   
         block_ptr => block_ptr % next
         
      END DO

      DEALLOCATE (interfaceLocations)

   END SUBROUTINE ocn_init_setup_baroclinic_channel
   !------------------------------------------------------------------------------------!
   !                                                                                    !
   !------------------------------------------------------------------------------------!
   !
   !  routine ocn_init_validate_baroclinic_channel
   !
   !> \brief   Validation for baroclinic channel test case
   !> \author  Doug Jacobsen
   !> \date    02/20/2014
   !> \details
   !> This routine validates the configuration options for the baroclinic channel test          
   !  case.
   !
   SUBROUTINE ocn_init_validate_baroclinic_channel(configPool, packagePool, iocontext, &
                                                   iErr)
   
      IMPLICIT NONE
   
      TYPE(mpas_pool_type), INTENT(INOUT) :: configPool, packagePool
      TYPE(mpas_io_context_type), INTENT(INOUT) :: iocontext
      INTEGER, INTENT(OUT) :: iErr
      
      ! Local variables
      CHARACTER(LEN=StrKIND), POINTER :: config_init_configuration
      INTEGER, POINTER :: config_vert_levels, config_baroclinic_channel_vert_levels

      iErr = 0

      CALL mpas_pool_get_config(configPool, 'config_init_configuration', &
                                config_init_configuration)

      IF (config_init_configuration /= TRIM('baroclinic_channel')) RETURN

      CALL mpas_pool_get_config(configPool, 'config_vert_levels', config_vert_levels)
      CALL mpas_pool_get_config(configPool, 'config_baroclinic_channel_vert_levels', &
                                config_baroclinic_channel_vert_levels)

      IF (config_vert_levels <= 0 .AND. config_baroclinic_channel_vert_levels > 0) THEN
         config_vert_levels = config_baroclinic_channel_vert_levels
      ELSE IF (config_vert_levels <= 0) THEN
         CALL mpas_log_write('Validation failed for baroclinic channel. Not given a &
                              &usable value for vertical levels.', MPAS_LOG_CRIT)
         iErr = 1
      END IF

   END SUBROUTINE ocn_init_validate_baroclinic_channel

END MODULE ocn_init_baroclinic_channel