!
! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!  ocn_init_overflow
!
!> \brief MPAS ocean initialize case -- Overflow
!> \author Doug Jacobsen
!> \date   02/18/2014
!> \details
!> This module contains the routines for initializing the the overflow test case.
!
MODULE ocn_init_overflow

   USE mpas_kind_types
   USE mpas_io_units
   USE mpas_derived_types
   USE mpas_pool_routines
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

   PUBLIC :: ocn_init_setup_overflow, ocn_init_validate_overflow

   !------------------------------------------------------------------------------------!
   !
   ! Private module variables
   !
   !------------------------------------------------------------------------------------!

   CONTAINS
   !
   !  routine ocn_init_setup_overflow
   !
   !> \brief   Setup for overflow test case
   !> \author  Doug Jacobsen
   !> \date    02/18/2014
   !> \details
   !> This routine sets up the initial conditions for the overflow test case.
   !
   SUBROUTINE ocn_init_setup_overflow(domain, iErr)
   
      IMPLICIT NONE

      TYPE(domain_type), intent(inout) :: domain
      INTEGER, INTENT(OUT) :: iErr
      
      ! Local variables      
      REAL(KIND=RKIND) :: yMin, yMax, dcEdgeMin
      REAL(KIND=RKIND) :: yMinGlobal, yMaxGlobal, dcEdgeMinGlobal
      REAL(KIND=RKIND) :: plugWidth
      REAL(KIND=RKIND) :: slopeCenter, slopeWidth

      TYPE(block_type), POINTER :: block_ptr

      TYPE(mpas_pool_type), POINTER :: meshPool
      TYPE(mpas_pool_type), POINTER :: statePool
      TYPE(mpas_pool_type), POINTER :: verticalMeshPool
      TYPE(mpas_pool_type), POINTER :: tracersPool

      INTEGER :: iCell, k

      ! Define dimensions
      INTEGER, POINTER :: nCellsSolve, nEdgesSolve, nVertLevels, nVertLevelsP1
      INTEGER, POINTER :: index_temperature, index_salinity, index_tracer1

      ! Define arrays
      INTEGER, DIMENSION(:), POINTER :: maxLevelCell
      REAL(KIND=RKIND), DIMENSION(:), POINTER :: yCell, refBottomDepth, bottomDepth, &
                                                 vertCoordMovementWeights, dcEdge
      REAL(KIND=RKIND), DIMENSION(:,:), POINTER :: layerThickness, restingThickness
      REAL(KIND=RKIND), DIMENSION(:,:,:), POINTER :: activeTracers, debugTracers

      ! Define configs
      CHARACTER(LEN=StrKIND), POINTER :: config_init_configuration, &
                                         config_vertical_grid, config_overflow_layer_type
      LOGICAL, POINTER :: config_overflow_use_distances
      REAL(KIND=RKIND), POINTER :: config_overflow_plug_width_dist, &
                                   config_overflow_slope_center_dist, &
                                   config_overflow_slope_width_dist, &
                                   config_overflow_plug_width_frac, &
                                   config_overflow_slope_center_frac, &
                                   config_overflow_slope_width_frac, &
                                   config_overflow_bottom_depth, &
                                   config_overflow_ridge_depth, &
                                   config_overflow_plug_temperature, &
                                   config_overflow_domain_temperature, &
                                   config_overflow_salinity, &
                                   config_overflow_isopycnal_min_thickness

      REAL(KIND=RKIND), DIMENSION(:), POINTER :: interfaceLocations

      iErr = 0

      CALL mpas_pool_get_config(ocnConfigs, 'config_init_configuration', &
                                config_init_configuration)

      IF (config_init_configuration .ne. TRIM('overflow')) RETURN

      CALL mpas_pool_get_config(ocnConfigs, 'config_vertical_grid', config_vertical_grid)

      CALL mpas_pool_get_config(ocnConfigs, 'config_overflow_use_distances', &
                                config_overflow_use_distances)
      CALL mpas_pool_get_config(ocnConfigs, 'config_overflow_plug_width_dist', &
                                config_overflow_plug_width_dist)
      CALL mpas_pool_get_config(ocnConfigs, 'config_overflow_slope_center_dist', &
                                config_overflow_slope_center_dist)
      CALL mpas_pool_get_config(ocnConfigs, 'config_overflow_slope_width_dist', &
                                config_overflow_slope_width_dist)
      CALL mpas_pool_get_config(ocnConfigs, 'config_overflow_plug_width_frac', &
                                config_overflow_plug_width_frac)
      CALL mpas_pool_get_config(ocnConfigs, 'config_overflow_slope_center_frac', &
                                config_overflow_slope_center_frac)
      CALL mpas_pool_get_config(ocnConfigs, 'config_overflow_slope_width_frac', &
                                config_overflow_slope_width_frac)
      CALL mpas_pool_get_config(ocnConfigs, 'config_overflow_bottom_depth', &
                                config_overflow_bottom_depth)
      CALL mpas_pool_get_config(ocnConfigs, 'config_overflow_ridge_depth', &
                                config_overflow_ridge_depth)
      CALL mpas_pool_get_config(ocnConfigs, 'config_overflow_plug_temperature', &
                                config_overflow_plug_temperature)
      CALL mpas_pool_get_config(ocnConfigs, 'config_overflow_domain_temperature', &
                                config_overflow_domain_temperature)
      CALL mpas_pool_get_config(ocnConfigs, 'config_overflow_salinity', &
                                config_overflow_salinity)
      CALL mpas_pool_get_config(ocnConfigs, 'config_overflow_layer_type', &
                                config_overflow_layer_type)
      CALL mpas_pool_get_config(ocnConfigs, 'config_overflow_isopycnal_min_thickness', &
                                config_overflow_isopycnal_min_thickness)

      CALL mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
      CALL mpas_pool_get_dimension(meshPool, 'nVertLevelsP1', nVertLevelsP1)

      ALLOCATE (interfaceLocations(nVertLevelsP1))
      CALL ocn_generate_vertical_grid(config_vertical_grid, interfaceLocations)

      ! Initalize y values to large positive and negative values
      yMin = 1.0E10_RKIND
      yMax = -1.0E10_RKIND
      dcEdgeMin = 1.0E10_RKIND

      ! Determine local min and max y value.
      block_ptr => domain % blocklist
      DO WHILE (ASSOCIATED(block_ptr))
        CALL mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)

        CALL mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
        CALL mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)

        CALL mpas_pool_get_array(meshPool, 'yCell', yCell)
        CALL mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)

        yMin = MIN(yMin, MINVAL(yCell(1:nCellssolve)))
        yMax = MAX(yMax, MAXVAL(yCell(1:nCellssolve)))
        dcEdgeMin = MIN(dcEdgeMin, MINVAL(dcEdge(1:nEdgessolve)))

        block_ptr => block_ptr % next
      END DO

      ! Determine global min and max y value. This is so the domain can be split into 
      ! north and south.
      CALL mpas_dmpar_min_real(domain % dminfo, yMin, yMinGlobal)
      CALL mpas_dmpar_max_real(domain % dminfo, yMax, yMaxGlobal)
      CALL mpas_dmpar_min_real(domain % dminfo, dcEdgeMin, dcEdgeMinGlobal)

      IF (config_overflow_use_distances) THEN
         plugWidth = config_overflow_plug_width_dist
         slopeCenter = yMinGlobal + config_overflow_slope_center_dist
         slopeWidth = config_overflow_slope_width_dist
      ELSE
         plugWidth = (yMaxGlobal - yMinGlobal) * config_overflow_plug_width_frac
         slopeCenter = yMinGlobal + (yMaxGlobal - yMinGlobal) &
                                    *config_overflow_slope_center_frac
         slopeWidth = (yMaxGlobal - yMinGlobal) * config_overflow_slope_width_frac
      END IF

      block_ptr => domain % blocklist
      DO WHILE (ASSOCIATED(block_ptr))
         CALL mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
         CALL mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
         CALL mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', &
                                    verticalMeshPool)
         CALL mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

         CALL mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
         CALL mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

         CALL mpas_pool_get_dimension(tracersPool, 'index_temperature', &
                                      index_temperature)
         CALL mpas_pool_get_dimension(tracersPool, 'index_salinity', index_salinity)
         CALL mpas_pool_get_dimension(tracersPool, 'index_tracer1', index_tracer1)

         CALL mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
         CALL mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
         CALL mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
         CALL mpas_pool_get_array(meshPool, 'vertCoordMovementWeights', &
                                  vertCoordMovementWeights)

         CALL mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)
         CALL mpas_pool_get_array(tracersPool, 'debugTracers', debugTracers, 1)
         CALL mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)

         CALL mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness)

         CALL ocn_mark_north_boundary(meshPool, yMaxGlobal, dcEdgeMinGlobal, iErr)
         CALL ocn_mark_south_boundary(meshPool, yMinGlobal, dcEdgeMinGlobal, iErr)

         ! Set refBottomDepth, bottomDepth, and maxLevelCell
         DO k = 1, nVertLevels
            refBottomDepth(k) = config_overflow_bottom_depth*interfaceLocations(k+1)
         END DO

         DO iCell = 1, nCellsSolve
            ! From Mehmet Ilicak:
            ! depth=2000
            ! val1 = 500 is top of ridge
            ! h(i,j) = val1 + 0.5*(depth-val1) * (1.0+TANH((lon(i,j)-40000.0)/7000.0))
            bottomDepth(iCell) &
            = config_overflow_ridge_depth &
              + 0.5_RKIND*(config_overflow_bottom_depth - config_overflow_ridge_depth) &
                         *(1.0_RKIND+tanh((yCell(iCell) - slopeCenter)/slopeWidth))

            IF (TRIM(config_overflow_layer_type) == 'sigma' .OR. &
                TRIM(config_overflow_layer_type) == 'isopycnal') THEN
               maxLevelCell(iCell) = nVertLevels
            ELSE IF (TRIM(config_overflow_layer_type) == 'z-level') THEN
               maxLevelCell(iCell) = -1
               DO k = 1, nVertLevels
                  IF (bottomDepth(iCell) .le. refBottomDepth(k) .AND. &
                     maxLevelCell(iCell) == -1) THEN

                     maxLevelCell(iCell) = k
                  END IF
               END DO
            END IF
         END DO

         DO iCell = 1, nCellsSolve
            ! Set temperature
            IF (ASSOCIATED(activeTracers)) THEN
               IF (TRIM(config_overflow_layer_type) == 'sigma' .OR. &
                   TRIM(config_overflow_layer_type) == 'z-level') THEN
                  DO k = 1, maxLevelCell(iCell)
                     IF (yCell(iCell) < yMinGlobal + plugWidth) THEN
                        activeTracers(index_temperature, k, iCell) &
                        = config_overflow_plug_temperature
                     ELSE
                        activeTracers(index_temperature, k, iCell) &
                        = config_overflow_domain_temperature
                     END IF
                  END DO
               ELSE IF (TRIM(config_overflow_layer_type) == 'isopycnal') THEN
                  activeTracers(index_temperature, 1, :) &
                  = config_overflow_domain_temperature
                  activeTracers(index_temperature, 2:nVertLevels, :) &
                  = config_overflow_plug_temperature
               END IF
            END IF

            ! Set layerThickness and restingThickness
            IF (TRIM(config_overflow_layer_type) == 'z-level') THEN
               DO k = 1, maxLevelCell(iCell)
                  layerThickness(k, iCell) &
                  = config_overflow_bottom_depth &
                    *(interfaceLocations(k+1) - interfaceLocations(k))
                  restingThickness(k, iCell) = layerThickness(k, iCell)
               END DO
            ELSE IF (TRIM(config_overflow_layer_type) == 'sigma') THEN
               DO k = 1, nVertLevels
                  layerThickness(k, iCell) &
                  = bottomDepth(iCell)/REAL(nVertLevels,KIND=RKIND)
                  restingThickness(k, iCell) = layerThickness(k, iCell)
               END DO
            ELSE IF (TRIM(config_overflow_layer_type) == 'isopycnal') THEN
               ! Set layerThickness. Normally isopycnal overflow has only two layers.
               IF (yCell(iCell) < yMinGlobal + plugWidth) THEN
                  layerThickness(1, iCell) = config_overflow_isopycnal_min_thickness
                  layerThickness(2:nVertLevels, iCell) &
                  = bottomDepth(iCell) - config_overflow_isopycnal_min_thickness
                  restingThickness(:, iCell) = layerThickness(:, iCell)
               ELSE
                  layerThickness(1, iCell) &
                  = bottomDepth(iCell) - config_overflow_isopycnal_min_thickness
                  layerThickness(2:nVertLevels, iCell) &
                  = config_overflow_isopycnal_min_thickness
                  restingThickness(:, iCell) = layerThickness(:, iCell)
               END IF
            END IF

            ! Set salinity
            IF (ASSOCIATED(activeTracers)) THEN
               activeTracers(index_salinity, :, iCell) = config_overflow_salinity
            END IF

            ! Set debug tracer
            IF (ASSOCIATED(debugTracers)) THEN
               DO k = 1, nVertLevels
                  debugTracers(index_tracer1, k, iCell) = 1.0_RKIND
               END DO
            END IF

         END DO

         ! Set vertCoordMovementWeights
         vertCoordMovementWeights(:) = 1.0_RKIND

         block_ptr => block_ptr % next
      END DO

      DEALLOCATE (interfaceLocations)

   END SUBROUTINE ocn_init_setup_overflow
   !------------------------------------------------------------------------------------!
   !                                                                                    !
   !------------------------------------------------------------------------------------!
   !
   !  routine ocn_init_validate_overflow
   !
   !> \brief   Validation for overflow test case
   !> \author  Doug Jacobsen
   !> \date    02/20/2014
   !> \details
   !> This routine validates the configuration options for the overflow test case.
   !
   SUBROUTINE ocn_init_validate_overflow(configPool, packagePool, iocontext, iErr)
   
      IMPLICIT NONE 

      TYPE(mpas_pool_type), intent(inout) :: configPool, packagePool
      TYPE(mpas_io_context_type), intent(inout) :: iocontext

      INTEGER, INTENT(OUT) :: iErr

      CHARACTER(LEN=StrKIND), POINTER :: config_init_configuration
      INTEGER, POINTER :: config_overflow_vert_levels, config_vert_levels

      iErr = 0

      CALL mpas_pool_get_config(configPool, 'config_init_configuration', &
                                config_init_configuration)

      IF (config_init_configuration .ne. TRIM('overflow')) RETURN

      CALL mpas_pool_get_config(configPool, 'config_vert_levels', config_vert_levels)
      CALL mpas_pool_get_config(configPool, 'config_overflow_vert_levels', &
                                config_overflow_vert_levels)

      IF (config_vert_levels <= 0 .AND. config_overflow_vert_levels > 0) THEN
         config_vert_levels = config_overflow_vert_levels
      ELSE IF (config_vert_levels <= 0) THEN
         CALL mpas_log_write('Validation failed for overflow test case. Not given a &
                              &usable value for vertical levels.', MPAS_LOG_CRIT)
         iErr = 1
      END IF

   END SUBROUTINE ocn_init_validate_overflow

END MODULE ocn_init_overflow