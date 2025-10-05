const SOLID_CELL = 0;
const FLUID_CELL = 1;
const EMPTY_CELL = 2;

/**
 * Implements a 2D FLIP fluid solver.
 */
export default class FlipFluidSimulation {
    /**
     * Creates a new fluid simulation with given particle and grid parameters.
     *
     * @param {number} particleCount - Number of fluid particles.
     * @param {number} particleRadius - Interaction radius of each particle.
     * @param {number} cellCountX - Number of grid cells along the x-axis.
     * @param {number} cellCountY - Number of grid cells along the y-axis.
     * @param {number} cellSize - Physical size of each grid cell (width).
     */
    constructor(particleCount, particleRadius, cellCountX, cellCountY, cellSize) {
        this.particleCount = particleCount;
        this.particleRadius = particleRadius;
        this.particlePositions = new Float32Array(this.particleCount * 2);  // (x, y).
        this.particleVelocities = new Float32Array(this.particleCount * 2);  // (vx, vy).
        this.particleColours = new Float32Array(this.particleCount * 3);  // (r, g, b).

        this.cellCountX = cellCountX;
        this.cellCountY = cellCountY;
        this.cellCount = cellCountX * cellCountY;

        this.cellSize = cellSize;
        this.halfCellSize = cellSize / 2.0;
        this.inverseCellSize = 1.0 / cellSize;

        this.cellType = new Int32Array(this.cellCount); 
        this.solidCells = new Int32Array(this.cellCount);

        this.width = cellCountX * cellSize;
        this.height = cellCountY * cellSize;

        this.uGrid = new Float32Array(this.cellCount);
        this.vGrid = new Float32Array(this.cellCount);
        this.uGridPrevious = new Float32Array(this.cellCount);
        this.vGridPrevious = new Float32Array(this.cellCount);
        this.uGridWeights = new Float32Array(this.cellCount);
        this.vGridWeights = new Float32Array(this.cellCount);

        this.restDensity = 0.0;
        this.densityGrid = new Float32Array(this.cellCount);

        this.spatialCellSize = 2.0 * this.particleRadius;
        this.inverseSpatialCellSize = 1.0 / this.spatialCellSize;
        this.spatialCellCountX = Math.ceil(this.width * this.inverseSpatialCellSize);
        this.spatialCellCountY = Math.ceil(this.height * this.inverseSpatialCellSize);
        this.spatialCellCount = this.spatialCellCountX * this.spatialCellCountY;
        this.particleCountPerCell = new Int32Array(this.spatialCellCount);
        this.partialSums = new Int32Array(this.spatialCellCount + 1)
        this.cellParticleIndices = new Int32Array(this.particleCount);
    }

    // -----------------------------------------------------------------------------
    // Main simulation methods.
    // -----------------------------------------------------------------------------

    /**
     * Updates particle positions and velocities using semi-implicit Euler integration.
     *
     * @param {number} dt - Time step size.
     * @param {number} gravity - Acceleration due to gravity.
     */
    integrateParticles(dt, gravity) {
        let positions = this.particlePositions;
        let velocities = this.particleVelocities;

        for (let i = 0; i < this.particleCount; i++) {
            let xi = 2 * i;
            let yi = 2 * i + 1;

            velocities[yi] += gravity * dt;

            positions[xi] += velocities[xi] * dt;
            positions[yi] += velocities[yi] * dt;
        }
    }

    /**
     * Resolves particle overlaps by iteratively separating nearby overlapping pairs.
     * Uses a spatial hash grid to check overlaps with neighbours in a 3×3 cell region.
     *
     * @param {number} numberOfIterations - Number of separation iterations to perform.
     */
    separateParticles(numberOfIterations) {
        let inverseSpatialCellSize = this.inverseSpatialCellSize;
        let spatialCellCountX = this.spatialCellCountX;
        let spatialCellCountY = this.spatialCellCountY;
        let particleCountPerCell = this.particleCountPerCell;
        let partialSums = this.partialSums;
        let cellParticleIndices = this.cellParticleIndices;
        let positions = this.particlePositions;

        particleCountPerCell.fill(0)

        // Count the number of particles in each grid cell.
        for (let i = 0; i < this.particleCount; i++) {
            let xi = 2 * i;
            let yi = 2 * i + 1;

            let gridX = Math.floor(positions[xi] * inverseSpatialCellSize);
            let gridY = Math.floor(positions[yi] * inverseSpatialCellSize);
            let gridIndex = gridX + gridY * spatialCellCountX;

            particleCountPerCell[gridIndex]++;
        }

        // Compute partial sums.
        let currentSum = 0;
        for (let i = 0; i < this.spatialCellCount; i++) {
            currentSum += particleCountPerCell[i];
            partialSums[i] = currentSum;
        }
        partialSums[this.spatialCellCount] = currentSum;  // Guard.

        // Store ordered particle indices.
        for (let i = 0; i < this.particleCount; i++) {
            let xi = 2 * i;
            let yi = 2 * i + 1;

            let gridX = Math.floor(positions[xi] * inverseSpatialCellSize);
            let gridY = Math.floor(positions[yi] * inverseSpatialCellSize);
            let gridIndex = gridX + gridY * spatialCellCountX;
            
            // Move the insertion index back one step.
            partialSums[gridIndex]--;

            // Add the particle index at the insertion index.
            cellParticleIndices[partialSums[gridIndex]] = i;
        }

        let minDistance = 2 * this.particleRadius;
        let minDistanceSquared = minDistance * minDistance;
        for (let iteration = 0; iteration < numberOfIterations; iteration++) {
            // For each particle, check and resolve overlaps with neighbors in adjacent spatial cells (3×3 region).
            for (let i = 0; i < this.particleCount; i++) {
                let xi = 2 * i;
                let yi = 2 * i + 1;
  
                let gridX = Math.floor(positions[xi] * inverseSpatialCellSize);
                let gridY = Math.floor(positions[yi] * inverseSpatialCellSize);

                let minGridX = Math.max(0, gridX - 1);
                let maxGridX = Math.min(spatialCellCountX - 1, gridX + 1)
                let minGridY = Math.max(0, gridY - 1);
                let maxGridY = Math.min(spatialCellCountY - 1, gridY + 1)
                
                for (let gridX = minGridX; gridX <= maxGridX; gridX++) {
                    for (let gridY = minGridY; gridY <= maxGridY; gridY++) {
                        let gridIndex = gridX + gridY * spatialCellCountX;

                        let startIndex = partialSums[gridIndex];
                        let endIndex = partialSums[gridIndex + 1];

                        for (let j = startIndex; j < endIndex; j++) {
                            let iNeighbour = cellParticleIndices[j];

                            let xiNeighbour = 2 * iNeighbour;
                            let yiNeighbour = 2 * iNeighbour + 1;

                            // Don't resolve collisions between the same two particles twice.
                            if (iNeighbour <= i) continue;

                            let dx = positions[xiNeighbour] - positions[xi];
                            let dy = positions[yiNeighbour] - positions[yi];
                            let distanceSquared = dx * dx + dy * dy;

                            // If particles overlap, push them apart.
                            if (distanceSquared < minDistanceSquared && distanceSquared != 0.0) {
                                let distance = Math.sqrt(distanceSquared);

                                let ndx = dx / distance;
                                let ndy = dy / distance;

                                let overlap = 0.5 * (minDistance - distance);

                                positions[xi] -= ndx * overlap;
                                positions[yi] -= ndy * overlap;
                                positions[xiNeighbour] += ndx * overlap;
                                positions[yiNeighbour] += ndy * overlap;
                            }
                        }
                    }
                }
            }
        }
    }

    /**
     * Clamps particles within the simulation bounds and zeroes their velocity on impact.
     */
    handleWallCollisions() {
        let positions = this.particlePositions;
        let velocities = this.particleVelocities;

        // Simulation bounds (outermost layer of grid cells are solid walls).
        let minX = this.particleRadius + this.cellSize;
        let maxX = this.width - this.particleRadius - this.cellSize;
        let minY = this.particleRadius + this.cellSize;
        let maxY = this.height - this.particleRadius - this.cellSize;

        for (let i = 0; i < this.particleCount; i++) {
            let xi = 2 * i;
            let yi = 2 * i + 1;
            
            // Collision with left wall.
            if (positions[xi] < minX) {
                positions[xi] = minX;
                velocities[xi] = 0;
            }

            // Collision with right wall.
            if (positions[xi] > maxX) {
                positions[xi] = maxX;
                velocities[xi] = 0;
            }

            // Collision with bottom wall.
            if (positions[yi] < minY) {
                positions[yi] = minY;
                velocities[yi] = 0;
            }

            // Collision with top wall.
            if (positions[yi] > maxY) {
                positions[yi] = maxY; 
                velocities[yi] = 0;
            }
        }
    }

    /**
     * Transfers particle velocities to the staggered MAC grid by using bilinear interpolation to sum each 
     * particle's weighted velocity to surrounding grid points and normalising by the sum of weights.
     */
    transferVelocitiesToGrid() {
        let positions = this.particlePositions;
        let velocities = this.particleVelocities;
        let cellSize = this.cellSize;
        let inverseCellSize = this.inverseCellSize;
        let cellCountX = this.cellCountX
        
        // Reset grid data.
        this.uGrid.fill(0.0);
        this.vGrid.fill(0.0);
        this.uGridWeights.fill(0.0);
        this.vGridWeights.fill(0.0);

        this.markCellTypes();

        // 0 = x-component (uGrid), 1 = y-component (vGrid).
        for (let component = 0; component < 2; component++) {
            let dx = component === 0 ? 0 : this.halfCellSize;
            let dy = component === 0 ? this.halfCellSize : 0;
            let velocityGrid = component === 0 ? this.uGrid : this.vGrid;
            let velocityGridWeights = component === 0 ? this.uGridWeights : this.vGridWeights;

            for (let i = 0; i < this.particleCount; i++) {
                let xi = 2 * i;
                let yi = 2 * i + 1;

                // Particle position in this component's grid coordinate system.
                let x = positions[xi] - dx;
                let y = positions[yi] - dy;

                // Cell position containing the current particle.
                let gridX = Math.floor(x * inverseCellSize);
                let gridY = Math.floor(y * inverseCellSize);

                // Local position within the cell.
                let deltaX = x - gridX * cellSize;
                let deltaY = y - gridY * cellSize;

                // Bilinear interpolation weights and their corresponding array indices 
                // (w1/i1 = bottom-left corner, increasing counter-clockwise).
                let w1 = (1 - deltaX * inverseCellSize) * (1 - deltaY * inverseCellSize);
                let w2 = deltaX * inverseCellSize * (1 - deltaY * inverseCellSize);
                let w3 = (deltaX * inverseCellSize) * (deltaY * inverseCellSize);
                let w4 = (1 - deltaX * inverseCellSize) * (deltaY * inverseCellSize);

                let i1 = gridX + gridY * cellCountX;
                let i2 = (gridX + 1) + gridY * cellCountX;
                let i3 = (gridX + 1) + (gridY + 1) * cellCountX;
                let i4 = gridX + (gridY + 1) * cellCountX;
                
                // Particle velocity for this axis.
                let v = velocities[xi + component];

                // Accumulate weighted velocity and weights.
                velocityGrid[i1] += w1 * v;    
                velocityGrid[i2] += w2 * v;    
                velocityGrid[i3] += w3 * v;    
                velocityGrid[i4] += w4 * v; 
                
                velocityGridWeights[i1] += w1;
                velocityGridWeights[i2] += w2;
                velocityGridWeights[i3] += w3;
                velocityGridWeights[i4] += w4;
            }

            // Normalise velocities at each grid vertex.
            for (let i = 0; i < this.cellCount; i++) {
                if (velocityGridWeights[i] > 0) {
                    velocityGrid[i] /= velocityGridWeights[i];
                }
            }
        }
    }

    /**
     * Sets the current type of each cell in the simulation grid.
     */
    markCellTypes() {
        let positions = this.particlePositions;
        let inverseCellSize = this.inverseCellSize;
        let cellCountX = this.cellCountX
        let cellType = this.cellType;
        let solidCells = this.solidCells;

        // Reset cell types and mark solid cells.
        for (let i = 0; i < this.cellCount; i++) {
            cellType[i] = solidCells[i] === 0 ? SOLID_CELL : EMPTY_CELL;
        }

        // Mark cells containing particles as fluid cells.
        for (let i = 0; i < this.particleCount; i++) {
            let xi = 2 * i;
            let yi = 2 * i + 1;

            let gridX = Math.floor(positions[xi] * inverseCellSize);
            let gridY = Math.floor(positions[yi] * inverseCellSize);
            let gridIndex = gridX + gridY * cellCountX;

            if (cellType[gridIndex] === EMPTY_CELL) {
                cellType[gridIndex] = FLUID_CELL;
            }
        }
    }

    /**
     * Transfers velocities from the MAC grid back to particles using bilinear interpolation.
     * Blends PIC and FLIP velocities based on the given flipRatio to update particle velocities.
     *
     * @param {number} flipRatio - Blending factor between PIC (0.0) and FLIP (1.0) velocity updates.
     */
    transferVelocitiesToParticles(flipRatio) {
        let positions = this.particlePositions;
        let velocities = this.particleVelocities;
        let cellSize = this.cellSize;
        let inverseCellSize = this.inverseCellSize;
        let cellCountX = this.cellCountX;
        let cellType = this.cellType;

        for (let component = 0; component < 2; component++) {
            let dx = component === 0 ? 0 : this.halfCellSize;
            let dy = component === 0 ? this.halfCellSize : 0;
            let velocityGrid = component === 0 ? this.uGrid : this.vGrid;
            let velocityGridPrevious = component === 0 ? this.uGridPrevious : this.vGridPrevious;

            for (let i = 0; i < this.particleCount; i++) {
                let xi = 2 * i;
                let yi = 2 * i + 1;

                let x = positions[xi] - dx;
                let y = positions[yi] - dy;

                let gridX = Math.floor(x * inverseCellSize);
                let gridY = Math.floor(y * inverseCellSize);

                let deltaX = x - gridX * cellSize;
                let deltaY = y - gridY * cellSize;

                let w1 = (1 - deltaX * inverseCellSize) * (1 - deltaY * inverseCellSize);
                let w2 = deltaX * inverseCellSize * (1 - deltaY * inverseCellSize);
                let w3 = (deltaX * inverseCellSize) * (deltaY * inverseCellSize);
                let w4 = (1 - deltaX * inverseCellSize) * (deltaY * inverseCellSize);

                let i1 = gridX + gridY * cellCountX;
                let i2 = (gridX + 1) + gridY * cellCountX;
                let i3 = (gridX + 1) + (gridY + 1) * cellCountX;
                let i4 = gridX + (gridY + 1) * cellCountX;

                let v = velocities[xi + component];

                // Validate each weight based on whether the face it's stored at is adjacent to a fluid cell.
                let offset = component === 0 ? 1 : cellCountX;
                let w1Valid = cellType[i1] === FLUID_CELL || cellType[i1 - offset] === FLUID_CELL ? w1 : 0.0;
                let w2Valid = cellType[i2] === FLUID_CELL || cellType[i2 - offset] === FLUID_CELL ? w2 : 0.0;
                let w3Valid = cellType[i3] === FLUID_CELL || cellType[i3 - offset] === FLUID_CELL ? w3 : 0.0;
                let w4Valid = cellType[i4] === FLUID_CELL || cellType[i4 - offset] === FLUID_CELL ? w4 : 0.0;

                let weightSum = w1Valid + w2Valid + w3Valid + w4Valid;
                if (weightSum > 0.0) {
                    let weightedVelocitySum = w1Valid * velocityGrid[i1] + w2Valid * velocityGrid[i2] + w3Valid * velocityGrid[i3] + w4Valid * velocityGrid[i4];
                    let weightedDeltaVelocitySum = 
                        w1Valid * (velocityGrid[i1] - velocityGridPrevious[i1]) + 
                        w2Valid * (velocityGrid[i2] - velocityGridPrevious[i2]) + 
                        w3Valid * (velocityGrid[i3] - velocityGridPrevious[i3]) + 
                        w4Valid * (velocityGrid[i4] - velocityGridPrevious[i4]);

                    let picVelocity = weightedVelocitySum / weightSum;
                    let flipVelocity = v + (weightedDeltaVelocitySum / weightSum);

                    velocities[xi + component] = (1.0 - flipRatio) * picVelocity + flipRatio * flipVelocity;
                }  
            }
        }
    }

    /**
     * Iteratively adjusts grid velocities to enforce fluid incompressibility by zeroing divergence.
     * Applies over-relaxation and compensates for velocity drift by reducing divergence in dense regions.
     *
     * @param {number} numberOfIterations - Number of solver iterations to perform.
     * @param {number} overRelaxation - Factor to accelerate convergence of grid incompressibility.
     * @param {number} stiffnessConstant - Strength of density-based divergence correction.
     */
    solveIncompressibility(numberOfIterations, overRelaxation, stiffnessConstant) {
        let cellCountX = this.cellCountX;
        let cellCountY = this.cellCountY;
        let restDensity = this.restDensity;
        let densityGrid = this.densityGrid;
        let uGrid = this.uGrid;
        let vGrid = this.vGrid;
        let cellType = this.cellType;
        let solidCells = this.solidCells;

        // Store a copy of the current velocity grids for FLIP.
        this.uGridPrevious.set(uGrid);
        this.vGridPrevious.set(vGrid);

        for (let iteration = 0; iteration < numberOfIterations; iteration++) {

            // Loop over simulation grid, ignoring the outermost layer of cells (as these are walls).
            for (let gridX = 1; gridX < cellCountX - 1; gridX++) {
                for (let gridY = 1; gridY < cellCountY - 1; gridY++) {
                    let gridIndex = gridX + gridY * cellCountX;

                    // Skip non-fluid cells.
                    if (cellType[gridIndex] !== FLUID_CELL) continue;

                    // Get the array index of the current fluid cell.
                    let centre = gridIndex;

                    // Get the array indices of the surrounding cells.
                    let left = centre - 1;
                    let right = centre + 1;
                    let bottom = centre - cellCountX;
                    let top = centre + cellCountX;

                    // Get the solid status of the surrounding cells (0 = solid).
                    let leftSolid = solidCells[left];
                    let rightSolid = solidCells[right];
                    let bottomSolid = solidCells[bottom];
                    let topSolid = solidCells[top];

                    let surroundingSolid = leftSolid + rightSolid + bottomSolid + topSolid;

                    // Calculate divergence of the cell and apply over relaxation.
                    let divergence = uGrid[right] - uGrid[centre] + vGrid[top] - vGrid[centre];
                    divergence *= overRelaxation;

                    // Reduce divergence in dense regions to compensate for velocity drift.
                    if (restDensity > 0.0) {
                        let compression = densityGrid[gridIndex] - restDensity;
                        if (compression > 0.0) {
                            divergence -= stiffnessConstant * compression;
                        }       
                    }

                    let divergenceCorrection = divergence / surroundingSolid;

                    // Update grid velocities to zero the current cells divergence.
                    uGrid[centre] += leftSolid * divergenceCorrection;
                    uGrid[right] -= rightSolid * divergenceCorrection;
                    vGrid[centre] += bottomSolid * divergenceCorrection;
                    vGrid[top] -= topSolid * divergenceCorrection;
                }
            }
        }
    }

    /**
     * Computes a smoothed particle density per grid cell using bilinear weighting.
     * Each particle contributes to its 4 nearest non-solid cells. If rest density
     * is unset, it's initialised as the average density of all fluid cells.
     */
    updateDensityGrid() {
        let positions = this.particlePositions;
        let cellSize = this.cellSize;
        let inverseCellSize = this.inverseCellSize;
        let cellCountX = this.cellCountX;
        let densityGrid = this.densityGrid;
        let cellType = this.cellType;
        
        densityGrid.fill(0.0);

        let dx = this.halfCellSize;
        let dy = this.halfCellSize;

        // Particles contribute to the 4 nearest cell centres using bilinear weights,
        // producing a smooth estimate of particle density per grid cell.
        for (let i = 0; i < this.particleCount; i++) {
            let xi = 2 * i;
            let yi = 2 * i + 1;

            let x = positions[xi] - dx;
            let y = positions[yi] - dy;

            let gridX = Math.floor(x * inverseCellSize);
            let gridY = Math.floor(y * inverseCellSize);

            let deltaX = x - gridX * cellSize;
            let deltaY = y - gridY * cellSize;

            let w1 = (1 - deltaX * inverseCellSize) * (1 - deltaY * inverseCellSize);
            let w2 = deltaX * inverseCellSize * (1 - deltaY * inverseCellSize);
            let w3 = (deltaX * inverseCellSize) * (deltaY * inverseCellSize);
            let w4 = (1 - deltaX * inverseCellSize) * (deltaY * inverseCellSize);

            let i1 = gridX + gridY * cellCountX;
            let i2 = (gridX + 1) + gridY * cellCountX;
            let i3 = (gridX + 1) + (gridY + 1) * cellCountX;
            let i4 = gridX + (gridY + 1) * cellCountX;

            if (cellType[i1] !== SOLID_CELL) densityGrid[i1] += w1;
            if (cellType[i2] !== SOLID_CELL) densityGrid[i2] += w2;
            if (cellType[i3] !== SOLID_CELL) densityGrid[i3] += w3;
            if (cellType[i4] !== SOLID_CELL) densityGrid[i4] += w4;
        }

        // If unset, calculate rest density as the average density per fluid cell.
        if (this.restDensity === 0.0) {
            let densitySum = 0.0;
            let fluidCellCount = 0;

            for (let i = 0; i < this.cellCount; i++) {
                if (cellType[i] === FLUID_CELL) {
                    densitySum += densityGrid[i];
                    fluidCellCount++;
                }
            }

            if (fluidCellCount > 0) {
                this.restDensity = densitySum / fluidCellCount;
            }
        }
    }

    stepSimulation(dt, gravity, flipRatio, overRelaxation, particleSeparationIterations, projectionIterations, stiffnessConstant) {
        this.integrateParticles(dt, gravity);
        this.separateParticles(particleSeparationIterations);
        this.handleWallCollisions();
        this.transferVelocitiesToGrid();
        this.markCellTypes();
        this.updateDensityGrid();
        this.solveIncompressibility(projectionIterations, overRelaxation, stiffnessConstant);
        this.transferVelocitiesToParticles(flipRatio);
    }

    // -----------------------------------------------------------------------------
    // Additional methods.
    // -----------------------------------------------------------------------------

    /**
     * Applies a circular repulsive force to particles within a given radius of a point.
     *
     * @param {number} x - X-coordinate of the point.
     * @param {number} y - Y-coordinate of the point.
     * @param {number} radius - Radius of influence.
     * @param {number} strength - Strength of the repulsive force.
     */
    repelParticles(x, y, radius, strength) {
        let positions = this.particlePositions;
        let velocities = this.particleVelocities;

        let radiusSquared = radius * radius;
        
        // Apply repulsive force to each particle within the radius of influence.
        for (let i = 0; i < this.particleCount; i++) {
            let xi = 2 * i;
            let yi = 2 * i + 1;
            
            // Vector from the point to the particle.
            let dx = positions[xi] - x;
            let dy = positions[yi] - y;

            let distanceSquared = dx * dx + dy * dy;
            
            if (distanceSquared < radiusSquared && distanceSquared != 0) {
                let distance = Math.sqrt(distanceSquared);

                let ndx = dx / distance;
                let ndy = dy / distance;

                // Quadratic falloff (stronger near the center, weaker at the edge).
                let falloff = 1 - distance / radius;

                let interactionForceX = (ndx * strength - velocities[xi]) * falloff;
                let interactionForceY = (ndy * strength - velocities[yi]) * falloff;
                
                velocities[xi] += interactionForceX;
                velocities[yi] += interactionForceY;
            }
        }
    }

    /**
     * Updates each particle’s colour by fading it toward a target determined by density and/or velocity. By default
     * particles move toward `baseColour`. If `colourLowDensityParticles` is enabled and the particle’s grid cell 
     * density falls below a threshold, it takes on `lowDensityColour`. Otherwise, if `colourParticlesBySpeed` is 
     * enabled, the particle is coloured using `speedColourMap` based on its normalised speed. Low-density colouring 
     * takes precedence over speed colouring when both are enabled.
     *
     * @param {Array<number>} baseColour - RGB colour for normal particles.
     * @param {Array<number>} lowDensityColour - RGB colour for low density particles.
     * @param {boolean} colourLowDensityParticles - Enable low density colouring.
     * @param {boolean} colourParticlesBySpeed - Enable speed-based colouring.
     * @param {function} speedColourMap - Function mapping normalised speed (0-1) to RGB colour.
     * @param {number} [baseColourFade] - Fade speed toward baseColour per update (0 = no change, 1 = instant).
     * @param {number} [lowDensityColourFade] - Fade speed toward lowDensityColour.
     * @param {number} [particleSpeedColourFade] - Fade speed toward speed colour.
     */
    updateParticleColours(
        baseColour, 
        lowDensityColour, 
        colourLowDensityParticles, 
        colourParticlesBySpeed, 
        speedColourMap, 
        baseColourFade = 0.01, 
        lowDensityColourFade = 0.5, 
        particleSpeedColourFade = 0.05,
    ) {
        let positions = this.particlePositions;
        let velocities = this.particleVelocities;
        let colours = this.particleColours;
        let densityGrid = this.densityGrid;
        let inverseCellSize = this.inverseCellSize;
        let cellCountX = this.cellCountX;
        let restDensity = this.restDensity;

        let lowDensityThreshold = 0.70;
        let maxSpeed = 800;
        let maxSpeedSquared = maxSpeed * maxSpeed;

        for (let i = 0; i < this.particleCount; i++) {
            let xi = 2 * i;
            let yi = 2 * i + 1;

            let ri = 3 * i;
            let gi = 3 * i + 1;
            let bi = 3 * i + 2;

            // Default to base colour and fade speed.
            let targetColour = baseColour;
            let fadeSpeed = baseColourFade;

            // Low density logic.
            let isLowDensity = false;
            if (colourLowDensityParticles && restDensity > 0.0) {
                let gridX = Math.floor(positions[xi] * inverseCellSize);
                let gridY = Math.floor(positions[yi] * inverseCellSize);
                let gridIndex = gridX + gridY * cellCountX;
                let relativeDensity = densityGrid[gridIndex] / restDensity;
                if (relativeDensity < lowDensityThreshold) {
                    targetColour = lowDensityColour;
                    isLowDensity = true;
                    fadeSpeed = lowDensityColourFade;
                }
            }

            // Speed logic (only if not low density).
            if (colourParticlesBySpeed && !isLowDensity) {
                let vx = velocities[xi]
                let vy = velocities[yi];
                let t = Math.sqrt((vx * vx + vy * vy) / maxSpeedSquared);
                t = Math.min(1.0, Math.max(0.0, t));
                targetColour = speedColourMap(t);
                fadeSpeed = particleSpeedColourFade;
            }

            // Fade current colour toward target.
            colours[ri] = FlipFluidSimulation.fadeTowards(colours[ri], targetColour[0], fadeSpeed);
            colours[gi] = FlipFluidSimulation.fadeTowards(colours[gi], targetColour[1], fadeSpeed);
            colours[bi] = FlipFluidSimulation.fadeTowards(colours[bi], targetColour[2], fadeSpeed);
        }
    }

    // -----------------------------------------------------------------------------
    // Helper methods.
    // -----------------------------------------------------------------------------

    /**
     * Moves a value one step toward a target at the given rate, without overshooting.
     *
     * @param {number} value - The current value.
     * @param {number} target - The value to move toward.
     * @param {number} rate - The maximum change allowed per step.
     * @returns {number} The updated value after applying one rate-limited step.
     */
    static fadeTowards(value, target, rate) {
        let delta = target - value;
        if (Math.abs(delta) < rate) return target;
        return value + Math.sign(delta) * rate;
    }

    /**
     * Linear helpers for custom gradients.
     */
    static _lerp(a, b, t) { return a + (b - a) * t; }
    static _lerp3(c0, c1, t) {
        return [
            FlipFluidSimulation._lerp(c0[0], c1[0], t),
            FlipFluidSimulation._lerp(c0[1], c1[1], t),
            FlipFluidSimulation._lerp(c0[2], c1[2], t),
        ];
    }

    /**
     * Build a gradient function from stops: [[t0,[r,g,b]],[t1,[r,g,b]],...]
     * t must be ascending in [0,1]. Returns (t)=>[r,g,b] with linear interpolation.
     */
    static makeGradient(stops) {
        return (t) => {
            t = Math.min(1, Math.max(0, t));
            for (let i = 1; i < stops.length; i++) {
                let [t1, c1] = stops[i];
                let [t0, c0] = stops[i - 1];
                if (t <= t1) {
                    let u = (t - t0) / Math.max(1e-8, (t1 - t0));
                    return FlipFluidSimulation._lerp3(c0, c1, u);
                }
            }
            return stops[stops.length - 1][1];
        };
    }

    /** 
     * Ready-to-use colour maps.
     */
    static colourMaps = {
        // Jet (blue → cyan → green → yellow → red).
        jet: FlipFluidSimulation.makeGradient([
            [0.00,  [0.00, 0.00, 0.50]],
            [0.125, [0.00, 0.00, 1.00]],
            [0.375, [0.00, 1.00, 1.00]],
            [0.625, [1.00, 1.00, 0.00]],
            [0.875, [1.00, 0.00, 0.00]],
            [1.00,  [0.50, 0.00, 0.00]],
        ]),

        // Fire (black → red → orange → yellow → white).
        fire: FlipFluidSimulation.makeGradient([
            [0.00, [0.00, 0.00, 0.00]],
            [0.30, [0.50, 0.00, 0.00]],
            [0.60, [1.00, 0.30, 0.00]],
            [0.80, [1.00, 0.80, 0.00]],
            [1.00, [1.00, 1.00, 1.00]],
        ]),

        // Ice (navy → blue → cyan → white).
        ice: FlipFluidSimulation.makeGradient([
            [0.00, [0.02, 0.05, 0.25]],
            [0.35, [0.05, 0.20, 0.75]],
            [0.70, [0.10, 0.80, 0.90]],
            [1.00, [1.00, 1.00, 1.00]],
        ]),

        // Greyscale (black → white).
        greyscale: (t) => [t, t, t],

        // Plasma (purple → magenta → orange → yellow).
        plasma: FlipFluidSimulation.makeGradient([
            [0.00, [0.05, 0.00, 0.35]],
            [0.25, [0.35, 0.00, 0.68]],
            [0.50, [0.69, 0.13, 0.68]],
            [0.75, [0.99, 0.57, 0.38]],
            [1.00, [0.99, 0.90, 0.15]],
        ]),

        // Magma (black → deep purple → red-orange → peach).
        magma: FlipFluidSimulation.makeGradient([
            [0.00, [0.00, 0.00, 0.00]],
            [0.20, [0.10, 0.05, 0.20]],
            [0.45, [0.38, 0.10, 0.40]],
            [0.75, [0.82, 0.30, 0.20]],
            [1.00, [0.99, 0.80, 0.60]],
        ]),

        // Coolwarm (blue → white → red).
        coolwarm: FlipFluidSimulation.makeGradient([
            [0.00, [0.23, 0.30, 0.75]],
            [0.50, [0.95, 0.95, 0.95]],
            [1.00, [0.80, 0.20, 0.19]],
        ]),

        // Viridis (dark blue → green → yellow).
        viridis: FlipFluidSimulation.makeGradient([
            [0.00, [0.267, 0.005, 0.329]],
            [0.15, [0.283, 0.141, 0.458]],
            [0.35, [0.254, 0.265, 0.530]],
            [0.55, [0.164, 0.471, 0.558]],
            [0.75, [0.128, 0.566, 0.551]],
            [1.00, [0.993, 0.906, 0.144]],
        ]),

        // Terrain (deep water → shallow → sand → grass → dense forest).
        terrain: FlipFluidSimulation.makeGradient([
            [0.00, [0.01, 0.09, 0.30]],
            [0.12, [0.00, 0.45, 0.60]],
            [0.20, [0.90, 0.84, 0.60]],
            [0.35, [0.40, 0.70, 0.30]],
            [0.55, [0.10, 0.45, 0.15]],
            [1.00, [0.00, 0.20, 0.00]],
        ]),

        // Ocean (deep blue → light turquoise).
        ocean: FlipFluidSimulation.makeGradient([
            [0.00, [0.00, 0.08, 0.15]],
            [0.25, [0.00, 0.22, 0.35]],
            [0.50, [0.00, 0.45, 0.55]],
            [0.70, [0.15, 0.75, 0.65]],
            [1.00, [0.70, 0.95, 0.85]]
        ]),

        // Flamingo (deep pink → light pink).
        flamingo: FlipFluidSimulation.makeGradient([
            [0.00, [0.40, 0.04, 0.20]],
            [0.25, [0.90, 0.20, 0.45]],
            [0.50, [1.00, 0.45, 0.60]],
            [0.75, [1.00, 0.70, 0.75]],
            [1.00, [1.00, 0.90, 0.88]],
        ]),

        // Frog (deep green → light green).
        frog: FlipFluidSimulation.makeGradient([
            [0.00, [0.02, 0.06, 0.02]],
            [0.20, [0.06, 0.20, 0.06]],
            [0.45, [0.10, 0.45, 0.10]],
            [0.70, [0.35, 0.80, 0.12]],
            [0.85, [0.70, 0.95, 0.20]],
            [1.00, [0.93, 0.98, 0.80]],
        ])
    };
}
