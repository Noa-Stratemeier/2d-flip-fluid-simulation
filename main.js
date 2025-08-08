const SOLID_CELL = 0;
const FLUID_CELL = 1;
const EMPTY_CELL = 2;

class FlipFluidSimulation {
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

        this.spatialCellSize = 2.2 * this.particleRadius;
        this.inverseSpatialCellSize = 1.0 / this.spatialCellSize;
        this.spatialCellCountX = Math.ceil(this.width * this.inverseSpatialCellSize);
        this.spatialCellCountY = Math.ceil(this.height * this.inverseSpatialCellSize);
        this.spatialCellCount = this.spatialCellCountX * this.spatialCellCountY;
        this.particleCountPerCell = new Int32Array(this.spatialCellCount);
        this.partialSums = new Int32Array(this.spatialCellCount + 1)
        this.cellParticleIndices = new Int32Array(this.particleCount);
    }

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

                            // Skip self-collision: don't resolve collisions between the particle and itself.
                            if (i === iNeighbour) continue;

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

        // Simulation bounds (outermost layer of grid cells are walls).
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
                    let weightedDeltaSum = 
                        w1Valid * (velocityGrid[i1] - velocityGridPrevious[i1]) + 
                        w2Valid * (velocityGrid[i2] - velocityGridPrevious[i2]) + 
                        w3Valid * (velocityGrid[i3] - velocityGridPrevious[i3]) + 
                        w4Valid * (velocityGrid[i4] - velocityGridPrevious[i4]);

                    let picVelocity = weightedVelocitySum / weightSum;
                    let flipVelocity = v + (weightedDeltaSum / weightSum);

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
                let falloff = Math.pow(1 - distance / radius, 2);

                let interactionForceX = (ndx * strength - velocities[xi]) * falloff;
                let interactionForceY = (ndy * strength - velocities[yi]) * falloff;
                
                velocities[xi] += interactionForceX;
                velocities[yi] += interactionForceY;
            }
        }
    }

    updateParticleColours(baseColour, lowDensityColour, fadeSpeed = 0.01, lowDensityThreshold = 0.7) {
        let positions = this.particlePositions;
        let particleColours = this.particleColours;
        let densityGrid = this.densityGrid;
        let inverseCellSize = this.inverseCellSize;
        let cellCountX = this.cellCountX;
        let restDensity = this.restDensity;

        for (let i = 0; i < this.particleCount; i++) {
            let xi = 2 * i;
            let yi = 2 * i + 1;

            let ri = 3 * i;
            let gi = 3 * i + 1;
            let bi = 3 * i + 2;

            let gridX = Math.floor(positions[xi] * inverseCellSize);
            let gridY = Math.floor(positions[yi] * inverseCellSize);
            let gridIndex = gridX + gridY * cellCountX;

            if (restDensity > 0.0) {
                let relativeDensity = densityGrid[gridIndex] / restDensity;
                if (relativeDensity < lowDensityThreshold) {
                    particleColours[ri] = lowDensityColour[0];
                    particleColours[gi] = lowDensityColour[1];
                    particleColours[bi] = lowDensityColour[2];
                }
            }

            // Fade particles towards the given base colour.
            particleColours[ri] = FlipFluidSimulation.fadeTowards(particleColours[ri], baseColour[0], fadeSpeed);
            particleColours[gi] = FlipFluidSimulation.fadeTowards(particleColours[gi], baseColour[1], fadeSpeed);
            particleColours[bi] = FlipFluidSimulation.fadeTowards(particleColours[bi], baseColour[2], fadeSpeed);
        }
    }

    /**
     * Colours particles using a rainbow gradient based on their speed.
     *
     * @param {number} maxSpeed - Speed mapped to the highest colour value.
     */
    updateParticleColoursBySpeed(maxSpeed = 500) {
        let velocities = this.particleVelocities;
        let colours = this.particleColours;

        let maxSpeedSquared = maxSpeed * maxSpeed;

        for (let i = 0; i < this.particleCount; i++) {
            let vx = velocities[2 * i];
            let vy = velocities[2 * i + 1];
            let speedSquared = vx * vx + vy * vy;

            // Normalised speed.
            let t = Math.sqrt(speedSquared / maxSpeedSquared);
            t = Math.min(1.0, Math.max(0.0, t));  // Clamp

            let [r, g, b] = FlipFluidSimulation.rainbowColourMap(t);

            colours[3 * i + 0] = FlipFluidSimulation.fadeTowards(colours[3 * i + 0], r, 0.05);
            colours[3 * i + 1] = FlipFluidSimulation.fadeTowards(colours[3 * i + 1], g, 0.05);
            colours[3 * i + 2] = FlipFluidSimulation.fadeTowards(colours[3 * i + 2], b, 0.05);
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

    ///////////////////////////////////////////////////////////////////////////
    // HELPER FUNCTIONS
    ///////////////////////////////////////////////////////////////////////////

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
     * Maps a value in [0, 1] to a rainbow colour (blue → cyan → green → yellow → red).
     * Based on MATLAB-style "jet" colormap.
     * 
     * @param {number} t - Normalised value between 0.0 and 1.0.
     * @returns {Array<number>} RGB colour in [0, 1].
     */
    static rainbowColourMap(t) {
        let r = Math.min(Math.max(1.5 - Math.abs(4.0 * t - 3.0), 0.0), 1.0);
        let g = Math.min(Math.max(1.5 - Math.abs(4.0 * t - 2.0), 0.0), 1.0);
        let b = Math.min(Math.max(1.5 - Math.abs(4.0 * t - 1.0), 0.0), 1.0);

        return [r, g, b];
    }
}









let scene = {   
    // Tank.
    tankWidth: window.innerWidth - 30,
    tankHeight: window.innerHeight - 30,
    resolution: 90,

    // Fluid.
    relativeFluidWidth: 0.6,
    relativeFluidHeight: 0.8,

    baseColour: [0.0, 0.0, 1.0],
    lowDensityColour: [1.0, 1.0, 1.0],
    fadeSpeed: 0.01,
    lowDensityThreshold: 0.7,

    particleDisplaySize: 2.0,

    // Interaction.
    cursorRepelRadius: 100,
    cursorRepelStrength: 1000,

    // Simulation.
    gravity: -1100,
    dt: 1.0 / 60.0,
    flipRatio: 0.95,
    projectionIterations: 50,
    particleSeparationIterations: 1,
    overRelaxation: 1.5,
    stiffnessConstant: 500.0,

    flipFluidSimulation: null,
};

function initialiseScene() {
    // Calculate simulation parameters.
    let cellSize = scene.tankHeight / scene.resolution;
    let cellCountX = Math.ceil(scene.tankWidth / cellSize)
    let cellCountY = scene.resolution;
    
    let particleRadius = 0.3 * cellSize;

    let particleSpacingX = 2 * particleRadius;
    let particleSpacingY = Math.sqrt(3) * particleRadius;

    let simulationWidth = cellCountX * cellSize;
    let simulationHeight = cellCountY * cellSize;

    let particleCountX = Math.floor((scene.relativeFluidWidth * simulationWidth - 2.0 * cellSize - 2.0 * particleRadius) / particleSpacingX);
	let particleCountY = Math.floor((scene.relativeFluidHeight * simulationHeight - 2.0 * cellSize - 2.0 * particleRadius) / particleSpacingY);
	let particleCount = particleCountX * particleCountY;
    
    // Create simulation.
    let f = scene.flipFluidSimulation = new FlipFluidSimulation(particleCount, particleRadius, cellCountX, cellCountY, cellSize);

    // Initialise particles.
    let i = 0;
    for (let px = 0; px < particleCountX; px++) {
        for (let py = 0; py < particleCountY; py++) {
            let xi = 2 * i;
            let yi = 2 * i + 1;

            let offset = py % 2 === 0 ? 0 : particleRadius;
            f.particlePositions[xi] = cellSize + particleRadius + px * particleSpacingX + offset;
            f.particlePositions[yi] = cellSize + particleRadius + py * particleSpacingY;

            i++;
        }
    }

    // Initialise particle colours.
    for (let i = 0; i < particleCount; i++) {
        let ri = i * 3;
        let gi = i * 3 + 1;
        let bi = i * 3 + 2;

        f.particleColours[ri] = scene.baseColour[0];
        f.particleColours[gi] = scene.baseColour[1];
        f.particleColours[bi] = scene.baseColour[2];
    }

    // Initialise solid cells.
    for (let gridX = 0; gridX < cellCountX; gridX++) {
        for (let gridY = 0; gridY < cellCountY; gridY++) {
            if (gridX === 0 || gridY === 0 || gridX === cellCountX - 1 || gridY === cellCountY - 1) {
                f.solidCells[gridX + gridY * cellCountX] = 0.0;  // Solid.
            } else {
                f.solidCells[gridX + gridY * cellCountX] = 1.0;  // Not solid.
            }
        }  
    }
}










initialiseScene();

// Define canvas and webgl2 context.
let canvas = document.getElementById('canvas');
canvas.width = scene.flipFluidSimulation.width;
canvas.height = scene.flipFluidSimulation.height;
let gl = canvas.getContext('webgl2');






// BUFFERS.
let particlePositionsBuffer;
let particleColoursBuffer;






let vertexShaderSource = `#version 300 es
precision mediump float;

uniform vec2 uSimulationDomainSize;
uniform float uPointSize;
in vec2 aPosition;
in vec3 aColour;

out vec3 vColour;

void main() {
    vColour = aColour;

    vec2 positionNDC = (aPosition / uSimulationDomainSize) * 2.0 - 1.0;

    gl_Position = vec4(positionNDC, 0.0, 1.0);
    gl_PointSize = uPointSize;
}`;

let fragmentShaderSource = `#version 300 es
precision mediump float;

in vec3 vColour;

out vec4 fragColour;

void main() {
    fragColour = vec4(vColour, 1.0);

    // Convert square point into a circle.
    vec2 circleCoord = 2.0 * gl_PointCoord - 1.0;
    if (dot(circleCoord, circleCoord) > 1.0) {
        discard;
    }
}`;


function createShaderProgram(gl, vsSource, fsSource) {
    let program = gl.createProgram();

    let vertexShader = gl.createShader(gl.VERTEX_SHADER);
    gl.shaderSource(vertexShader, vsSource);
    gl.compileShader(vertexShader);
    gl.attachShader(program, vertexShader);

    let fragmentShader = gl.createShader(gl.FRAGMENT_SHADER);
    gl.shaderSource(fragmentShader, fsSource);
    gl.compileShader(fragmentShader);
    gl.attachShader(program, fragmentShader);

    gl.linkProgram(program);

    return program;
}


function initialiseGL() {
    let program = createShaderProgram(gl, vertexShaderSource, fragmentShaderSource);
    gl.useProgram(program);

    // Get and set uniform locations.
    let uPointSizeLocation = gl.getUniformLocation(program, 'uPointSize');
    let uSimulationDomainSizeLocation = gl.getUniformLocation(program, 'uSimulationDomainSize');

    gl.uniform1f(uPointSizeLocation, 2.0 * scene.flipFluidSimulation.particleRadius * scene.particleDisplaySize);
    gl.uniform2f(uSimulationDomainSizeLocation, scene.flipFluidSimulation.width, scene.flipFluidSimulation.height);

    // Attribute setup.
    let aPositionLocation = gl.getAttribLocation(program, 'aPosition');
    let aColourLocation = gl.getAttribLocation(program, 'aColour');

    gl.enableVertexAttribArray(aPositionLocation);
    gl.enableVertexAttribArray(aColourLocation);

    particlePositionsBuffer = gl.createBuffer();
    particleColoursBuffer = gl.createBuffer();

    gl.bindBuffer(gl.ARRAY_BUFFER, particlePositionsBuffer);
    gl.vertexAttribPointer(aPositionLocation, 2, gl.FLOAT, false, 2 * 4, 0);

    gl.bindBuffer(gl.ARRAY_BUFFER, particleColoursBuffer);
    gl.vertexAttribPointer(aColourLocation, 3, gl.FLOAT, false, 3 * 4, 0);
}




initialiseGL();
function animate() {
    scene.flipFluidSimulation.stepSimulation(
        scene.dt, 
        scene.gravity, 
        scene.flipRatio, 
        scene.overRelaxation, 
        scene.particleSeparationIterations, 
        scene.projectionIterations, 
        scene.stiffnessConstant
    );

    // scene.flipFluidSimulation.updateParticleColours(
    //     scene.baseColour, 
    //     scene.lowDensityColour, 
    //     scene.fadeSpeed, 
    //     scene.lowDensityThreshold
    // );

    scene.flipFluidSimulation.updateParticleColoursBySpeed(1000);


    // Update positions buffer.
    gl.bindBuffer(gl.ARRAY_BUFFER, particlePositionsBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, scene.flipFluidSimulation.particlePositions, gl.DYNAMIC_DRAW);
    
    // Update colours buffer.
    gl.bindBuffer(gl.ARRAY_BUFFER, particleColoursBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, scene.flipFluidSimulation.particleColours, gl.STATIC_DRAW);

    // Clear and draw.
    gl.clearColor(0.0, 0.0, 0.0, 1.0);
    gl.clear(gl.COLOR_BUFFER_BIT);
    gl.drawArrays(gl.POINTS, 0, scene.flipFluidSimulation.particleCount);

    requestAnimationFrame(animate);
}
animate();









// Set up mouse interaction.
canvas.addEventListener('mousedown', handleMouseInteraction);
canvas.addEventListener('mousemove', handleMouseInteraction);

function handleMouseInteraction(e) {
    if (e.buttons !== 1) return;  // Left mouse button.
    
    // Get mouse position in simulation coordinates.
    let rect = canvas.getBoundingClientRect();
    let mouseX = e.clientX - rect.left;
    let mouseY = rect.height - (e.clientY - rect.top);  // Flip Y axis.
    
    // Push particles away from mouse.
    scene.flipFluidSimulation.repelParticles(mouseX, mouseY, scene.cursorRepelRadius, scene.cursorRepelStrength);
}


// Set up touch interaction (mobile).
canvas.addEventListener('touchstart', handleTouchInteraction, { passive: false });
canvas.addEventListener('touchmove', handleTouchInteraction, { passive: false });

function handleTouchInteraction(e) {
    e.preventDefault(); // Prevent scrolling

    if (e.touches.length === 0) return;

    let touch = e.touches[0];

    let rect = canvas.getBoundingClientRect();
    let touchX = touch.clientX - rect.left;
    let touchY = rect.height - (touch.clientY - rect.top); // Flip Y axis

    scene.flipFluidSimulation.repelParticles(touchX, touchY, 50, scene.cursorRepelStrength);
}

