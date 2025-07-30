// SIMULATION PARAMETERS.
const FLIP_RATIO = 0;  // Ratio of FLIP to PIC in velocity transfer (i.e. 0.95 means 95% FLIP, 5% PIC)

const GRAVITY = -1100;
const dt = 1.0 / 120.0;
const numberOfDivergenceIterations = 50;
const overRelaxation = 1.9;

const NUMBER_OF_PARTICLES = 14000;
const PARTICLE_RADIUS = 3;

const CELL_SPACING = 8;
const HALF_CELL_SPACING = CELL_SPACING / 2;
const INVERSE_CELL_SPACING = 1 / CELL_SPACING;
const X_CELLS = 220;
const Y_CELLS = 100;
const TOTAL_CELLS = X_CELLS * Y_CELLS;

const CANVAS_WIDTH = X_CELLS * CELL_SPACING;
const CANVAS_HEIGHT = Y_CELLS * CELL_SPACING;

// Vertices for the velocity grids (note, these grids are one cell smaller than the simulation grid in x and y, they are shifted by half a cell spacing).
// The vertical velocity grid is shifted to the right, and the horizontal velocity grid is shifted up.
const X_VERTICES = X_CELLS;
const Y_VERTICES = Y_CELLS;
const TOTAL_VERTICES = TOTAL_CELLS;
let uGrid = new Float32Array(TOTAL_VERTICES);  // Note, these store the cell vertices of these grids, not the cells
let vGrid = new Float32Array(TOTAL_VERTICES);
let uGridWeights = new Float32Array(TOTAL_VERTICES);  
let vGridWeights = new Float32Array(TOTAL_VERTICES);

let uGridPrevious = new Float32Array(TOTAL_VERTICES);  // For FLIP
let vGridPrevious = new Float32Array(TOTAL_VERTICES);  // For FLIP

const SOLID_CELL = 0;
const FLUID_CELL = 1;
const EMPTY_CELL = 2;
let cellType = new Float32Array(TOTAL_CELLS); 
let solidCells = new Float32Array(TOTAL_CELLS);  // 0 = solid, 1 = not


let particles = [];


class Particle {
    constructor(x, y, vx, vy, color) {
        this.x = x;
        this.y = y;
        this.vx = vx;
        this.vy = vy;
        this.color = color;
    }
}



const SPATIAL_CELL_SPACING = 2.2 * PARTICLE_RADIUS;  // Spacing of the spatial grid cells for particle separation
const INVERSE_SPATIAL_CELL_SPACING = 1.0 / SPATIAL_CELL_SPACING;
const SPATIAL_X_CELLS = Math.ceil(CANVAS_WIDTH * INVERSE_SPATIAL_CELL_SPACING);
const SPATIAL_Y_CELLS = Math.ceil(CANVAS_HEIGHT * INVERSE_SPATIAL_CELL_SPACING);
const SPATIAL_TOTAL_CELLS = SPATIAL_X_CELLS * SPATIAL_Y_CELLS;

let particleCountPerCell = new Int32Array(SPATIAL_TOTAL_CELLS);  // Number of particles in each spatial cell
let partialSums = new Int32Array(SPATIAL_TOTAL_CELLS + 1)  // +1 for the guard
let cellParticleIndices = new Int32Array(NUMBER_OF_PARTICLES);


function separateParticles(numberOfIterations) {
    particleCountPerCell.fill(0)

    // Count the number of particles in each grid cell.
    for (let particle of particles) {
        let gridX = Math.floor(particle.x * INVERSE_SPATIAL_CELL_SPACING);
        let gridY = Math.floor(particle.y * INVERSE_SPATIAL_CELL_SPACING);
        let gridIndex = gridX + gridY * SPATIAL_X_CELLS;

        particleCountPerCell[gridIndex]++;
    }

    // Compute partial sums.
    let currentSum = 0;
    for (let i = 0; i < SPATIAL_TOTAL_CELLS; i++) {
        currentSum += particleCountPerCell[i];
        partialSums[i] = currentSum;
    }
    partialSums[SPATIAL_TOTAL_CELLS] = currentSum;  // guard

    // Store ordered particle indices.
    for (let [i, particle] of particles.entries()) {
        let gridX = Math.floor(particle.x * INVERSE_SPATIAL_CELL_SPACING);
        let gridY = Math.floor(particle.y * INVERSE_SPATIAL_CELL_SPACING);
        let gridIndex = gridX + gridY * SPATIAL_X_CELLS;
        
        // Move the insertion index back one step.
        partialSums[gridIndex]--;

        // Add the particle index at the insertion index.
        cellParticleIndices[partialSums[gridIndex]] = i;
    }

    let minDist = 2 * PARTICLE_RADIUS;
    let minDistSquared = minDist * minDist;
    for (let iteration = 0; iteration < numberOfIterations; iteration++) {
        for (let [i, particle] of particles.entries()) {
            // Find the cell it's in and the 8 direct neighbours

            // for each particle in the 9 cells compute the distance between them
            // and push overlapping particles apart. Remember to skip the
            // current particle in this loop

            let gridX = Math.floor(particle.x * INVERSE_SPATIAL_CELL_SPACING);
            let gridY = Math.floor(particle.y * INVERSE_SPATIAL_CELL_SPACING);

            let minGridX = Math.max(0, gridX - 1);
            let maxGridX = Math.min(SPATIAL_X_CELLS - 1, gridX + 1)
            let minGridY = Math.max(0, gridY - 1);
            let maxGridY = Math.min(SPATIAL_Y_CELLS - 1, gridY + 1)

            for (let currentGridX = minGridX; currentGridX < maxGridX; currentGridX++) {
                for (let currentGridY = minGridY; currentGridY < maxGridY; currentGridY++) {
                    let currentGridIndex = currentGridX + currentGridY * SPATIAL_X_CELLS;

                    let startIndex = partialSums[currentGridIndex];
                    let endIndex = partialSums[currentGridIndex + 1];

                    for (let j = startIndex; j < endIndex; j++) {
                        let particleIndex = cellParticleIndices[j];

                        // Skip the current particle.
                        if (i === particleIndex) {
                            continue;
                        }

                        let neighbour = particles[particleIndex];

                        // Compute the displacement between the particles
                        let dx = neighbour.x - particle.x;
                        let dy = neighbour.y - particle.y;
                        let distSq = dx * dx + dy * dy;

                        if (distSq < minDistSquared && distSq > 0.0001) {
                            let dist = Math.sqrt(distSq);

                            // Normalized displacement vector
                            let nx = dx / dist;
                            let ny = dy / dist;

                            // Overlap amount
                            let overlap = 0.5 * (minDist - dist);

                            // Push each particle away from each other
                            particle.x -= nx * overlap;
                            particle.y -= ny * overlap;

                            neighbour.x += nx * overlap;
                            neighbour.y += ny * overlap;
                        }
                    }
                }
            }
        }
    }
}



/**
 * Sets the current type of each cell in the simulation grid (solid, fluid, or empty).
 */
function markCellTypes() {
    // Reset cell types and mark solid cells.
    for (let i = 0; i < TOTAL_CELLS; i++) {
        cellType[i] = solidCells[i] === 0 ? SOLID_CELL : EMPTY_CELL;
    }

    // Mark cells containing particles as fluid cells.
    for (let particle of particles) {
        let gridX = Math.floor(particle.x * INVERSE_CELL_SPACING);
        let gridY = Math.floor(particle.y * INVERSE_CELL_SPACING);
        let gridIndex = gridX + gridY * X_CELLS;

        if (cellType[gridIndex] === EMPTY_CELL) {
            cellType[gridIndex] = FLUID_CELL;
        }
    }
}

/**
 * Initializes the solidCells array to define which grid cells are solid (walls) and which are not.
 */
function setupSolidCells() {
    solidCells.fill(1, 0);

    // Mark the outermost layer of cells as solid.
    for (let gridX = 0; gridX < X_CELLS; gridX++) {
        for (let gridY = 0; gridY < Y_CELLS; gridY++) {
            if (gridX === 0 || gridY === 0 || gridX === X_CELLS - 1 || gridY === Y_CELLS - 1) {
                solidCells[gridX + gridY * X_CELLS] = 0;
            }
        }
    }
}

function solveIncompressibility(numberOfIterations, dt, overRelaxation) {
    // Store a copy of the current velocity grids for FLIP.
    uGridPrevious.set(uGrid);
    vGridPrevious.set(vGrid);

    for (let iteration = 0; iteration < numberOfIterations; iteration++) {

        // Loop over simulation grid, ignoring the outermost layer of cells (as these are walls).
        for (let gridX = 1; gridX < X_CELLS - 1; gridX++) {
            for (let gridY = 1; gridY < Y_CELLS - 1; gridY++) {

                // Skip non-fluid cells.
                if (cellType[gridX + gridY * X_CELLS] !== FLUID_CELL) {
                    continue;
                }

                // Get the array index of the current fluid cell.
                let centre = gridX + gridY * X_CELLS

                // Get the array indices of the surrounding cells.
                let left = centre - 1;
                let right = centre + 1;
                let bottom = centre - X_CELLS;
                let top = centre + X_CELLS;

                // Get the solid status of the surrounding cells (0 = solid).
                let leftSolid = solidCells[left];
                let rightSolid = solidCells[right];
                let bottomSolid = solidCells[bottom];
                let topSolid = solidCells[top];

                let surroundingSolid = leftSolid + rightSolid + bottomSolid + topSolid;

                // Calculate the divergence of the current fluid cell.
                let divergence = uGrid[right] - uGrid[centre] + vGrid[top] - vGrid[centre];
                let relaxedDivergence = divergence * overRelaxation;
                let divergenceCorrection = relaxedDivergence / surroundingSolid;

                // Update grid velocities to zero the current cells divergence.
                uGrid[centre] += leftSolid * divergenceCorrection;
                uGrid[right] -= rightSolid * divergenceCorrection;
                vGrid[centre] += bottomSolid * divergenceCorrection;
                vGrid[top] -= topSolid * divergenceCorrection;
            }
        }
    }
}

function transferVelocitiesToParticles() {
    for (let component of ['x', 'y']) {
        let { dx, dy, grid } = getComponentData(component);

        let prevF = component === 'x' ? uGridPrevious : vGridPrevious;  // FIX ME

        for (let particle of particles) {
            let { w1, w2, w3, w4, i1, i2, i3, i4 } = getWeightsAndIndices(particle, dx, dy);

            // FIX ME
            var offset = component === 'x' ? X_CELLS : 1;
            var valid0 = cellType[i1] != EMPTY_CELL || cellType[i1 - offset] != EMPTY_CELL ? 1.0 : 0.0;
            var valid1 = cellType[i2] != EMPTY_CELL || cellType[i2 - offset] != EMPTY_CELL ? 1.0 : 0.0;
            var valid2 = cellType[i3] != EMPTY_CELL || cellType[i3 - offset] != EMPTY_CELL ? 1.0 : 0.0;
            var valid3 = cellType[i4] != EMPTY_CELL || cellType[i4 - offset] != EMPTY_CELL ? 1.0 : 0.0;

            var d = valid0 * w1 + valid1 * w2 + valid2 * w3 + valid3 * w4;

            var v = component === 'x' ? particle.vx : particle.vy;

            if (d > 0.0) {
                var picV = (valid0 * w1 * grid[i1] + valid1 * w2 * grid[i2] + valid2 * w3 * grid[i3] + valid3 * w4 * grid[i4]) / d;
                var corr = (valid0 * w1 * (grid[i1] - prevF[i1]) + valid1 * w2 * (grid[i2] - prevF[i2])
								+ valid2 * w3 * (grid[i3] - prevF[i3]) + valid3 * w4 * (grid[i4] - prevF[i4])) / d;
				var flipV = v + corr;

                if (component === 'x') {
                    particle.vx = (1.0 - FLIP_RATIO) * picV + FLIP_RATIO * flipV;
                } else {
                    particle.vy = (1.0 - FLIP_RATIO) * picV + FLIP_RATIO * flipV;
                }
            }
            // END FIX ME
        }
    }
}

function transferVelocitiesToGrid() {
    // Reset velocities and weights.
    uGrid.fill(0, 0);
    vGrid.fill(0, 0);
    uGridWeights.fill(0, 0);
    vGridWeights.fill(0, 0);

    markCellTypes();

    for (let component of ['x', 'y']) {
        let { dx, dy, grid } = getComponentData(component);

        let gridWeights = component === 'x' ? uGridWeights : vGridWeights;

        for (let particle of particles) {
            let { w1, w2, w3, w4, i1, i2, i3, i4 } = getWeightsAndIndices(particle, dx, dy);
            
            let particleVelocity = component === 'x' ? particle.vx : particle.vy;

            // Update the sum of weighted grid velocities at the four corners of the grid cell.
            grid[i1] += w1 * particleVelocity;
            grid[i2] += w2 * particleVelocity;
            grid[i3] += w3 * particleVelocity;
            grid[i4] += w4 * particleVelocity;

            // Update the sum of grid weights at the four corners of the grid cell.
            gridWeights[i1] += w1;
            gridWeights[i2] += w2;
            gridWeights[i3] += w3;
            gridWeights[i4] += w4;
        }

        for (let i = 0; i < TOTAL_VERTICES; i++) {
            if (gridWeights[i] > 0) {
                // Normalise the sum of weighted grid velocities by the sum of weights at each grid cell.
                grid[i] /= gridWeights[i];
            }
        }
    }
}

function getComponentData(component) {
    return {
        dx: component === 'x' ? 0 : HALF_CELL_SPACING,
        dy: component === 'x' ? HALF_CELL_SPACING : 0,
        grid: component === 'x' ? uGrid : vGrid
    }
}

function getWeightsAndIndices(particle, dx, dy) {
    // Transform particle position into the given velocity grid’s coordinate system.
    let x = particle.x - dx;
    let y = particle.y - dy;

    // Find the grid cell axis indices for the cell containing the current particle.
    let gridX = Math.floor(x * INVERSE_CELL_SPACING);
    let gridY = Math.floor(y * INVERSE_CELL_SPACING);

    // Remainders of the particle's position in the grid cell.
    let deltaX = x - gridX * CELL_SPACING;
    let deltaY = y - gridY * CELL_SPACING;

    // Bilinear interpolation weights and array indices (w1/i1 = bottom-left corner, counter-clockwise).
    return {
        w1: (1 - deltaX * INVERSE_CELL_SPACING) * (1 - deltaY * INVERSE_CELL_SPACING),
        w2: deltaX * INVERSE_CELL_SPACING * (1 - deltaY * INVERSE_CELL_SPACING),
        w3: (deltaX * INVERSE_CELL_SPACING) * (deltaY * INVERSE_CELL_SPACING),
        w4: (1 - deltaX * INVERSE_CELL_SPACING) * (deltaY * INVERSE_CELL_SPACING),

        i1: gridX + gridY * X_VERTICES,
        i2: (gridX + 1) + gridY * X_VERTICES,
        i3: (gridX + 1) + (gridY + 1) * X_VERTICES,
        i4: gridX + (gridY + 1) * X_VERTICES
    }
}


/**
 * Detects and resolves collisions between particles and walls in the simulation 
 * by clamping positions to stay within bounds and zeroing velocity along the 
 * collision axis. Note that the outermost layer of grid cells in the simulation
 * are the walls considered here.
 */
function handleWallCollisions() {
    for (let particle of particles) {
        // Define the simulation bounds.
        let minX = PARTICLE_RADIUS + CELL_SPACING;
        let maxX = CANVAS_WIDTH - PARTICLE_RADIUS - CELL_SPACING;
        let minY = PARTICLE_RADIUS + CELL_SPACING;
        let maxY = CANVAS_HEIGHT - PARTICLE_RADIUS - CELL_SPACING;

        // Collision with left wall.
        if (particle.x < minX) {
            particle.x = minX
            particle.vx = 0;
        }

        // Collision with right wall.
        if (particle.x > maxX) {
            particle.x = maxX;
            particle.vx = 0;
        }

        // Collision with bottom wall.
        if (particle.y < minY) {
            particle.y = minY;
            particle.vy = 0;
        }

        // Collision with top wall.
        if (particle.y > maxY) {
            particle.y = maxY; 
            particle.vy = 0;
        }
    }
}

/**
 * Updates particle positions and velocities using Euler integration.
 *
 * @param {number} dt - Time step for integration.
 * @param {number} gravity - Vertical acceleration due to gravity.
 */
function integrateParticles(dt, gravity) {
    for (let particle of particles) {
        // Apply gravity.
        particle.vy += gravity * dt;

        // Update particle position based on velocity.
        particle.x += particle.vx * dt;
        particle.y += particle.vy * dt;
    }
}

/**
 * Initialises the particle array by placing particles in a staggered grid
 * on the lower-left half of the canvas. 
 */
function createParticles() {
    const spacing = 2 * PARTICLE_RADIUS * 1.1;
    const particlesPerRow = Math.floor(0.5 * CANVAS_WIDTH / spacing);
    const rows = Math.ceil(NUMBER_OF_PARTICLES / particlesPerRow);

    let count = 0;
    for (let row = 0; row < rows; row++) {
        const xOffset = (row % 2 === 0) ? 0 : spacing * 0.5;
        for (let col = 0; col < particlesPerRow; col++) {
            if (count >= NUMBER_OF_PARTICLES) break;
            const x = col * spacing + PARTICLE_RADIUS + CELL_SPACING + xOffset;
            const y = row * spacing + PARTICLE_RADIUS + CELL_SPACING;
            const vx = 0.0;
            const vy = 0.0;
            const color = [0.0, 0.0, 1.0, 1];
            particles.push(new Particle(x, y, vx, vy, color));
            count++;
        }
    }
}

/**
 * Packs each particle's position (x, y) and color (r, g, b, a) into a Float32Array.
 *
 * @returns {Float32Array} Flat buffer of particle data.
 */
function particlesToBuffer() {
    let bufferData = new Float32Array(NUMBER_OF_PARTICLES * 6);
    for (let i = 0; i < NUMBER_OF_PARTICLES; i++) {
        let particle = particles[i];
        let particlePositionNDC = simulationCoordinatesToNDC(particle.x, particle.y);

        bufferData[i * 6] = particlePositionNDC.x;
        bufferData[i * 6 + 1] = particlePositionNDC.y;
        bufferData[i * 6 + 2] = particle.color[0];
        bufferData[i * 6 + 3] = particle.color[1];
        bufferData[i * 6 + 4] = particle.color[2];
        bufferData[i * 6 + 5] = particle.color[3];
    }
    return bufferData;
}

/**
 * Converts simulation coordinates (origin bottom-left) to Normalised Device Coordinates 
 * (NDC), where the origin is the centre of the screen and axes range from -1 to 1.
 *
 * @param {number} x - X position in simulation space.
 * @param {number} y - Y position in simulation space.
 * @returns {{x: number, y: number}} An object containing x and y in NDC space.
 */
function simulationCoordinatesToNDC(x, y) {
    return {
        x: (x / CANVAS_WIDTH) * 2 - 1,
        y: (y / CANVAS_HEIGHT) * 2 - 1
    };
}









const vertexShaderSource = `#version 300 es
precision mediump float;

uniform float uPointSize;
in vec2 aPosition;
in vec4 aColor;

out vec4 vColor;

void main() {
    vColor = aColor;
    gl_PointSize = uPointSize;
    gl_Position = vec4(aPosition, 0.0, 1.0);
}`;

const fragmentShaderSource = `#version 300 es
precision mediump float;

in vec4 vColor;

out vec4 fragColor;

void main() {
    fragColor = vec4(vColor);

    // Convert square point into a circle
    vec2 circCoord = 2.0 * gl_PointCoord - 1.0;
    if (dot(circCoord, circCoord) > 1.0) {
        discard;
    }
}`;

const canvas = document.getElementById('canvas');
canvas.width = CANVAS_WIDTH;
canvas.height = CANVAS_HEIGHT;
const gl = canvas.getContext('webgl2');

const program = gl.createProgram();

const vertexShader = gl.createShader(gl.VERTEX_SHADER);
gl.shaderSource(vertexShader, vertexShaderSource);
gl.compileShader(vertexShader);
gl.attachShader(program, vertexShader);

const fragmentShader = gl.createShader(gl.FRAGMENT_SHADER);
gl.shaderSource(fragmentShader, fragmentShaderSource);
gl.compileShader(fragmentShader);
gl.attachShader(program, fragmentShader);

gl.linkProgram(program);

gl.useProgram(program);

createParticles();  // Create particles
setupSolidCells();  // Setup solid cells
const particlesBufferData = particlesToBuffer();

const uPointSizeLocation = gl.getUniformLocation(program, 'uPointSize');
gl.uniform1f(uPointSizeLocation, 2 * PARTICLE_RADIUS);

const aPositionLocation = gl.getAttribLocation(program, 'aPosition');
const aColorLocation = gl.getAttribLocation(program, 'aColor');

gl.enableVertexAttribArray(aPositionLocation);
gl.enableVertexAttribArray(aColorLocation);

const particlesBuffer = gl.createBuffer();
gl.bindBuffer(gl.ARRAY_BUFFER, particlesBuffer);
gl.bufferData(gl.ARRAY_BUFFER, particlesBufferData, gl.DYNAMIC_DRAW);

// Each particle has 6 components. The first two are the particle's position (x, y),
// the next four are the particle's RGBA color (r, g, b, a). Each component is a 32 
// bit float, thus 6 * 4 = 24 bytes per particle (recall 8 bits = 1 byte).
gl.vertexAttribPointer(aPositionLocation, 2, gl.FLOAT, false, 6 * 4, 0);
gl.vertexAttribPointer(aColorLocation, 4, gl.FLOAT, false, 6 * 4, 2 * 4);




function animate() {
    integrateParticles(dt, GRAVITY);
    separateParticles(2);
    handleWallCollisions();
    transferVelocitiesToGrid();
    solveIncompressibility(numberOfDivergenceIterations, dt, overRelaxation);
    transferVelocitiesToParticles();

    const bufferData = particlesToBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, particlesBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, bufferData, gl.DYNAMIC_DRAW);

    gl.clearColor(0, 0, 0, 1);  // Black background
    gl.clear(gl.COLOR_BUFFER_BIT);
    gl.drawArrays(gl.POINTS, 0, particles.length);

    requestAnimationFrame(animate);
}
requestAnimationFrame(animate);