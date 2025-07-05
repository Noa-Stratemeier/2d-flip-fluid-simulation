// SIMULATION PARAMETERS.
const NUMBER_OF_PARTICLES = 100;
const PARTICLE_RADIUS = 5;

const CELL_SPACING = 10;
const HALF_CELL_SPACING = CELL_SPACING / 2;
const INVERSE_CELL_SPACING = 1 / CELL_SPACING;
const X_CELLS = 40;
const Y_CELLS = 30;
const TOTAL_CELLS = X_CELLS * Y_CELLS;

const CANVAS_WIDTH = X_CELLS * CELL_SPACING;
const CANVAS_HEIGHT = Y_CELLS * CELL_SPACING;

// Vertices for the velocity grids (note, these grids are one cell smaller than the simulation grid in x and y, they are shifted by half a cell spacing).
// The vertical velocity grid is shifted to the right, and the horizontal velocity grid is shifted up.
X_VERTICES = X_CELLS;
Y_VERTICES = Y_CELLS;
TOTAL_VERTICES = TOTAL_CELLS;
let horizontalGridVelocities = new Float32Array(TOTAL_VERTICES);  // Note, these store the cell vertices of these grids, not the cells
let verticalGridVelocities = new Float32Array(TOTAL_VERTICES);
let horizontalGridWeights = new Float32Array(TOTAL_VERTICES);  
let verticalGridWeights = new Float32Array(TOTAL_VERTICES);

const SOLID_CELL = 1;
const FLUID_CELL = 2;
const EMPTY_CELL = 3;
let cellType = new Float32Array(TOTAL_CELLS); 
let solidCells = new Float32Array(TOTAL_CELLS);  // 1 = solid, 0 = not


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


function solveIncompressibility(numberOfIterations, dt, overRelaxation) {

    for (let iteration = 0; iteration < numberOfIterations; iteration++) {

        // Loop over simulation grid, ignoring the outermost layer of cells (as these are walls).
        for (let gridX = 1; gridX < X_CELLS - 1; gridX++) {
            for (let gridY = 1; gridY < Y_CELLS - 1; gridX++) {
                


            }

        }


    }

}

function transferVelocitiesToGrid() {
    // Reset velocities and weights
    horizontalGridVelocities.fill(0, 0);
    verticalGridVelocities.fill(0, 0);
    horizontalGridWeights.fill(0, 0);
    verticalGridWeights.fill(0, 0);

    let components = ['x', 'y'];

    for (let component of components) {
        let dx = component === 'x' ? 0 : HALF_CELL_SPACING;
        let dy = component === 'x' ? HALF_CELL_SPACING : 0;

        let grid = component === 'x' ? horizontalGridVelocities : verticalGridVelocities;
        let gridWeights = component === 'x' ? horizontalGridWeights : verticalGridWeights;

        for (let i = 0; i < NUMBER_OF_PARTICLES; i++) {
            let particle = particles[i];

            // Particle's position in the current velocity grid's coordinate system.
            let x = particle.x - dx;
            let y = particle.y - dy;

            // Find the grid cell indices for the cell containing the current particle.
            let gridX = Math.floor(x * INVERSE_CELL_SPACING);
            let gridY = Math.floor(y * INVERSE_CELL_SPACING);

            // Remainders of the particle's position in the grid cell.
            let deltaX = x - gridX * CELL_SPACING;
            let deltaY = y - gridY * CELL_SPACING;

            // Bilinear interpolation weights for the four corners of the grid cell.
            let w1 = (1 - deltaX * INVERSE_CELL_SPACING) * (1 - deltaY * INVERSE_CELL_SPACING);
            let w2 = deltaX * INVERSE_CELL_SPACING * (1 - deltaY * INVERSE_CELL_SPACING);
            let w3 = (deltaX * INVERSE_CELL_SPACING) * (deltaY * INVERSE_CELL_SPACING);
            let w4 = (1 - deltaX * INVERSE_CELL_SPACING) * (deltaY * INVERSE_CELL_SPACING);

            // Calculate the array indices of the four corners of the grid cell.
            let bottomLeft = gridX + gridY * X_VERTICES;
            let bottomRight = (gridX + 1) + gridY * X_VERTICES;
            let topLeft = gridX + (gridY + 1) * X_VERTICES;
            let topRight = (gridX + 1) + (gridY + 1) * X_VERTICES;

            let particleVelocity = component === 'x' ? particle.vx : particle.vy;

            // Update the sum of weighted grid velocities at the four corners of the grid cell.
            grid[bottomLeft] += w1 * particleVelocity;
            grid[bottomRight] += w2 * particleVelocity;
            grid[topRight] += w3 * particleVelocity;
            grid[topLeft] += w4 * particleVelocity;

            // Update the sum of grid weights at the four corners of the grid cell.
            gridWeights[bottomLeft] += w1;
            gridWeights[bottomRight] += w2;
            gridWeights[topRight] += w3;
            gridWeights[topLeft] += w4;
        }

        for (let i = 0; i < TOTAL_VERTICES; i++) {
            if (gridWeights[i] > 0) {
                // Normalise the sum of weighted grid velocities by the sum of weights at each grid cell.
                grid[i] /= gridWeights[i];
            }
        }
    }
}


/**
 * Detects and resolves collisions between particles and walls in the simulation 
 * by clamping positions to stay within bounds and zeroing velocity along the 
 * collision axis. Note that the outermost layer of grid cells in the simulation
 * are the walls considered here.
 */
function handleWallCollisions() {
    for (let i = 0; i < NUMBER_OF_PARTICLES; i++) {
        let particle = particles[i];

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

        // Collision with top wall
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
 * @param {number} gravity - Vertical acceleration due to gravity (negative = down).
 */
function integrateParticles(dt, gravity) {
    for (let i = 0; i < NUMBER_OF_PARTICLES; i++) {
        let particle = particles[i];

        // Apply gravity.
        particle.vy += gravity * dt;

        // Update particle position based on velocity.
        particle.x += particle.vx * dt;
        particle.y += particle.vy * dt;
    }
}

/**
 * Initialises the particle array by placing particles in a non-overlapping grid
 * on the lower-left half of the canvas. Each particle is initialised with zero 
 * velocity and a solid blue color.
 */
function createParticles() {
    const spacing = 2 * PARTICLE_RADIUS * 1.1;
    const particlesPerRow = Math.floor(0.5 * CANVAS_WIDTH / spacing);
    const rows = Math.ceil(NUMBER_OF_PARTICLES / particlesPerRow);

    let count = 0;
    for (let row = 0; row < rows; row++) {
        for (let col = 0; col < particlesPerRow; col++) {
            if (count >= NUMBER_OF_PARTICLES) break;
            const x = col * spacing + PARTICLE_RADIUS + CELL_SPACING;
            const y = row * spacing + PARTICLE_RADIUS + CELL_SPACING;
            const vx = 0.0;
            const vy = 0.0;
            const color = [0.0, 0.0, 1.0, 1.0];
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
        let p = particles[i];
        let pPositionNDC = simulationCoordinatesToNDC(p.x, p.y);

        bufferData[i * 6] = pPositionNDC.x;
        bufferData[i * 6 + 1] = pPositionNDC.y;
        bufferData[i * 6 + 2] = p.color[0];
        bufferData[i * 6 + 3] = p.color[1];
        bufferData[i * 6 + 4] = p.color[2];
        bufferData[i * 6 + 5] = p.color[3];
    }
    return bufferData;
}

/**
 * Converts simulation coordinates (origin bottom-left) to Normalised Device Coordinates 
 * (NDC), where the origin is the center of the screen and axes range from -1 to 1. WebGL 
 * expects coordinates in NDC space, so this conversion is necessary for rendering.
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
const particlesBufferData = particlesToBuffer();

const uPointSizeLocation = gl.getUniformLocation(program, 'uPointSize');
gl.uniform1f(uPointSizeLocation, 2 * PARTICLE_RADIUS);

const aPositionLocation = gl.getAttribLocation(program, 'aPosition');
const aColorLocation = gl.getAttribLocation(program, 'aColor');

gl.enableVertexAttribArray(aPositionLocation);
gl.enableVertexAttribArray(aColorLocation);

const particlesBuffer = gl.createBuffer();
gl.bindBuffer(gl.ARRAY_BUFFER, particlesBuffer);
gl.bufferData(gl.ARRAY_BUFFER, particlesBufferData, gl.STATIC_DRAW);

// Each particle has 6 components. The first two are the particle's position (x, y),
// the next four are the particle's RGBA color (r, g, b, a). Each component is a 32 
// bit float, thus 6 * 4 = 24 bytes per particle (recall 8 bits = 1 byte).
gl.vertexAttribPointer(aPositionLocation, 2, gl.FLOAT, false, 6 * 4, 0);
gl.vertexAttribPointer(aColorLocation, 4, gl.FLOAT, false, 6 * 4, 2 * 4);

gl.drawArrays(gl.POINTS, 0, NUMBER_OF_PARTICLES);








// const dt = 1 / 60;
// function animate() {
//     integrateParticles(dt, -9.81);
//     handleWallCollisions();

//     const bufferData = particlesToBuffer();
//     gl.bindBuffer(gl.ARRAY_BUFFER, particlesBuffer);
//     gl.bufferData(gl.ARRAY_BUFFER, bufferData, gl.STATIC_DRAW);

//     gl.clear(gl.COLOR_BUFFER_BIT);
//     gl.drawArrays(gl.POINTS, 0, particles.length);

//     requestAnimationFrame(animate);
// }
// requestAnimationFrame(animate);


