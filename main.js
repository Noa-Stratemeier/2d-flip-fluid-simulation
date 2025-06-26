// SIMULATION PARAMETERS
const NUMBER_OF_PARTICLES = 500;
const PARTICLE_DIAMETER = 10.0; 
const PARTICLE_RADIUS = PARTICLE_DIAMETER / 2.0;
const CANVAS_WIDTH = 800;
const CANVAS_HEIGHT = 600;

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

/**
 * Detects and resolves collisions between particles and walls in
 * the simulation by clamping positions to stay within bounds and 
 * zeroing velocity along the collision axis.
 */
function handleWallCollisions() {
    for (let i = 0; i < NUMBER_OF_PARTICLES; i++) {
        let p = particles[i];

        // Collision with left wall
        if (p.x - PARTICLE_RADIUS < 0) {
            p.x = PARTICLE_RADIUS;
            p.vx = 0;
        }

        // Collision with right wall
        if (p.x + PARTICLE_RADIUS > CANVAS_WIDTH) {
            p.x = CANVAS_WIDTH - PARTICLE_RADIUS;
            p.vx = 0;
        }

        // Collision with bottom wall
        if (p.y - PARTICLE_RADIUS < 0) {
            p.y = PARTICLE_RADIUS;
            p.vy = 0;
        }

        // Collision with top wall
        if (p.y + PARTICLE_RADIUS > CANVAS_HEIGHT) {
            p.y = CANVAS_HEIGHT - PARTICLE_RADIUS; 
            p.vy = 0;
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
    for (let i = 0; i < NUMBER_OF_PARTICLES; i++) {
        let p = particles[i];

        // Update particle velocity based on gravity
        p.vy += gravity * dt;

        // Update particle position based on velocity
        p.x += p.vx * dt;
        p.y += p.vy * dt;
    }
}

/**
 * Initialises the particle array by placing particles in a non-overlapping grid
 * on the lower-left half of the canvas. Each particle is initialised with zero 
 * velocity and a solid blue color.
 */
function createParticles() {
    const spacing = PARTICLE_DIAMETER * 1.1;
    const particlesPerRow = Math.floor(0.5 * CANVAS_WIDTH / spacing);
    const rows = Math.ceil(NUMBER_OF_PARTICLES / particlesPerRow);

    let count = 0;
    for (let row = 0; row < rows; row++) {
        for (let col = 0; col < particlesPerRow; col++) {
            if (count >= NUMBER_OF_PARTICLES) break;
            const x = col * spacing + PARTICLE_RADIUS;
            const y = row * spacing + PARTICLE_RADIUS;
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

createParticles();
const particlesBufferData = particlesToBuffer();

const uPointSizeLocation = gl.getUniformLocation(program, 'uPointSize');
gl.uniform1f(uPointSizeLocation, PARTICLE_DIAMETER);

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








// dt = 1 / 60;
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


