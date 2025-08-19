import FlipFluidSimulation from "./FlipFluidSimulation.js";
import InputHandler from "./InputHandler.js";
import FPSCounter from "./FPSCounter.js";
import * as render from "./render.js";



let canvas = document.getElementById('fluid-canvas');
let gl;
let input = new InputHandler(canvas);
let fpsCounter = new FPSCounter("fps");



let scene = {
    // Colours.
    colourParticlesByDensity: true,
    baseColour: [0.0, 0.0, 1.0],
    lowDensityColour: [1.0, 1.0, 1.0],

    colourParticlesBySpeed: true,

    // Tank.
    tankWidth: window.innerWidth - 10,
    tankHeight: window.innerHeight - 50,
    resolution: 80,

    // Fluid.
    relativeFluidWidth: 0.6,
    relativeFluidHeight: 0.8,

    particleDisplaySize: 1.8,

    // Interaction.
    cursorRepelRadius: 80,
    cursorRepelStrength: 5000,

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
    let cellCountX = Math.floor(scene.tankWidth / cellSize); 
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
    let fluid = scene.flipFluidSimulation = new FlipFluidSimulation(particleCount, particleRadius, cellCountX, cellCountY, cellSize);

    // Match the canvas resolution to the simulation's domain size
    // before creating the WebGL2 context, so rendering coordinates
    // align exactly with the simulation grid.
    canvas.width = fluid.width;
    canvas.height = fluid.height;
    gl = canvas.getContext('webgl2');

    // Initialise particles.
    let i = 0;
    for (let px = 0; px < particleCountX; px++) {
        for (let py = 0; py < particleCountY; py++) {
            let xi = 2 * i;
            let yi = 2 * i + 1;

            let offset = py % 2 === 0 ? 0 : particleRadius;
            fluid.particlePositions[xi] = cellSize + particleRadius + px * particleSpacingX + offset;
            fluid.particlePositions[yi] = cellSize + particleRadius + py * particleSpacingY;

            i++;
        }
    }

    // Initialise particle colours.
    for (let i = 0; i < particleCount; i++) {
        let ri = i * 3;
        let gi = i * 3 + 1;
        let bi = i * 3 + 2;

        fluid.particleColours[ri] = scene.baseColour[0];
        fluid.particleColours[gi] = scene.baseColour[1];
        fluid.particleColours[bi] = scene.baseColour[2];
    }

    // Initialise solid cells.
    for (let gridX = 0; gridX < cellCountX; gridX++) {
        for (let gridY = 0; gridY < cellCountY; gridY++) {
            if (gridX === 0 || gridY === 0 || gridX === cellCountX - 1 || gridY === cellCountY - 1) {
                fluid.solidCells[gridX + gridY * cellCountX] = 0;  // Solid.
            } else {
                fluid.solidCells[gridX + gridY * cellCountX] = 1;  // Not solid.
            }
        }  
    }
}



function animate() {
    fpsCounter.tick();

    scene.flipFluidSimulation.stepSimulation(
        scene.dt, 
        scene.gravity, 
        scene.flipRatio, 
        scene.overRelaxation, 
        scene.particleSeparationIterations, 
        scene.projectionIterations, 
        scene.stiffnessConstant
    );

    // Update particle colours.
    if (scene.colourParticlesByDensity) scene.flipFluidSimulation.updateParticleColoursByDensity(scene.baseColour, scene.lowDensityColour);
    if (scene.colourParticlesBySpeed) scene.flipFluidSimulation.updateParticleColoursBySpeed(800);

    // Apply user interaction.
    if (input.isPointerDown) scene.flipFluidSimulation.repelParticles(input.x, input.y, scene.cursorRepelRadius, scene.cursorRepelStrength);

    render.draw(gl, scene.flipFluidSimulation);

    requestAnimationFrame(animate);
}



initialiseScene();
render.initialise(gl, scene.flipFluidSimulation, scene.particleDisplaySize);
animate();