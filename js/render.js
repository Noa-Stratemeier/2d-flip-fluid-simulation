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

export function initialise(gl, simulation, pointScale) {
    let program = createShaderProgram(gl, vertexShaderSource, fragmentShaderSource);
    gl.useProgram(program);

    // Get and set uniform locations.
    let uPointSizeLocation = gl.getUniformLocation(program, 'uPointSize');
    let uSimulationDomainSizeLocation = gl.getUniformLocation(program, 'uSimulationDomainSize');

    gl.uniform1f(uPointSizeLocation, 2.0 * simulation.particleRadius * pointScale);
    gl.uniform2f(uSimulationDomainSizeLocation, simulation.width, simulation.height);

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

export function draw(gl, simulation) {
    // Update positions buffer.
    gl.bindBuffer(gl.ARRAY_BUFFER, particlePositionsBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, simulation.particlePositions, gl.DYNAMIC_DRAW);
    
    // Update colours buffer.
    gl.bindBuffer(gl.ARRAY_BUFFER, particleColoursBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, simulation.particleColours, gl.DYNAMIC_DRAW);

    // Clear and draw.
    gl.clearColor(0.0, 0.0, 0.0, 1.0);
    gl.clear(gl.COLOR_BUFFER_BIT);
    gl.drawArrays(gl.POINTS, 0, simulation.particleCount);
}