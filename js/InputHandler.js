export default class InputHandler {
    constructor(canvas, scene) {
        this.canvas = canvas;
        this.scene = scene;

        // Pointer input.
        this.isPointerDown = false;
        this.pointerX = 0;
        this.pointerY = 0;

        canvas.addEventListener("pointerdown", e => this.onDown(e));
        canvas.addEventListener("pointermove", e => this.onMove(e));
        window.addEventListener("pointerup", () => this.onUp());

        // UI input.
        this.bindSlider("pointerInteractionRadius", v => v);
        this.bindSlider("pointerInteractionStrength", v => v);

        this.bindCheckbox("colourLowDensityParticles");
        this.bindCheckbox("colourParticlesBySpeed");
        this.bindSlider("particleDisplaySize", v => v);

        this.bindSlider("gravity", v => v);
        this.bindSlider("dt", v => v);
        this.bindSlider("flipRatio", v => v);
        this.bindSlider("projectionIterations", v => v);
        this.bindSlider("particleSeparationIterations", v => v);
        this.bindSlider("overRelaxation", v => v);
        this.bindSlider("stiffness", v => v);
    }

    // -----------------------------------------------------------------------------
    // Pointer input.
    // -----------------------------------------------------------------------------

    onDown(e) {
        if (e.button === 0) {
            this.isPointerDown = true;
        }
        this.updatePointerPosition(e);
    }

    onMove(e) {
        if (this.isPointerDown) this.updatePointerPosition(e);
    }

    onUp() {
        this.isPointerDown = false;
    }

    updatePointerPosition(e) {
        let {x, y} = this.getSimulationCoordinates(e.clientX, e.clientY);
        this.pointerX = x;
        this.pointerY = y;
    }

    getSimulationCoordinates(clientX, clientY) {
        let rect = this.canvas.getBoundingClientRect();
        return {
            x: clientX - rect.left,
            y: rect.height - (clientY - rect.top)  // Flip Y axis.
        };
    }

    // -----------------------------------------------------------------------------
    // UI input.
    // -----------------------------------------------------------------------------

    bindSlider(sceneKey, format = v => v) {
        let slider  = document.getElementById(`${sceneKey}-slider`);
        let valueElement = document.getElementById(`${sceneKey}-value`);

        // initial UI from scene.
        slider.value = this.scene[sceneKey];
        valueElement.textContent = format(this.scene[sceneKey]);

        // Keep scene + UI in sync.
        slider.addEventListener("input", () => {
            this.scene[sceneKey] = Number(slider.value);
            valueElement.textContent = format(slider.value);
        });
    }

    bindCheckbox(sceneKey) {
        let checkbox = document.getElementById(`${sceneKey}-toggle`);
        if (!checkbox) return;
        checkbox.checked = !!this.scene[sceneKey];
        checkbox.addEventListener("input", () => {
            this.scene[sceneKey] = checkbox.checked;
        });
    }
}