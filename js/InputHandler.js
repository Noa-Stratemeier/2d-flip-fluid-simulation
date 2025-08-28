export default class InputHandler {
    constructor(canvas, scene, onDomainChange) {
        this.canvas = canvas;
        this.scene = scene;
        this.onDomainChange = onDomainChange;

        // Pointer input.
        this.isPointerDown = false;
        this.pointerX = 0;
        this.pointerY = 0;

        canvas.addEventListener("pointerdown", e => this.onDown(e));
        canvas.addEventListener("pointermove", e => this.onMove(e));
        window.addEventListener("pointerup", () => this.onUp());

        // UI open/close.
        this.uiToggleCheckbox = document.getElementById("ui-toggle");
        this.uiToggleLabel = document.getElementById("ui-toggle-label");
        this.uiToggleCheckbox.addEventListener("change", () => this.onUITogglePress());

        // UI input.
        this.bindSlider("pointerInteractionRadius");
        this.bindSlider("pointerInteractionStrength");

        this.bindCheckbox("colourLowDensityParticles");
        this.bindCheckbox("colourParticlesBySpeed");
        this.bindSelect("speedColourMap");
        this.bindSlider("particleDisplaySize");

        this.bindSlider("gravity");
        this.bindSlider("dt");
        this.bindSlider("flipRatio");
        this.bindSlider("projectionIterations");
        this.bindSlider("particleSeparationIterations");
        this.bindSlider("overRelaxation");
        this.bindSlider("stiffness");

        this.bindSlider("relativeFluidWidth", true);
        this.bindSlider("relativeFluidHeight", true);
        this.bindSlider("resolution", true);
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
    // UI open/close.
    // -----------------------------------------------------------------------------

    onUITogglePress() {
        this.uiToggleLabel.textContent = this.uiToggleCheckbox.checked ? "CLOSE" : "OPEN";
    }

    // -----------------------------------------------------------------------------
    // UI input.
    // -----------------------------------------------------------------------------

    bindSlider(sceneKey, triggersDomainRebuild) {
        let slider  = document.getElementById(`${sceneKey}-slider`);
        let valueElement = document.getElementById(`${sceneKey}-value`);

        // Initial UI from scene.
        slider.value = this.scene[sceneKey];
        valueElement.textContent = this.scene[sceneKey];

        // Keep scene + UI in sync.
        slider.addEventListener("input", () => {
            this.scene[sceneKey] = Number(slider.value);
            valueElement.textContent = slider.value;

            if (triggersDomainRebuild) {
                this.onDomainChange();
            }
        });
    }

    bindCheckbox(sceneKey) {
        let checkbox = document.getElementById(`${sceneKey}-toggle`);

        checkbox.checked = !!this.scene[sceneKey];
        checkbox.addEventListener("input", () => {
            this.scene[sceneKey] = checkbox.checked;
        });
    }

    bindSelect(sceneKey) {
        let select = document.getElementById(`${sceneKey}-select`);

        select.value = this.scene[sceneKey];
        select.addEventListener("change", () => {
            this.scene[sceneKey] = select.value;
        });
    }
}