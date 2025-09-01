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
        this.uiOpenCloseCheckbox = document.getElementById("ui-open-close");
        this.uiOpenCloseLabel = document.getElementById("ui-open-close-label");
        this.uiElement = document.getElementById("ui");
        this.uiOpenCloseCheckbox.addEventListener("change", () => this.onUIOpenClosePress());

        // UI input.
        this.bindSlider("pointerInteractionRadius");
        this.bindSlider("pointerInteractionStrength");

        this.bindCheckbox("colourLowDensityParticles");
        this.bindCheckbox("colourParticlesBySpeed");
        this.bindSelect("speedColourMap");
        this.bindSlider("particleDisplaySize");
        this.bindColor("baseColour");
        this.bindColor("lowDensityColour");

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

    onUIOpenClosePress() {
        this.uiElement.classList.remove("animation-open", "animation-close");

        if (this.uiOpenCloseCheckbox.checked) {
            this.uiOpenCloseLabel.textContent = "CLOSE";
            this.uiOpenCloseLabel.classList.add("open");
            this.uiElement.classList.add("animation-open");
        } else {
            this.uiOpenCloseLabel.textContent = "OPEN";
            this.uiOpenCloseLabel.classList.remove("open");
            this.uiElement.classList.add("animation-close");
        }
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
        let checkbox = document.getElementById(`${sceneKey}-checkbox`);

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

    bindColor(sceneKey) {
        let colourPicker = document.getElementById(`${sceneKey}-input`);
        let valueElement = document.getElementById(`${sceneKey}-value`);

        // Convert initial RGB array to hex
        colourPicker.value = rgbArrayToHex(this.scene[sceneKey]);
        valueElement.textContent = rgbArrayToHex(this.scene[sceneKey]);

        colourPicker.addEventListener("input", () => {
            this.scene[sceneKey] = hexToRgbArray(colourPicker.value);
            valueElement.textContent = colourPicker.value;
        });
    }  
}

// Utility functions
function rgbArrayToHex(arr) {
    return "#" + arr.map(x => Math.round(x * 255).toString(16).padStart(2, "0")).join("");
}
function hexToRgbArray(hex) {
    hex = hex.replace("#", "");
    return [
        parseInt(hex.substring(0,2), 16) / 255,
        parseInt(hex.substring(2,4), 16) / 255,
        parseInt(hex.substring(4,6), 16) / 255
    ];
}