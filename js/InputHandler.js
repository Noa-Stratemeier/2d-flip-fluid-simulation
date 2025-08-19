export default class InputHandler {
    constructor(canvas) {
        this.canvas = canvas;
        this.isPointerDown = false;
        this.x = 0;
        this.y = 0;

        canvas.addEventListener("pointerdown", e => this.onDown(e));
        canvas.addEventListener("pointermove", e => this.onMove(e));
        window.addEventListener("pointerup", () => this.onUp());
    }

    onDown(e) {
        this.isPointerDown = true;
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
        this.x = x;
        this.y = y;
    }

    getSimulationCoordinates(clientX, clientY) {
        let rect = this.canvas.getBoundingClientRect();
        return {
            x: clientX - rect.left,
            y: rect.height - (clientY - rect.top)  // Flip Y axis.
        };
    }
}