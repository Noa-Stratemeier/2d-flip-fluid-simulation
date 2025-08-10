/**
 * A simple FPS (frames per second) counter that updates a DOM element.
 */
export default class FPSCounter {
    /**
     * Creates an FPSCounter.
     * 
     * @param {string} elementId - The ID of the DOM element where the FPS will be displayed.
     */
    constructor(elementId) {
        this.fpsElement = document.getElementById(elementId);
        this.timeBefore = 0;
        this.frameCount = 0;
    }

    /**
     * Should be called once per rendered frame.
     * Updates the FPS display approximately once per second.
     */
    tick() {
        this.frameCount++;
        let timeNow = performance.now();

        if (timeNow - this.timeBefore >= 1000) {
            this.fpsElement.textContent = this.frameCount;

            this.timeBefore = timeNow;
            this.frameCount = 0;
        }
    }
}