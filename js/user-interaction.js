export default function initialiseUserInteraction(canvas, scene) {
    
    // Mouse interaction.
    canvas.addEventListener('mousedown', handleMouseInteraction);
    canvas.addEventListener('mousemove', handleMouseInteraction);

    // Touch interaction (mobile).
    canvas.addEventListener('touchstart', handleTouchInteraction, { passive: false });
    canvas.addEventListener('touchmove', handleTouchInteraction, { passive: false });

    // -----------------------------------------------------------------------------
    // Handle user interaction functions.
    // -----------------------------------------------------------------------------

    function handleMouseInteraction(e) {
        if (e.buttons !== 1) return;  // Left mouse button.
        
        // Get mouse position in simulation coordinates.
        let { x, y } = getSimulationCoordinates(e.clientX, e.clientY);
        
        // Push particles away from mouse.
        scene.flipFluidSimulation.repelParticles(x, y, scene.cursorRepelRadius, scene.cursorRepelStrength);
    }

    function handleTouchInteraction(e) {
        e.preventDefault();  // Prevent scrolling.

        if (e.touches.length === 0) return;

        let touch = e.touches[0];
        let { x, y } = getSimulationCoordinates(touch.clientX, touch.clientY);

        scene.flipFluidSimulation.repelParticles(x, y, 50, scene.cursorRepelStrength);
    }

    // -----------------------------------------------------------------------------
    // Helper functions.
    // -----------------------------------------------------------------------------

    function getSimulationCoordinates(clientX, clientY) {
        let rect = canvas.getBoundingClientRect();
        return {
            x: clientX - rect.left,
            y: rect.height - (clientY - rect.top)  // Flip Y axis.
        };
    }
}
