window.pychaste_plot = function pychaste_plot(three_container, file_name, width, height) {
    var renderer = new THREE.WebGLRenderer({ antialias: true });
    renderer.setSize(width, height);
    var canvas = renderer.domElement;
    three_container.appendChild(canvas);
    canvas.style.cursor = "move";

    var camera, controls, scene;

    camera = new THREE.PerspectiveCamera(60, window.innerWidth / window.innerHeight, 0.01, 1e10);
    camera.position.z = 6;

    controls = new THREE.OrbitControls(camera, canvas);

    // Only respond to the mouse while the render window is "active" - hovered
    // or clicked - and show a border to indicate this. Deactivate when the
    // mouse leaves or clicks outside, so the plot does not hijack scrolling or
    // dragging while navigating the rest of the notebook.
    function pychaste_set_active(active) {
        controls.enabled = active;
        three_container.style.outline = active ? "1px solid lightgray" : "none";
    }
    pychaste_set_active(false);

    function pychaste_focus_handler(event) {
        if (!document.body.contains(three_container)) {
            document.removeEventListener("mousedown", pychaste_focus_handler);
            return;
        }
        pychaste_set_active(three_container.contains(event.target));
    }
    document.addEventListener("mousedown", pychaste_focus_handler);

    three_container.addEventListener("mouseenter", function () {
        pychaste_set_active(true);
    });
    three_container.addEventListener("mouseleave", function () {
        pychaste_set_active(false);
    });

    scene = new THREE.Scene();
    scene.add(camera);

    // light
    var dirLight = new THREE.DirectionalLight(0xffffff);
    dirLight.position.set(200, 200, 1000).normalize();
    camera.add(dirLight);
    camera.add(dirLight.target);

    var loader = new THREE.VRMLLoader();
    loader.load(
        file_name,
        function (object) {
            scene.add(object);
        },
        undefined,
        function (error) {
            console.error("VRMLLoader error:", error);
        },
    );

    // add a small amount of background ambient light
    var ambientLight = new THREE.AmbientLight(0x444444);
    scene.add(ambientLight);

    function animate() {
        if (camera instanceof THREE.Camera === false || !document.body.contains(three_container)) {
            console.log("Animation loop failed: stopping");
            return;
        }
        requestAnimationFrame(animate);
        controls.update();
        renderer.render(scene, camera);
    }

    animate();
};
