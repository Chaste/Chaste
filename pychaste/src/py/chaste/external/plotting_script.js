function pychaste_plot(three_container, file_name, width, height) {

    var renderer;
    if (Detector.webgl)
        renderer = new THREE.WebGLRenderer( {antialias:true} );
    else
        renderer = new THREE.CanvasRenderer(); 

    renderer.setSize(width, height);
    var canvas = renderer.domElement;
    three_container.appendChild(canvas);
    canvas.style.cursor = "move";

    var camera, controls, scene;

    camera = new THREE.PerspectiveCamera( 60, window.innerWidth / window.innerHeight, 0.01, 1e10 );
    camera.position.z = 6;

    controls = new THREE.OrbitControls( camera );
    controls.rotateSpeed = 5.0;
    controls.zoomSpeed = 5;
    controls.noZoom = false;
    controls.noPan = false;

    scene = new THREE.Scene();
    scene.add( camera );

    // light
    var dirLight = new THREE.DirectionalLight( 0xffffff );
    dirLight.position.set( 200, 200, 1000 ).normalize();
    camera.add( dirLight );
    camera.add( dirLight.target);

    var loader = new THREE.VRMLLoader();
    loader.load(file_name,  function(object){ scene.add(object)});

    // add a small amount of background ambient light
    var ambientLight = new THREE.AmbientLight(0x444444);
    scene.add(ambientLight);    

    function animate() {
        if ( camera instanceof THREE.Camera === false || ! document.body.contains(three_container)) {
            console.log("Animation loop failed: stopping");
            return;
        }
        requestAnimationFrame( animate );
        controls.update();
        renderer.render( scene, camera );
    }

    animate();
}
