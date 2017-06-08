function [ scene, mappings ] = MakeCoralMexximpRemodeller( scene, mappings, names, conditionValues, conditionNumber )

distance = rtbGetNamedNumericValue(names,conditionValues,'cameraDistance',[]);
lookAt = eval(rtbGetNamedValue(names,conditionValues,'lookAt',[]));
viewAxis = eval(rtbGetNamedValue(names,conditionValues,'viewAxis',[]));


hasAmbient = rtbGetNamedValue(names,conditionValues,'illumination',[]);



% Point the camera towards the scene
[scene, camera, cameraNode] = mexximpCentralizeCamera(scene,'viewAxis',viewAxis,...
                                                            'viewUp',[0.1 0.1 1]);


% And move to the specified distance
cameraPosition = mexximpApplyTransform(camera.position,cameraNode.transformation);
currentDistance = sqrt(sum(cameraPosition.^2));
ratio = (distance/currentDistance);

newPosition = cameraPosition*ratio;

cameraSelector = strcmp('Camera',{scene.rootNode.children.name});
scene.rootNode.children(cameraSelector).transformation = mexximpLookAt(newPosition,lookAt,[0 0 1]');

%% Point source

light = mexximpConstants('light');
light.position = [0 0 0]';
light.type = 'point';
light.name = 'Flash';
light.lookAtDirection = [0 0 0]';
light.ambientColor = 100*[1 1 1]';
light.diffuseColor = 100*[1 1 1]';
light.specularcolor = 100*[1 1 1]';
light.constantAttenuation = 1;
light.linearAttenuation = 0;
light.quadraticAttenuation = 1;
light.innerConeAngle = 0;
light.outerConeAngle = 0;

scene.lights = light;

lightNode = mexximpConstants('node');
lightNode.name = light.name;
lightNode.transformation = mexximpTranslate(newPosition + [50; 50; 50]);

scene.rootNode.children = [scene.rootNode.children, lightNode];


%% Ambient illumination
% If we decide that there is ambient illumination, we need to add it.
if strcmp(hasAmbient,'ambient')


    ambient = mexximpConstants('light');
    ambient.position = [0 0 0]';
    ambient.type = 'directional';
    ambient.name = 'SunLight';
    ambient.lookAtDirection = [0 0 -1]';
    ambient.ambientColor = 10*[1 1 1]';
    ambient.diffuseColor = 10*[1 1 1]';
    ambient.specularcolor = 100*[1 1 1]';
    ambient.constantAttenuation = 1;
    ambient.linearAttenuation = 0;
    ambient.quadraticAttenuation = 1;
    ambient.innerConeAngle = 0;
    ambient.outerConeAngle = 0;
    
    scene.lights = [scene.lights, ambient];
    
    ambientNode = mexximpConstants('node');
    ambientNode.name = ambient.name;
    ambientNode.transformation = eye(4);
    
    scene.rootNode.children = [scene.rootNode.children, ambientNode];

    
    ambient = mexximpConstants('light');
    ambient.position = [0 0 0]';
    ambient.type = 'directional';
    ambient.name = 'SunLight2';
    ambient.lookAtDirection = [1 1 -0.5]';
    ambient.ambientColor = 10*[1 1 1]';
    ambient.diffuseColor = 10*[1 1 1]';
    ambient.specularcolor = 100*[1 1 1]';
    ambient.constantAttenuation = 1;
    ambient.linearAttenuation = 0;
    ambient.quadraticAttenuation = 1;
    ambient.innerConeAngle = 0;
    ambient.outerConeAngle = 0;
    
    scene.lights = [scene.lights, ambient];
    
    ambientNode = mexximpConstants('node');
    ambientNode.name = ambient.name;
    ambientNode.transformation = eye(4);
    
    scene.rootNode.children = [scene.rootNode.children, ambientNode];
    
    ambient = mexximpConstants('light');
    ambient.position = [0 0 0]';
    ambient.type = 'directional';
    ambient.name = 'SunLight3';
    ambient.lookAtDirection = [1 -1 -0.5]';
    ambient.ambientColor = 10*[1 1 1]';
    ambient.diffuseColor = 10*[1 1 1]';
    ambient.specularcolor = 100*[1 1 1]';
    ambient.constantAttenuation = 1;
    ambient.linearAttenuation = 0;
    ambient.quadraticAttenuation = 1;
    ambient.innerConeAngle = 0;
    ambient.outerConeAngle = 0;
    
    scene.lights = [scene.lights, ambient];
    
    ambientNode = mexximpConstants('node');
    ambientNode.name = ambient.name;
    ambientNode.transformation = eye(4);
    
    scene.rootNode.children = [scene.rootNode.children, ambientNode];
    
    ambient = mexximpConstants('light');
    ambient.position = [0 0 0]';
    ambient.type = 'directional';
    ambient.name = 'SunLight4';
    ambient.lookAtDirection = [-1 1 -0.5]';
    ambient.ambientColor = 10*[1 1 1]';
    ambient.diffuseColor = 10*[1 1 1]';
    ambient.specularcolor = 100*[1 1 1]';
    ambient.constantAttenuation = 1;
    ambient.linearAttenuation = 0;
    ambient.quadraticAttenuation = 1;
    ambient.innerConeAngle = 0;
    ambient.outerConeAngle = 0;
    
    scene.lights = [scene.lights, ambient];
    
    ambientNode = mexximpConstants('node');
    ambientNode.name = ambient.name;
    ambientNode.transformation = eye(4);
    
    scene.rootNode.children = [scene.rootNode.children, ambientNode];
    
    ambient = mexximpConstants('light');
    ambient.position = [0 0 0]';
    ambient.type = 'directional';
    ambient.name = 'SunLight5';
    ambient.lookAtDirection = [-1 -1 -0.5]';
    ambient.ambientColor = 10*[1 1 1]';
    ambient.diffuseColor = 10*[1 1 1]';
    ambient.specularcolor = 100*[1 1 1]';
    ambient.constantAttenuation = 1;
    ambient.linearAttenuation = 0;
    ambient.quadraticAttenuation = 1;
    ambient.innerConeAngle = 0;
    ambient.outerConeAngle = 0;
    
    scene.lights = [scene.lights, ambient];
    
    ambientNode = mexximpConstants('node');
    ambientNode.name = ambient.name;
    ambientNode.transformation = eye(4);
    
    scene.rootNode.children = [scene.rootNode.children, ambientNode];
end

end

