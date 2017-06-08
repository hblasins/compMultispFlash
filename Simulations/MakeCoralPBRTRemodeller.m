function [ nativeScene ] = MakeCoralPBRTRemodeller( parentScene, nativeScene, mappings, names, conditionValues, conditionNumbers )

mode = rtbGetNamedValue(names,conditionValues,'mode',[]);
ledID = rtbGetNamedNumericValue(names,conditionValues,'ledID',[]);

filmDiag = rtbGetNamedNumericValue(names,conditionValues,'filmDiag',[]);
filmDist = rtbGetNamedNumericValue(names,conditionValues,'filmDist',[]);
pixelSamples =  rtbGetNamedNumericValue(names,conditionValues,'pixelSamples',[]);
volumeStep = rtbGetNamedNumericValue(names,conditionValues,'volumeStep',[]);
depth = rtbGetNamedNumericValue(names,conditionValues,'depth',[]);

chlConc = rtbGetNamedNumericValue(names,conditionValues,'chlConc',[]);
cdomConc = rtbGetNamedNumericValue(names,conditionValues,'cdomConc',[]);
smallPartConc = rtbGetNamedNumericValue(names,conditionValues,'smallPartConc',[]);
largePartConc = rtbGetNamedNumericValue(names,conditionValues,'largePartConc',[]);

hasAmbient = rtbGetNamedValue(names,conditionValues,'illumination',[]);



camera = nativeScene.overall.find('Camera');
camera.type = 'pinhole';
camera.setParameter('filmdiag','float',filmDiag);
camera.setParameter('filmdistance','float',filmDist);
        
sampler = nativeScene.overall.find('Sampler');
sampler.setParameter('pixelsamples','integer',pixelSamples);

integrator = nativeScene.overall.find('SurfaceIntegrator');
integrator.type = 'path';

if strcmp(hasAmbient,'ambient')
    distantLight = nativeScene.world.find('LightSource','name','SunLight');
    distantLight.setParameter('L','spectrum','resources/SunLight.spd');
    
    distantLight = nativeScene.world.find('LightSource','name','SunLight2');
    distantLight.setParameter('L','spectrum','resources/SunLight.spd');
    
    distantLight = nativeScene.world.find('LightSource','name','SunLight3');
    distantLight.setParameter('L','spectrum','resources/SunLight.spd');
    
    distantLight = nativeScene.world.find('LightSource','name','SunLight4');
    distantLight.setParameter('L','spectrum','resources/SunLight.spd');
    
    distantLight = nativeScene.world.find('LightSource','name','SunLight5');
    distantLight.setParameter('L','spectrum','resources/SunLight.spd');
end

pointLight = nativeScene.world.find('LightSource','name','Flash');
pointLight.setParameter('I','spectrum',sprintf('resources/LED%i.spd',ledID));


switch mode
    case 'water'
        nativeScene.overall.find('SurfaceIntegrator','remove',true);
        volumeIntegrator = MPbrtElement('VolumeIntegrator','type','single');
        volumeIntegrator.setParameter('stepsize','float',volumeStep);
        nativeScene.overall.append(volumeIntegrator);
        
        
        waterVolume = MPbrtElement('Volume','type','water');
        waterVolume.setParameter('p0','point','-5000 -5000 -1000');
        waterVolume.setParameter('p1','point',sprintf('5000 5000 %i',depth));
        waterVolume.setParameter('absorptionCurveFile','spectrum',sprintf('resources/abs_%.3f_%.3f.spd',chlConc,cdomConc));
        waterVolume.setParameter('phaseFunctionFile','string',sprintf('resources/phase_%.3f_%.3f.spd',smallPartConc,largePartConc));
        waterVolume.setParameter('scatteringCurveFile','spectrum',sprintf('resources/scat_%.3f_%.3f.spd',smallPartConc,largePartConc));
        nativeScene.world.append(waterVolume);
        
    case 'depth'
                
        sampler = nativeScene.overall.find('Sampler');        
        sampler.type = 'stratified';
        sampler.setParametr('xsamples','integer',1);
        sampler.setParamter('ysamples','integer',1);
        sampler.setParamter('jitter','bool','false');

                
        nativeScene.overall.find('PixelFilter');
        filter = MPbrtElement('PixelFilter','type','box');
        filter.setParameter('alpha','float',2);
        filter.setParameter('xwidth','float',0.5);
        filter.setParameter('ywidth','float',0.5);
        
    otherwise
        
        filter = nativeScene.overall.find('PixelFilter');
        filter.type = 'box';
        filter.setParameter('alpha','float',2);
        filter.setParameter('xwidth','float',0.5);
        filter.setParameter('ywidth','float',0.5);

end
        
        

      

end


