function [ oi, referenceOi ] = getOiData( varargin )

p = inputParser;
p.KeepUnmatched = true;
p.addOptional('targetDistance',1000);
p.addOptional('depth',5000);
p.addOptional('chlConc',0);
p.addOptional('cdomConc',0);
p.addOptional('smallPartConc',0);
p.addOptional('largePartConc',0);
p.addOptional('target','Macbeth');
p.addOptional('nLeds',7);
p.addOptional('rtbResultFolder',fullfile('/','home','hblasins','Documents','MATLAB','render_toolbox'));

p.parse(varargin{:});
inputs = p.Results;


switch inputs.target
    case 'Macbeth'
        path = fullfile(inputs.rtbResultFolder,'Macbeth','renderings','PBRT');
        scale = 10^14;
    case 'Objects'
        path = fullfile(inputs.rtbResultFolder,'Objects','renderings','PBRT');
        scale = 10^14;
    otherwise
        path = fullfile(inputs.rtbResultFolder,sprintf('Coral-%s',inputs.target),'renderings','PBRT');
        scale = 10^10;
end


%% First get the data we will work with
oi = cell(inputs.nLeds+1,1);  
for i=1:inputs.nLeds+1
    
    switch inputs.target
        case 'Objects'
            fName = fullfile(path,sprintf('Scene_ambient_LED%02i.mat',i-1));
        otherwise
            fName = fullfile(path,sprintf('%s_water_ambient_LED%i_%i_%i_%.3f_%.3f_%.3f_%.3f.mat',inputs.target,i-1,inputs.targetDistance,inputs.depth,inputs.chlConc,inputs.cdomConc,...
                inputs.smallPartConc,inputs.largePartConc));
    end

    try
        radianceData = load(fName);
    catch
        error('Image data not rendered for the specified water parameters.');
    end
        
    % Create an oi
    oiParams.lensType = 'pinhole';
    oiParams.filmDistance = 20;
    oiParams.filmDiag = 20;
    
    
    oi{i} = BuildOI(radianceData.multispectralImage/radianceData.radiometricScaleFactor*scale, [], oiParams);
    if i==1
        oi{i} = oiSet(oi{i},'name',sprintf('Ambient'));
    else
        oi{i} = oiSet(oi{i},'name',sprintf('Ambient + LED%i',i-1));
    end
    
end

%% Now get the reference data to compare to

referenceOi = cell(inputs.nLeds+2,1);

if strcmp(inputs.target,'Objects') == 0
for i=1:inputs.nLeds+2
   
    if i<=inputs.nLeds
        fName = fullfile(path,sprintf('%s_water_noAmbient_LED%i_%i_%i_%.3f_%.3f_%.3f_%.3f.mat',inputs.target,i,inputs.targetDistance,inputs.depth,inputs.chlConc,inputs.cdomConc,...
            inputs.smallPartConc,inputs.largePartConc));
        
        name = sprintf('Water noAmbient LED%i',i);
        
    elseif i==inputs.nLeds+1
        fName = fullfile(path,sprintf('%s_surface_ambient_LED0_%i_%i_%.3f_%.3f_%.3f_%.3f.mat',inputs.target,inputs.targetDistance,inputs.depth,inputs.chlConc,inputs.cdomConc,...
            inputs.smallPartConc,inputs.largePartConc));
        
        name = sprintf('Ambient (no water)');
        
    elseif i==inputs.nLeds+2
        fName = fullfile(path,sprintf('%s_surface_noAmbient_LED8_%i_%i_%.3f_%.3f_%.3f_%.3f.mat',inputs.target,inputs.targetDistance,inputs.depth,inputs.chlConc,inputs.cdomConc,...
            inputs.smallPartConc,inputs.largePartConc)); 
        
        name = sprintf('Ambient flash (no water)');
    end
    
    try
        radianceData = load(fName);
    catch
        error('Image data not rendered for the specified water parameters.');
    end
        
    % Create an oi
    oiParams.lensType = 'pinhole';
    oiParams.filmDistance = 20;
    oiParams.filmDiag = 20;
    
    
    referenceOi{i} = BuildOI(radianceData.multispectralImage/radianceData.radiometricScaleFactor*scale, [], oiParams);
    referenceOi{i} = oiSet(referenceOi{i},'name',name);
    
end
end


end

