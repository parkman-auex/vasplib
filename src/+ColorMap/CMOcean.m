classdef CMOcean
    %COLORMAP perceptually uniform colormaps for oceanography
    %Thyng, K.M., C.A. Greene, R.D. Hetland, H.M. Zimmerle, and S.F. DiMarco. 2016. True colors of oceanography: Guidelines for effective and accurate colormap selection. Oceanography 29(3):9–13. http://dx.doi.org/10.5670/oceanog.2016.66
    %   <a href="https://matplotlib.org/cmocean/">See here for names/images of the colormaps</a>
    
    properties (Constant)
        algae = CMOcean.get('algae');
        amp = CMOcean.get('amp');
        balance = CMOcean.get('balance');
        curl = CMOcean.get('curl');
        deep = CMOcean.get('deep');
        delta = CMOcean.get('delta');
        dense = CMOcean.get('dense');
        diff = CMOcean.get('diff');
        gray = CMOcean.get('gray');
        haline = CMOcean.get('haline');
        ice = CMOcean.get('ice');
        matter = CMOcean.get('matter');
        oxy = CMOcean.get('oxy');
        phase = CMOcean.get('phase');
        rain = CMOcean.get('rain');
        solar = CMOcean.get('solar');
        speed = CMOcean.get('speed');
        tarn = CMOcean.get('tarn');
        tempo = CMOcean.get('tempo');
        thermal = CMOcean.get('thermal');
        topo = CMOcean.get('topo');
        turbid = CMOcean.get('turbid');
    end
    
    methods (Static, Access = protected)
        function cmap = get(name)
            folder = fileparts(mfilename('fullpath'));
            C = load(fullfile(folder, 'CMOceanCMOceans.mat'), name);
            cmap = C.(name);
        end
    end
    
end

