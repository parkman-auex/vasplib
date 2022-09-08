%% An example for the visualization of 3D-bandstructure for your own by using MATLAB
% usage: figure_k = bandplot_3d(bandlist,Ecut,contourf_mode,wt_mode)
function varargout=bandplot_3d(varargin)
    [varargout{1:nargout}] = vasplib_plot.bandplot_3d(varargin{:});
end