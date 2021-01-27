%% femm_plot_manual.m

% Code intended to plot solution from a *.ans file, with the results
% calculated from FEMM.

close all
clear all

file_name = 'reference';
% file_name = 'mag_tutorial';

% Copy just to test updated extent_mfemm() function:
% file_name = 'reference - Copia';

FemmProblem = loadfemmfile([file_name '.fem']);
ansfile = [file_name '.ans'];

%%
%
% Finally we can take a look at the problem using Matlab's plotting
% commands

% plotfemmproblem(FemmProblem);

%% Extracting Results
%
% Having solved te problem, it is generally desired that we do some
% post-processing on the results. The file which is produced is perfectly
% compatible with FEMM, and if it is installed, we can open the file and
% use all the normal post-processing functions
% if exist('openfemm.m', 'file')
%     opendocument(ansfile);
% else
%     fprintf(1, 'Looks like femm isn''t installed, or at least its m-files aren''t on the path.\n');
% end
%
% However, some post-processing methods are also provided with mfemm. The
% main method is a class based method for loading and manipulating the
% data. The class is called fpproc, and also requires some compilation
% before use as it is an interface to more C++ code.
%
% So to manipulate the output, we first create an fpproc class
myfpproc = fpproc();

%%
%
% and then load the results file using the opendocument method
myfpproc.opendocument(ansfile);

%%

% The original code doesn't recognize arc segments when establishing the
% problem extents. This was updated.
% [x,y,w,h] = extent_mfemm(FemmProblem);


%%

% Plot using Cylindrical Coordinates
myfpproc.plotBfield(0, 0, 10, 90, 'Method', 1, 'PlotCylindrical', true);

% The vector potential can also be plotted as density plot
myfpproc.plotAfield();
