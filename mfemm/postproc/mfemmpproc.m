classdef mfemmpproc < handle
    % mfemmpproc - base class for xfemm post-processing C++ interfaces


    properties (SetAccess = protected, Hidden = true)
        
        objectHandle; % Handle to the underlying C++ class instance
        
    end
    
    properties (SetAccess = protected, Hidden = false)

        isdocopen = false; % flag denoting whether a document has been opened yet
        
        openfilename = ''; % contains the path of any currently open files
        
        FemmProblem = struct (); % The femmproblem structure representing the problem
        
    end
    
    methods (Access = protected)
          
          function hfig = plotvectorfield(this, datafcn, x, y, w, h, varargin)
              % creates a plot of a vector field on the FemmProblem
              %
              % Syntax
              %
              % hpproc.plotvectorfield(method, x, y, w, h, points, datafcn)
              %
              % Input
              %
              %   method - plot method use 0 for a vector field plot using
              %     coloured arrows. Use 1 for a contour plot of the
              %     magnitude of the magnetic field.
              %
              %   x - x (or r) coordinate lower left corner of region to be
              %     plotted
              % 
              %   y - y (or x) coordinate of  lower left corner of region to 
              %     be plotted
              %
              %   w - width of region to be plotted
              % 
              %   h - height of region to be plotted
              % 
              % Further options are supplied using parameter-value pairs,
              % the possible options are:
              %
              %   'Points' - determines the number of points that will be
              %     plotted using method 0. If points is a scalar, a grid of
              %     this number of points on both sides will be created. If
              %     points is a two element vector it will be the number of
              %     points in the x and y direction respectively.
              % 
              %    'PlotNodes' - determines if nodes are drawn when
              %      plotting the femmproblem
              %
              %   'Method' - plot method use 0 for a vector field plot using
              %     coloured arrows. Use 1 for a contour plot of the
              %     magnitude of the magnetic field.
              %
            
              Inputs.Points = 40;
              Inputs.Density = 250;
              Inputs.FluxLines = 19;
              Inputs.ShowFlux = true;
              Inputs.Method = 0;
              Inputs.PlotCylindrical = false;
              Inputs.PlotNodes = false;
              Inputs.AddLabels = false;
              
              Inputs = mfemmdeps.parseoptions (Inputs, varargin);
              
              % Romoving generic problem plot. Using simple plot after the
              % results visualization
              
%               if ~isempty (this.FemmProblem)
%                   [hfig, hax] = plotfemmproblem(this.FemmProblem, ...
%                                     'PlotNodes', Inputs.PlotNodes, ...
%                                     'InitialViewPort', [x,y,w,h], ...
%                                     'AddLabels', Inputs.AddLabels);
%               else
%                   hfig = figure;
%               end

              hfig = figure;
              
              if isscalar(Inputs.Points)
                  Inputs.Points = [Inputs.Points, Inputs.Points];
              end
              
              switch Inputs.Method
                  
                  case 0
                      
                      if Inputs.PlotCylindrical
                            % Since the vector plot works better with
                            % cartesian coordinates, lets translated
                            % This here is to adjust the plot extent with
                            % cylindrical inputs.

                            d_ang = h*pi/180;
                            ang = y*pi/180;
                            r = x + w;
                            
                            x = r*cos(ang + d_ang);

                            if x > 0
                                x = 0;
                            end

                            w = r*cos(ang) - x;

                            y = 0;
                            h = r*sin(ang + d_ang);

                            if d_ang > pi/2
                                h = r*sin(pi/2);
                            end
                          
                          xgv = linspace(x, x + w, Inputs.Points(1));
                          ygv = linspace(y, y + h, Inputs.Points(2));
                          [Xsample,Ysample] = meshgrid(xgv, ygv);
                          
                      else

                          xgv = linspace(x, x + w, Inputs.Points(1));
                          ygv = linspace(y, y + h, Inputs.Points(2));
                          [Xsample,Ysample] = meshgrid(xgv, ygv);
                          
                      end

                      data = feval(datafcn, Xsample, Ysample);
              
                      % plot a vector field using colored arrows
                      hax = mfemmdeps.cquiver( cat(3, reshape(data(1,:), size(Xsample)), reshape(data(2,:), size(Xsample))), ...
                                               'sx', xgv(2)-xgv(1), ...
                                               'sy', ygv(2)-ygv(1), ...
                                               'xshift', x, ...
                                               'yshift', y); %, ...
%                                'hax', hax );
                           
                       colorbar;
                       
                  case 1
                      
                      if Inputs.PlotCylindrical
                          
                          rgv = linspace(x, x + w, Inputs.Density);
                          thgv = linspace(y, y + h, Inputs.Density);
                          Xsample = rgv.'*cos(thgv*pi/180);
                          Ysample = rgv.'*sin(thgv*pi/180);
                          
                      else

                          xgv = linspace(x, x + w, Inputs.Density);
                          ygv = linspace(y, y + h, Inputs.Density);
                          [Xsample,Ysample] = meshgrid(xgv, ygv);
                          
                      end
                      
                      data = feval(datafcn, Xsample, Ysample);
                      
                      % Density plot                      
                      
                      u = reshape(data(1,:), size(Xsample));
                      v = reshape(data(2,:), size(Xsample));
%                       maxd = sqrt(max(u(:).^2 + v(:).^2));
                      % Remove outliers of the data;
                      maxd = prctile(sqrt(u(:).^2 + v(:).^2),99.9);
                      

                      pcolor(Xsample,Ysample,reshape(mfemmdeps.magn(data), size(Xsample)));
                      hax = hfig.CurrentAxes;

                      shading interp;
                      set(hax, 'CLim', [0 maxd]);
                      colorbar;
                      colormap jet;
                      
                      
%                       map = jet(50);
%                       colormap(map(20:50,:));
%                       map = jet(35);
%                       colormap(map(15:35,:))

%                       colormap parula;
                      
                      
                      
                      % Original code for contour
                      %{
%                       cont_lines = max(Inputs.Points);
                      cont_lines = 19;
                      
                      hold on
                      contour( Xsample, Ysample, reshape(mfemmdeps.magn(data), size(Xsample)), cont_lines );
                      set(hax, 'CLim', [0 maxd]);
                      colorbar;
                      %}
                      
                  otherwise
                          
              end
              
              hold off;
              
               %{
              % Going to leave this out of here.
              axis equal
              
              if Inputs.PlotCylindrical
                  set (hax, 'XLim', [0, x + w], 'YLim', [0, x + w]);
              else
                  set (hax, 'XLim', [x, x + w], 'YLim', [y, y + h]);
              end
              %}
              
          end
          
          function hfig = plotscalarfield(this, datafcn, x, y, w, h, varargin)
              % creates a plot of a scalar field on the FemmProblem
              %
              % Syntax
              %
              % hpproc.plotscalarfield(method, x, y, w, h, points, datafcn)
              %
              % Input
              %
              %
              %   x - x (or r) coordinate lower left corner of region to be
              %     plotted
              % 
              %   y - y (or x) coordinate of  lower left corner of region to 
              %     be plotted
              %
              %   w - width of region to be plotted
              % 
              %   h - height of region to be plotted
              %
              % Further options are supplied using parameter-value pairs,
              % the possible options are:
              %
              %   'Points' - determines the number of points that will be
              %     plotted using method 0. If points is a scalar, a grid of
              %     this number of points on both sides will be created. If
              %     points is a two element vector it will be the number of
              %     points in the x and y direction respectively.
              %
              %    'Method' - plot method use 0 for a filled contour plot.
              %      Method 1 is used only to add.
              %
              %    'FigureHandle' - using the method to add to another
              %      plot.
              %
            
              Inputs.Points = 40;
              Inputs.Density = 250;
              Inputs.FluxLines = 19;
              Inputs.ShowFlux = true;
              Inputs.Method = 0;
              Inputs.PlotCylindrical = false;
              Inputs.PlotNodes = false;
              Inputs.AddLabels = false;
              Inputs.FigureHandle = [];
              Inputs.AxesHandle = [];
              
              Inputs = mfemmdeps.parseoptions (Inputs, varargin);
              
              if isscalar(Inputs.Points)
                  Inputs.Points = [Inputs.Points, Inputs.Points];
              end

              % Add this code so that the function can be called to add
              % flux contours to a plot.
              if isempty (Inputs.FigureHandle)
                  hfig = figure;
                  hax = hfig.CurrentAxes;
              else
                  hfig = Inputs.FigureHandle;
                  hax = hfig.CurrentAxes;
              end
              
              if Inputs.PlotCylindrical

                  rgv = linspace(x, x + w, Inputs.Density);
                  thgv = linspace(y, y + h, Inputs.Density);
                  Xsample = rgv.'*cos(thgv*pi/180);
                  Ysample = rgv.'*sin(thgv*pi/180);
                  
              else

                  xgv = linspace(x, x + w, Inputs.Density);
                  ygv = linspace(y, y + h, Inputs.Density);
                  [Xsample,Ysample] = meshgrid(xgv, ygv);
                  
              end

              data = feval(datafcn, Xsample, Ysample);
              
              switch Inputs.Method
                  
                  case 0
                      
                      pcolor(Xsample,Ysample,data);

                      shading interp;
                      
                      if Inputs.ShowFlux
                      
                          hold on

                          contour ( Xsample, Ysample, data, Inputs.FluxLines, 'Color', 'white');

                          hold off
                      
                      end
                      
                      colorbar;
                      colormap jet;
                       
                  case 1
                      % contour plot
                        
                      figure(hfig)

                      hold on
                        
                      contour ( Xsample, Ysample, data, Inputs.FluxLines, 'Color', 'white', 'LineWidth', 0.1);
                                            
                      hold off
                      
                  otherwise
                          
              end
              
              %{
              % Going to leave this out of here.
              axis equal
              
              if Inputs.PlotCylindrical
                  set (hax, 'XLim', [0, x + w], 'YLim', [0, x + w]);
              else
                  set (hax, 'XLim', [x, x + w], 'YLim', [y, y + h]);
              end
              %}
              
          end
          
          
      end

end