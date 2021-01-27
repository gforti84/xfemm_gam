function [x,y,w,h] = extent_mfemm(FemmProblem)
% extent_mfemm: gets the extent of a problem described in an mfemm
% structure
%
% Syntax
%
% [x,y,w,h] = extent_mfemm(FemmProblem)
%
% Input
%
% FemmProblem - an mfemm problem structure
%
% Output
%
% If there are no nodes in the problem all outputs will be a single Nan
% value, otherwise:
%
% x - x coordinate of lower left corner of simulation region
% 
% y - y coordinate of lower left corner of simulation region
% 
% w - width of problem region
% 
% h - height of problem region
%

% Copyright 2012 Richard Crozier
% 
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
% 
%        http://www.apache.org/licenses/LICENSE-2.0
% 
%    Unless required by applicable law or agreed to in writing, software
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

    nodecoords = getnodecoords_mfemm(FemmProblem);
    
    num_arc = length(FemmProblem.ArcSegments);
    
    if ~isempty (nodecoords)
        
        if num_arc > 0
            
            %
            
            % Find the extreme node positions, not considering the arc
            % segments
            x_max = max(nodecoords(:,1));
            y_max = max(nodecoords(:,2));
            x_min = min(nodecoords(:,1));
            y_min = min(nodecoords(:,2));
            
%             xy_max = max(x_max, y_max);
            
            xy_set = [1*x_max, 1*y_max; ...
                      1*x_max, 0*y_max; ...
                      1*x_max, 1*y_min; ...
                      0*x_max, 1*y_max; ...
                      0*x_max, 1*y_min; ...
                      1*x_min, 1*y_max; ...
                      1*x_min, 0*y_max; ...
                      1*x_min, 1*y_min];
            
            result = mfemmdeps.ipdm(xy_set, nodecoords, 'Result', 'Structure', 'Subset', 'NearestNeighbor');
            % Results in two columns containing the node id and distance
            % from teste (x,y) position.
            
            % Arrange the arc elements in a 2xN matrix
            arc_elem = [cell2mat({FemmProblem.ArcSegments.n0}); ... 
                        cell2mat({FemmProblem.ArcSegments.n1})];
                    
            % arc_ang is the angle span (in radians)
            arc_ang = (pi/180)*cell2mat({FemmProblem.ArcSegments.ArcLength});
            
            id_node = [];
            arc_pos = [];
            
            for idx_res = 1:length(result.columnindex)
                
                res_id = result.columnindex(idx_res) - 1;
                res_x_max_dis = abs(FemmProblem.Nodes(res_id + 1).Coords(1) - x_max);
                res_y_max_dis = abs(FemmProblem.Nodes(res_id + 1).Coords(2) - y_max);
                res_x_min_dis = abs(FemmProblem.Nodes(res_id + 1).Coords(1) - x_min);
                res_y_min_dis = abs(FemmProblem.Nodes(res_id + 1).Coords(2) - y_min);
                
                % Testa se os nós encontrados estão em algum extremo
                if res_x_max_dis == 0 || res_y_max_dis == 0 || ... 
                   res_x_min_dis == 0 || res_y_min_dis == 0
                    
                    id_node = [id_node; res_id];
                    % Arc segments conected to these elements
                    arc_pos = [arc_pos; ceil(find(arc_elem == res_id)/2)];
                    
                end
                
            end
            
            id_node = unique(id_node);
            % Remove repeated results (for example, two segments conected to same nodes, forming a circle)
            count_arc = length(unique(arc_pos,'rows','first'));

            for idx_node = 1:length(id_node)
                % It runs all nodes listed in the current extremes.
                % Now, this checks if there is any arc segment extending
                % the extremes.
                
                % arc_elem_pos is the position of id_node in the arc
                % segment, n0 or n1. It sets the direction of the arc (countclockwise)
                arc_elem_pos = find(arc_elem == id_node(idx_node));

                for idx = 1:length(arc_elem_pos)
                    
                    if mod(arc_elem_pos(idx),2) == 1
                        
                        % odd case, node is in the first row.(n0)
                        arc_xy1 = FemmProblem.Nodes(arc_elem(arc_elem_pos(idx)) + 1).Coords;
                        arc_xy2 = FemmProblem.Nodes(arc_elem(arc_elem_pos(idx) + 1) + 1).Coords;

                        id_node = id_node(~(id_node == arc_elem(arc_elem_pos(idx) + 1)));
                        
                    elseif mod(arc_elem_pos(idx),2) == 0
                        
%                         % even case, node is in the second row.
                        arc_xy1 = FemmProblem.Nodes(arc_elem(arc_elem_pos(idx) - 1) + 1).Coords;
                        arc_xy2 = FemmProblem.Nodes(arc_elem(arc_elem_pos(idx)) + 1).Coords;
                        
                        id_node = id_node(~(id_node == arc_elem(arc_elem_pos(idx) - 1)));
                        
                    end
                        
                        dist =  sqrt((arc_xy1(1) - arc_xy2(1))^2 + (arc_xy1(2) - arc_xy2(2))^2);
                        radius = (dist/2) / (sin(arc_ang(ceil(arc_elem_pos(idx)/2))/2));
                        
                        % Distance from circle center to dist segment
                        dis_c = radius*cos(arc_ang(ceil(arc_elem_pos(idx)/2))/2);
                        
                        % Finding the 'Sagitta' - hight of the arc segment
                        % from the segment xy1 to xy2
                        hight = radius - dis_c;
                        
                        del_xy = (arc_xy2 - arc_xy1);
                        mean_xy = (arc_xy2 + arc_xy1)/2;
                        
                        dir = atan2(del_xy(2),del_xy(1))*180/pi;
                        dir = mod(dir, 360);  % Assure positive, from 0 to 360°
                                                
                        % For now, not calculating accurately. Lets gor for
                        % the wurst käse scenario. 
                        
                        % There are four quadrants:
%                              -45 to  45 -> increase y_min;
%                               45 to 135 -> increase x_max;
%                              135 to 225 -> increase y_max;
%                              225 to -45 -> increase x_min;

                        quad_dir = ceil((dir + 45)/90);
                        
                        switch quad_dir
                            case 1
                                temp_y_min = mean_xy(2) - hight;
                        
                                if temp_y_min < y_min
                                    y_min = temp_y_min;
                                end
                                    
                            case 2
                                temp_x_max = mean_xy(1) + hight;
                        
                                if temp_x_max > x_max
                                    x_max = temp_x_max;
                                end
                                
                            case 3
                                temp_y_max = mean_xy(2) + hight;
                        
                                if temp_y_max > y_max
                                    y_max = temp_y_max;
                                end
                                
                            case 4
                                temp_x_min = mean_xy(1) - hight;
                        
                                if temp_x_min < x_min
                                    x_min = temp_x_min;
                                end
                                
                        end

                        
                        count_arc = count_arc - 1;
                        
                end


                if count_arc == 0
                    break;
                end

            end
            
        
            x = x_min;

            y = y_min;


            w = x_max - x_min;

            h = y_max - y_min;
            
            
        else
        
            x = min(nodecoords(:,1));

            y = min(nodecoords(:,2));


            w = max(nodecoords(:,1)) - x;

            h = max(nodecoords(:,2)) - y;
            
        end
        
    else
        x = nan;
        y = nan;
        w = nan;
        h = nan;
    end

end