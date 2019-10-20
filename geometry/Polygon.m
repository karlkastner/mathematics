% 2015-09-29 14:36:03.575085215 +0200
% Bart Vermeulen
%% Simple 2D polygon class
%%
%%   Polygon properties:
%%       x - x coordinates of polygon
%%       y - y coordinates of polygon
%%       nnodes - number of nodes in the polygon
%%
%%   Polygon methods:
%%       in - checks whether given points lie inside, on the edge, or outside of the polygon
%%       area - returns the area of the polygon
%%       centerline - computes the centerline of the river
%%       iscw - check whether polygon is clockwise
%%       reverse - reverse the order of the polygon

%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

classdef Polygon < handle
    
    properties
        x % x coordinate of the polygon. should be a finite vector with at least three elements and equal first and last coordinates.
        y % y coordinate of the polygon. should be a finite vector with at least three elements and equal first and last coordinates.
    end % properties
    properties(Dependent, SetAccess=private)
        nnodes % number of nodes in the polygon
    end
    methods
        % Constructor
        function obj=Polygon(x,y)
            if nargin==0
                obj.x=[];
                obj.y=[];
            end
            if nargin==2
                obj.x=x;
                obj.y=y;
            end
        end
        % set and get
        function set.x(obj,val)
            assert(isempty(val) || (isvector(val) && isnumeric(val)))
            if isrow(val)
                val=val';
            end
            obj.x=val;
        end % set.x
        
        function val=get.x(obj)
            val=obj.x;
        end % get.x
        
        function set.y(obj,val)
            assert(isempty(val) || (isvector(val) && isnumeric(val)))
            if isrow(val)
                val=val';
            end
            obj.y=val;
        end % set.y
        
        function val=get.y(obj)
            val=obj.y;
        end % get.y
        
        function val=get.nnodes(obj)
            val=numel(obj.x);
        end % get.nnodes
        
        function varargout=isin(obj,x,y)
            %ISIN True for points inside the Polygon.
            %   IN = ISIN(Polygon,X,Y) returns a matrix IN the size of X and Y.
            %   IN(p,q) = 1 if the point (X(p,q), Y(p,q)) is either strictly inside or
            %   on the edge of Polygon
            %
            %   [IN ON] = ISIN(Polygon,X,Y) returns a second matrix, ON, which is
            %   the size of X and Y.  ON(p,q) = 1 if the point (X(p,q), Y(p,q)) is on
            %   the edge of Polygon; otherwise ON(p,q) = 0.
            %
            %   This function uses the matlab builtin function INPOLYGON
            [varargout{1}, out]=inpolygon(x,y,obj.x,obj.y);
            if nargout>1
                varargout{2}=out;
            end % if
        end % isin
        
        function val=area(obj)
            %AREA Area of the Polygon.
            %   area(Polygon) returns the area of Polygon.
            %
            %   The polygon edges must not intersect.  If they do, area
            %   returns the absolute value of the difference between the clockwise
            %   encircled areas and the counterclockwise encircled areas.
            %
            %   This function uses the matlab builtin function POLYAREA
            val=polyarea(obj.x,obj.y);
        end
        function C=centerline(objin,method)
            met='voronoi';
            if nargin>1
                assert(ischar(method),'Polygon:centreline_method_not_char','Second input variable method must be of type char')
                assert(any(strcmp(method,{'voronoi','skeleton','segment_voronoi'})),'Polygon:centreline_invalid_method','Second input must be ''voronoi'', ''skeleton'' or ''segment_voronoi''')
                met=method;
            end
            switch met
                case 'segment_voronoi'
                    if exist('segment_voronoi','file')~=3
                        warning('Polygon:segmentVoronoiUnavailable','Segment voronoi function not found. Falling back to voronoi method')
                        met='voronoi';
                    end
                case 'skeleton'
                    if exist('skeleton','file')~=3
                        warning('Polygon:SkeletonUnavailable','Skeleton function not found. Falling back to voronoi method')
                        met='voronoi';
                    end
            end
            C(numel(objin),1)=Centreline;
            C=reshape(C,size(objin));
            switch met
                case 'voronoi'
                    for co=1:numel(objin)
                        obj=objin(co);
                        dt=delaunayTriangulation(obj.x(2:end),obj.y(2:end),[(1:obj.nnodes-1)' [(2:obj.nnodes-1)';1]]); % Make a triangulation of the polygon points
                        [V,R]=dt.voronoiDiagram; % make a voronoi diagram
                        isIn=find(obj.isin(V(:,1),V(:,2))); % find diagram vertices within the polygon (centreline points)
                        R=cellfun(@(x) [x' x([2:end 1])'],R,'uniformOutput',false); % Make edges from polygon vertices
                        R=cellfun(@(x) x(all(ismember(x,isIn),2),:),R,'uniformOutput',false); % Remove edges not in polygon
                        R=vertcat(R{:}); % put all edges below each other
                        [C(co).conn,C(co).nodes]=unmesh([V(R(:,1),1),V(R(:,1),2),V(R(:,2),1),V(R(:,2),2)]); % Create graph
                    end
                case 'skeleton' % switch met
                    for co=1:numel(objin)
                        obj=objin(co);
                        if obj.iscw(), obj.reverse(); end;
                        pol.x=obj.x();
                        pol.y=obj.y();
                        [hedg,isbis]=skeleton(pol);
                        [C(co).conn,C(co).nodes]=unmesh(hedg(isbis,:));
                    end
                case 'segment_voronoi' % switch met
                    for co=1:numel(objin)
                        obj=objin(co);
                        sv=segment_voronoi(obj.x(2:end),obj.y(2:end),[1:obj.nnodes-1; [2:obj.nnodes-1 1]]);
                        sv(cellfun(@isstruct,sv))=[]; % remove rays from sv
                        [fin, fon]=cellfun(@(x) (obj.isin(x(1,:),x(2,:))),sv,'uniformoutput',false); % find sv segment with points in the polygon and on the edges
                        sv(cellfun(@(in,on) ~all(in) || any(on),fin,fon))=[]; % remove segments with not all points inside the polygon or any point on the polygon edge 
                        if isempty(sv)
                            return
                        end
                        sv=cellfun(@(x) [x(1,1:end-1); x(2,1:end-1);x(1,2:end); x(2,2:end)]',sv,'uniformoutput',false); % get endpoints of segments
                        sv=vertcat(sv{:}); % put the below each other
                        [C(co).conn,C(co).nodes]=unmesh(sv);
                    end
            end % switch met
        end % function centerline
        function out=iscw(obj)
            % iscw Checks whether polygon is clockwise
            %   TF=iscw(Polygon) returns true if polygon is clockwise and
            %   false if polygon is counterclockwise
            xl=obj.x;
            yl=obj.y;
            out=sum((xl(2:end)-xl(1:end-1)).*(yl(2:end)+yl(1:end-1)))>0;
        end
        function reverse(obj)
            % Changes the order of the polygon
            obj.x=flipud(obj.x);
            obj.y=flipud(obj.y);
        end
        function varargout=plot(obj,varargin)
            % Plot Plots the Polygon(s)
            allpl=reshape([reshape({obj(:).x},1,[]); reshape({obj(:).y},1,[])],[],1);
            hp=plot(allpl{:},varargin{:});
            if nargout==1
                varargout{1}=hp;
            end
            axis equal
        end % plot
    end % methods
end % Polygon
