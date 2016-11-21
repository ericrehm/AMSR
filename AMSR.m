classdef AMSR < matlab.mixin.Copyable
    
    properties
        resolution
        amsredir
        latlonfn
        latgrid
        longrid
        TreeRoot
        iceconc
        yyyymmdd
    end
    
    methods
        
        function obj = AMSR(resolution)
            
            obj.resolution = resolution;
            obj.TreeRoot = [];
            
            switch (obj.resolution)
                
                case 25
                    
                case 6
                    
                    % 6.25 km Ice: get grid and setup k-d tree
                    obj.amsredir = '/Volumes/taku-njall/AMSRE/';
                    obj.latlonfn = [obj.amsredir 'LongitudeLatitudeGrid-n6250-Arctic.hdf'];
                    obj.latgrid = double(hdfread(obj.latlonfn, 'Latitudes'));
                    obj.longrid = double(hdfread(obj.latlonfn, 'Longitudes'));

                    
                case 3
                    % 3.125 km Ice: get grid and setup k-d tree
%                     obj.amsredir = '/Volumes/taku-njall/AMSR2/';
%                     %         amsredir = '/Users/ericrehm/Documents/ICE/AMSR2/';
%                     obj.latlonfn = [obj.amsredir 'LongitudeLatitudeGrid_3.125km_Arctic.nc'];
%                     obj.latgrid = double(ncread(obj.latlonfn, 'latitude'));
%                     obj.longrid = double(ncread(obj.latlonfn, 'longitude'));
                    
                    % 3.125 km Ice: get grid and setup k-d tree
                    obj.amsredir = '/Volumes/taku-njall/AMSRE/';
                    obj.latlonfn = [obj.amsredir 'LongitudeLatitudeGrid-n3125-NorthWestPassage.hdf'];
                    obj.latgrid = double(hdfread(obj.latlonfn, 'Latitudes'));
                    obj.longrid = double(hdfread(obj.latlonfn, 'Longitudes'));
                    
            end
        end  % AMSR
        
        function [TreeRoot] = getKdTree(obj)
            
            if (obj.resolution == 25)
                TreeRoot = [];
            else
                ReferencePts = [obj.latgrid(:) obj.longrid(:)];
                [tmp, tmp, TreeRoot] = kdtree( ReferencePts, []);
            end
            obj.TreeRoot = TreeRoot;
            
        end  % getKdTree
        
        function Ice = getIce(obj, dateStr, lat, lon)
            
            if (isempty(obj.TreeRoot))
                obj.getKdTree();
            end
            
            [Ice bFound] = obj.getIceAny(dateStr, lat, lon);
%             disp([Ice, bFound]);
            
        end  % getIce
        
        function [Ice bFound] = getIceAny(obj, dateStr, lat, lon)
            
            % Get sea ice concentration
            bFound = true;
            
            % Get NSIDC sea ice concentration
            date2 = datestr(datenum(dateStr),30);
            yyyymmdd = date2(1:8);
            dvec = datevec(datenum(dateStr));
            yyyy = date2(1:4);
            
%             disp(obj.resolution);
            switch (obj.resolution)
                case 3
                    ice = obj.readIce(dvec(1), dvec(2), dvec(3));
                    
                case 6
                    ice = obj.readIce(dvec(1), dvec(2), dvec(3));
                    
                case 25
                    command = sprintf('~/bin/getIce %s %f %f', yyyymmdd, lat, lon);
                    [stat result] = system(command);
                    Ice = str2double(result)*100;
                    if (Ice == -999)
                        Ice = NaN;
                        qIce = 1;
                        bFound = false;
                    end
                    return
            end
            
            
            % If we found an 6 km AMSR ASI file, find the closest point
            if (bFound)
                
                if (lon < 0)
                    TestPoint = [lat lon+360];
                else
                    TestPoint = [lat lon];
                end
                % Find the row index of the closest point in ReferencePts
                [ ClosestPtIndex, DistB, TreeRoot ] = kdtreeidx([], TestPoint, ...
                    obj.TreeRoot);
                Ice = ice(ClosestPtIndex);
                if (~isnan(Ice))
                    qIce = 0;
                else
                    qIce = 1;
                end
            else
                Ice = NaN;
                qIce = 1;
            end
        end  %getIceAny
        
        function iceconc = readIce(obj, yyyy, mm, dd)
            
            yrday = yeardayEAD(dd,mm,yyyy,0,0,0);
            
            yyyymmdd = sprintf('%04d%02d%02d', yyyy, mm, dd);

            switch (obj.resolution)

                case 3
%                     icefn  = sprintf('%s%04d/Arc_%s_res3.125_pyres.nc', obj.amsredir,yyyy,yyyymmdd);
%                     if (~exist(icefn, 'file'))
%                         bFound = false;
%                         Ice = NaN;
%                         return
%                     end
%                     icevar =  'sea_ice_concentration';
%                     iceconc = double(ncread(icefn, icevar));

                    %asi-AMSR2-n3125-20161011.hdf
                    icefn  = sprintf('%s%04d/%03d/asi-AMSR2-n3125-%s.hdf', obj.amsredir,yyyy, yrday, yyyymmdd);
                    if (~exist(icefn, 'file'))
                        bFound = false;
                        Ice = NaN;
                        return
                    end
                    icevar =  'ASI Ice Concentration';
                    iceconc = double(hdfread(icefn, icevar));

                case 6
                    icefn  = sprintf('%s%04d/%03d/asi-n6250-%s-v5.hdf', obj.amsredir,yyyy, yrday, yyyymmdd);
                    if (~exist(icefn, 'file'))
                        icefn  = sprintf('%s%04d/%03d/asi-SSMIS17-n6250-%s-v5.hdf', obj.amsredir,yyyy, yrday, yyyymmdd);
                        if (~exist(icefn, 'file'))
                            icefn  = sprintf('%s%04d/%03d/asi-AMSR2-n6250-%s-v5.hdf', obj.amsredir,yyyy, yrday, yyyymmdd);
                            if (~exist(icefn, 'file'))
                                bFound = false;
                                iceconc = NaN;
                                return
                            end
                        end
                    end
                    icevar =  'ASI Ice Concentration';
                    iceconc = double(hdfread(icefn, icevar));
            end
            obj.iceconc = iceconc;
            obj.yyyymmdd = yyyymmdd;
            
        end  % getData
        
        function [cs1, h1] = mapIce(obj, level)
            if (nargin < 2)
                level = 10;
            end
            [cs1, h1] = m_contour(obj.longrid, obj.latgrid, obj.iceconc, [level level]);
            set(h1, 'color', 'r');
            set(h1, 'Fill', 'on', 'FaceColor', [255 234 234]/255)
            [cs, h2] = m_contour(obj.longrid, obj.latgrid, obj.iceconc, [50 50]);
            set(h2, 'color', 0.5*[1 1 1]);
            [cs, h3] = m_contour(obj.longrid, obj.latgrid, obj.iceconc, [100 100]);
            set(h3, 'color', 'k');
            set(h3, 'Fill', 'on', 'FaceColor', [255 214 214]/255)
            legend([h1 h2 h3], sprintf('%d%%', level), '50%', '100%');
           
        end  % mapContour

        function [h1] = mapIceSurf(obj)
            if (nargin < 2)
                level = 10;
            end
            h1 = m_pcolor(obj.longrid, obj.latgrid, obj.iceconc);
           
        end  % mapContour
    end
end  % classdef AMSR

