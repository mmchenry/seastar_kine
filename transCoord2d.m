function varargout = transCoord2d(trans_type,varargin)
% Transforms 2d data from one coordinate system to another

% If data_in are coordinates (n x 2):
%  coord_out =  transCoord2d(trans_type,coord_in,tform);
%  trans_type - type of transformation ('global to local' or 
%               'local to global' for coordinates or 'bw G2L' or 'bw L2G' 
%               for binary images)
%  tform - transformation structure for L system defined in G system, 
%          created by defineSystem2d
% If data_in is an image (n x m):
%  im_out =  transCoord2d(im_,trans_type,tform,im);
%  trans_type - type of transformation ('global to local' or 
%               'local to global' for coordinates or 'bw G2L' or 'bw L2G' 
%               for binary images)
%  tform - transformation structure for L system defined in G system, 
%          created by defineSystem2d
%
% Code developed by McHenryLab at UC Irvine

%% Translate inputs

% Inputs for coordinate transformation using tform
if strcmp(trans_type,'G2L') || ...
   strcmp(trans_type,'L2G')
    
    % First input needs to be tform
    tform = varargin{1};
    
    % Second input needs to be the coordinates in L frame
    coord_in = varargin{2};
 
    % Check coordinate dimensions
    if size(coord_in,2)~=2
        error('Coordinates need to given as a n x 2 matrix')
    end
    
    % Check dimensions of tform
    if tform.Dimensionality~=2
        error('Code only handles 2D transformations')
    end
    
% Inputs for image transformation (uses tform)
elseif strcmp(trans_type,'bw L2G')
   
   % First input needs to be tform
    tform = varargin{1};
    
    % Image in roi FOR
    bw_roi = varargin{2};
    
    % binary denoting the roi in global FOR
    bw_G = varargin{3};
    
    % Image in global FOR
    bw_roi_mask = varargin{4};
    
    % Check dimensions of tform
    if tform.Dimensionality~=2
        error('Code only handles 2D transformations')
    end
    
elseif strcmp(trans_type,'ang G2L') 
       
    % First input: origin (1 x 2)
    originG = varargin{1};
    
    % Second input: angle of x-axis (rad)
    angG = varargin{2};
    
    % Third input: global coordinates (n x 2)
    coordG = varargin{3};
     
elseif strcmp(trans_type,'ang L2G')
    
    % First input: origin (1x2)
    originG = varargin{1};
    
    % Second input: angle of x-axis (rad)
    angG = varargin{2};
    
    % Third input: local coordinates (n x 2)
    coordL = varargin{3};
    
    
elseif strcmp(trans_type,'xax G2L') 
       
    % First input: origin (1 x 2)
    originG = varargin{1};
    
    % Second input: point along x-axis (rad)
    xAxisG = varargin{2};
    
    % Third input: global coordinates (n x 2)
    coordG = varargin{3};
     
elseif strcmp(trans_type,'xax L2G')
    
    % First input: origin (1x2)
    originG = varargin{1};
    
    % Second input: angle of x-axis (rad)
    xAxisG = varargin{2};
    
    % Third input: local coordinates (n x 2)
    coordL = varargin{3};
    
else
    error('trans_type not recognized');
end    



%% Transformations
    
% Global to local transformation (coordinates)
if strcmp(trans_type,'G2L')
    
    % Translate
    coord_in(:,1) = coord_in(:,1) - tform.T(3,1);
    coord_in(:,2) = coord_in(:,2) - tform.T(3,2);
    
    % Rotate
    coord_out = [tform.T(1:2,1:2) * coord_in']';
    
    
% Local to global transformation (coordinates)
elseif strcmp(trans_type,'L2G')
    
    % Rotate points
    coord_out = (tform.T(1:2,1:2) \ coord_in')';
    
    % Translate global coordinates wrt origin
    coord_out(:,1) = coord_out(:,1) + tform.T(3,1);
    coord_out(:,2) = coord_out(:,2) + tform.T(3,2);
    
% Local to global transformation (coordinates)
elseif strcmp(trans_type,'bw L2G')
    
    % Coordinate system for im_roi
    R = imref2d(size(bw_roi)); 
    
    % Adjust WorldLimits to restrict transformation to just rotation
    % around center
%     R.XWorldLimits = R.XWorldLimits-mean(R.XWorldLimits);
%     R.YWorldLimits = R.YWorldLimits-mean(R.YWorldLimits);
    
    % Stablize image
    roiRot = imwarp(bw_roi,R,invert(tform),'OutputView',R,...
         'FillValues',0,'SmoothEdges',true);
%     
%     % White out beyond roi
    roiRot(~bw_roi_mask) = 0;
       
    % Black out region outside of roi
    imOut = logical(bw_G.*0);
    
    if sum(bw_G(:))~=(size(roiRot,1)*size(roiRot,2))
        error('Image dimensions do not match');
    end
    
    % Add image
    imOut(bw_G) = roiRot;

    
    
%     % Extract origin
%     origin = tform.T(3,1:2);
%     
%     % Remove translation component to transformation
%     tform.T(3,:) = [0 0 1];
%     
%     % Perform rotation
%     im_rot = imwarp(im_roi,tform,'OutputView',imref2d(size(im_roi)),...
%                   'FillValues',255,'SmoothEdges',true);
%     
%     % Perform translation
%     data_out = imtranslate(im_rot,-origin);
%     
%     
%     aaa=3;
    
  
    
elseif strcmp(trans_type,'ang G2L') 
    
    % Define local system
    S = localSystemAng(originG,angG);
    
    % Transformation
    pts(:,1)    = coordG(:,1)-originG(1);
    pts(:,2)    = coordG(:,2)-originG(2);   
    coordL     = [S'*pts']';
    
    
elseif strcmp(trans_type,'ang L2G') 
    
    % Define local system
    S = localSystemAng(originG,angG);
    
    % Transformation
    pts         = [inv(S)'*coordL']';
    coordG(:,1) = pts(:,1) + originG(1);
    coordG(:,2) = pts(:,2) + originG(2);   
    
elseif strcmp(trans_type,'xax G2L') 
       
    % Define local system
    S = localSystem(originG,xAxisG);

    % Transformation
    pts(:,1)    = coordG(:,1)-originG(1);
    pts(:,2)    = coordG(:,2)-originG(2);   
    coordL     = [S'*pts']';
     
elseif strcmp(trans_type,'xax L2G')
    
    % Define local system
    S = localSystem(originG,xAxisG);

    % Transformation
    pts         = [inv(S)'*coordL']';
    coordG(:,1) = pts(:,1) + originG(1);
    coordG(:,2) = pts(:,2) + originG(2);   
    
else
    
    error('Do not recognize requested transformation')
    
end

%% Set outputs

if strcmp(trans_type,'G2L') || ...
   strcmp(trans_type,'L2G')

    varargout{1} = coord_out;
    
elseif strcmp(trans_type,'bw L2G')
    
    varargout{1}  = imOut;
      
elseif strcmp(trans_type,'ang G2L') || strcmp(trans_type,'xax G2L')
    
    varargout{1}  = coordL;
    
elseif strcmp(trans_type,'ang L2G') || strcmp(trans_type,'xax L2G') 
    
    varargout{1}  = coordG;
    
end



function S = localSystem(pOrigin,pXaxis)
% Defines a transformation matrix for a local coordinate system in an
% inertial frame of reference. 

% Check dimensions of inputs
if size(pXaxis,1)~=1 || size(pXaxis,2)~=2 ||...
   size(pOrigin,1)~=1 || size(pOrigin,2)~=2 
    error('Coordinates must be 1x2 vectors');
end
 
% Define units vectors for x and y axes
xAxis   = [(pXaxis-pOrigin)./norm(pXaxis-pOrigin)];
yAxis   = [-xAxis(2) xAxis(1)];

% Define transformation matrix
S       = [xAxis' yAxis'];


function S = localSystemAng(originG,angG)
% originG   - Origin of coord system in global FOR
% angG      - Angle of x-axis of local FOR wrt to global (rad)
 
% Define axes
xAxis = [cos(angG) sin(angG)];
yAxis = [-xAxis(2) xAxis(1)];
 
% Define transformation matrix
S       = [xAxis' yAxis'];
 
 
% 
% 
% % FUNCTIONS --------------------------
% 
% function ptsT = globalToLocal(tform,coord_in)
% % Assumes column vectors for coordinates
% 
% 
% 
% %pts = [x y];
% 
% % Translate
% coord_in(:,1) = coord_in(:,1) - tform.T(3,1);
% coord_in(:,2) = coord_in(:,2) - tform.T(3,2);
% 
% % Rotate points
% data_out = [tform.T(1:2,1:2) * coord_in']';
% 
% % Extract columns of points
% %xT = ptsT(:,1);
% %yT = ptsT(:,2);
% end
% 
% function ptsT = localToGlobal(tform,pts)
% % Assumes columns vectors for coordinates
% 
% % Check dimensions
% if tform.Dimensionality~=2
%     error('Code only handles 2D transformations')
% end
% 
% % Loop thru columns of coordinates
% i = 1;
%     
%     %pts = [x(:,i) y(:,i)];
%     
%     % Rotate points
%     data_out = (tform.T(1:2,1:2) \ coord_in')';
%     
%     % Translate global coordinates wrt origin
%     data_out(:,1) = data_out(:,1) + tform.T(3,1);
%     data_out(:,2) = data_out(:,2) + tform.T(3,2);
% 
%     
%     clear ptsT pts
% %end
% end