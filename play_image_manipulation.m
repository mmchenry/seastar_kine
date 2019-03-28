function play_image_manipulation
% Playing with code for manipulating and anlyzing images



%% Very simple image
% Based on examples at https://blogs.mathworks.com/steve/2013/08/28/introduction-to-spatial-referencing/

if 1
    
% Create very simple image
I = reshape(9:-1:1,3,3);

% Crop image
I2 = imcrop(I,[1 1 1 1]);

RI = imref2d(size(I));
RI.XWorldLimits = [3 6];
RI.YWorldLimits = [6 9];

if 1
    %figure,
    subplot(2,2,1)
    image(I); axis square
    title('Intrinsic coordinates')
    
    subplot(2,2,2)
    imshow(I,RI,[],'InitialMagnification','fit');
    title('World coordinates')
    
    subplot(2,2,3)
    image(I2); axis square
    title('Intrinsic coordinates')
    
end

% Convert world coordinates to intrinsic coordinates
[intrinsicX,intrinsicY] = worldToIntrinsic(RI,4.6,7.3);

% This finds the closest row and column in the intrinsic system from world
% coordinates
[r,c] = worldToSubscript(RI,4.6,7.3);

% Here is the pixel value at that coordinate
I(r,c);

end


%% Tansformation of a checkerboard
% Based on: https://www.mathworks.com/help/images/matrix-representation-of-geometric-transformations.html

if 0

do_plot = 1;

%Create a checkerboard image that will undergo transformation.
cb = checkerboard(4,2);

%Also create a spatial reference object for the image.
cb_ref = imref2d(size(cb));

% create a flat background image
background = zeros(150);

% Visualize
figure
subplot(2,3,1)
imshowpair(cb,cb_ref,background,imref2d(size(background)))
subplot(2,3,4)
imshow(cb,cb_ref,[],'InitialMagnification','fit');


% Translation matrix: shift the image horizontally by 100 pixels.
T = [1 0 0;0 1 0;100 0 1];
tform_t = affine2d(T);

% Rotation matrix, to rotate by 30 deg (NOTE: Left-handed sign of rotation)
R = [cosd(30) sind(30) 0;-sind(30) cosd(30) 0;0 0 1];
tform_r = affine2d(R);


% Translation, followed by rotation
TR = T*R;
tform_tr = affine2d(TR);
[out,out_ref] = imwarp(cb,cb_ref,tform_tr);

% Display results
subplot(2,3,2)
imshowpair(out,out_ref,background,imref2d(size(background)))
title('Translation, then rotation')
subplot(2,3,5)
imshow(out,out_ref,[],'InitialMagnification','fit');

% Rotation, Followed by Translation
RT = R*T;
tform_rt = affine2d(RT);
[out,out_ref] = imwarp(cb,cb_ref,tform_rt);

% Display results
subplot(2,3,3)
imshowpair(out,out_ref,background,imref2d(size(background)))
title('Rotation, followed by translation')
subplot(2,3,6)
imshow(out,out_ref,[],'InitialMagnification','fit');

end


%% Test modification of image

if 0
    
    % Create a white background image
    bg     = ones(500);
    bg_ref = imref2d(size(bg));
    
    % Create a checkerboard image that will undergo transformation.
    cb = checkerboard(15,2);
    
    % Create a spatial reference object for the image.
    cb_ref = imref2d(size(cb));
    
    % Create transformation object
    tform = imMod('rot, then trans',15,[50 75]);
    
    % Transform image, output image in bg FOR
    [im_out,ref_out] = imwarp(cb,cb_ref,tform,'FillValues',1,'OutputView',bg_ref);
    
    % Display output
    figure;
    imshow(im_out,ref_out,[],'InitialMagnification','fit');
end


%% Create virtual movie of checkboard kinematics

% Rate of rotation (deg/frame)
rot_rate = 1.5*3;

% Translational velocity (pix/frame)
vel = [3.5 6].*2;

% Vectir of frame numbers
frames = 1:50;

% Create a white background image
bg     = ones(1000);
bg_ref = imref2d(size(bg));

% Create a checkerboard image that will undergo transformation.
cb = checkerboard(15,2);

% Create a spatial reference object for the image.
cb_ref = imref2d(size(cb));

% Offset of center 
cntr_off = [size(cb,1)/2 size(cb,1)/2];

% Loop thru frames
for i = 1:length(frames)
    
    if i==1
        angl = 0;
        Body.x = cntr_off(1)+ 2*size(cb,1);
        Body.y = cntr_off(2)+ 2*size(cb,1);
    else
        % Store data for comparison later
        angl(i,1)    = angl(end) + rot_rate;
        Body.x(i,1)  = Body.x(end) + vel(1);
        Body.y(i,1)  = Body.y(end) + vel(2);
    end
    
    % Create transformation object
    tform = imMod('rot, then trans',angl(end),[Body.x(i) Body.y(i)]-cntr_off);
    
    % Transform image, output image in bg FOR
    [im_out,ref_out] = imwarp(cb,cb_ref,tform,'FillValues',1,'OutputView',bg_ref);
    
    % Store image
    mov(:,:,i) = im_out;
    
end

clear tform im_out ref_out cb_ref frames


% Add noise to center point
nAmp = 3;
Body.x = Body.x+nAmp.*rand([length(Body.x) 1]);
Body.y = Body.y+nAmp.*rand([length(Body.y) 1]); 

% Play movie
if 0
    figure
    for i = 1:size(mov,3)
        imshow(mov(:,:,i),'InitialMagnification','fit');
        hold on
        plot(Body.x(i),Body.y(i),'ro')
        hold off
        title(['Frame ' num2str(i)])
        pause(0.1)
    end
end


%% Analyze virtual movie with image registration

 % Initialize image registration parameters
 [optimizer, metric]  = imregconfig('monomodal');
 optimizer.MaximumStepLength = 5e-4;
 optimizer.MaximumIterations = 1500;
 optimizer.RelaxationFactor  = 0.2;


% Radius of roi
roi_r = size(cb,1)/2*2.5;

% Reference frame
im0 = mov(:,:,1);

% roi rectangle
rect = [Body.x(1)-roi_r Body.y(1)-roi_r 2*roi_r 2*roi_r];

% x0 = 

% Reference roi image
im_roi0 = imcrop(im0,rect);

% Last transformation
Body.tform = imMod('rot',0);

figure;

% Loop thru movie
for i = 1:size(mov,3)
    
    % Current frame
    im = mov(:,:,i);
    
    % roi rectangle
    rect = [Body.x(i)-roi_r Body.y(i)-roi_r 2*roi_r 2*roi_r];
    
    % Current roi image
    im_roi = imcrop(im,rect);
    
    % Coordinate system for im_roi
    R = imref2d(size(im_roi));
    
    % Adjust image to last transformation
    im_roi0curr = imwarp(im_roi0,invert(Body.tform(end)),'OutputView',R,...
                      'FillValues',1,'SmoothEdges',true);
    
    % Transformation object for displacement since first frame
    tform_disp = imregtform(im_roi,im_roi0,'rigid',optimizer,metric,...
        'InitialTransformation',Body.tform(end));
        
    
% Bblwo code was an ill-fated attempt to replicate above without using
% 'InitialTransformation', which triggered a bug when run on lionfish
% analysis
%     % Trasnform im_roi to match last frame
%     im_roiM = imwarp(im_roi0,invert(tform_disp),'OutputView',R,...
%                       'FillValues',1,'SmoothEdges',true);
%     
%     % Transformation object for displacement since last frame
%     tform_dispM = imregtform(im_roi,im_roiM,'rigid',optimizer,metric);
    

    Body.tform(i) = tform_disp;
    
%     imStable = imwarp(im_roi,tform_disp,'OutputView',R,...
%                       'FillValues',1,'SmoothEdges',true);
                  

    % Transform coordinates
    [xL,yL] = transformPointsForward(invert(tform_disp),roi_r,roi_r);
    [xL2,yL2] = transformPointsForward(invert(tform_disp),roi_r,1.5*roi_r);

    % NOTE: Looks tformfwd does the same thing

    if 0
    % Update display
    subplot(2,3,1)
    imshow(im_roi0,imref2d(size(im_roi0)),[],'InitialMagnification','fit');
    title('Reference roi')
    
    subplot(2,3,4)
    imshow(im_roi,imref2d(size(im_roi)),[],'InitialMagnification','fit');
    title('Current roi')
    
    subplot(2,3,2)
    imshow(im_roi0curr,imref2d(size(im_roi0curr)),[],'InitialMagnification','fit');
    hold on
    plot(xL,yL,'r+')
    hold off
    title('Reference, adjusted to last')
    
    subplot(2,3,5)
    imshow(imStable,imref2d(size(imStable)),[],'InitialMagnification','fit');
    title('Current roi, stabilized')
    
    subplot(2,3,3)
    imshow(im,imref2d(size(im)),[],'InitialMagnification','fit')
    hold on
    plot([rect(1) rect(1)+rect(3) rect(1)+rect(3) rect(1) rect(1)],...
         [rect(2) rect(2) rect(2)+rect(4) rect(2)+rect(4) rect(2)],'g-');
     plot(Body.x(i)-roi_r+xL,Body.y(i)-roi_r+yL,'+r')
    hold off
    title('Whole frame')
    
    subplot(2,3,6)
    warning off
    imshowpair(im_roi0,imStable)
    warning on
    title('Image pair')
    end
    
    if 1
        imshow(im,imref2d(size(im)),[],'InitialMagnification','fit')
        hold on
        plot([rect(1) rect(1)+rect(3) rect(1)+rect(3) rect(1) rect(1)],...
            [rect(2) rect(2) rect(2)+rect(4) rect(2)+rect(4) rect(2)],'k-');
        plot(Body.x(i)-roi_r+xL,Body.y(i)-roi_r+yL,'or')
        plot([Body.x(i)-roi_r+xL Body.x(i)-roi_r+xL2],...
             [Body.y(i)-roi_r+yL Body.y(i)-roi_r+yL2],'-r')
        hold off
        title('Whole frame')
    end

    pause(0.3)
end



function tform = imMod(action,varargin)
% Modify an image, according to translation and rotation (in deg), 


% Image rotation
if strcmp(action,'rot')
    
    % Rotation in degrees
    rot_deg = varargin{1};
    
    % Rotation matrix
    R = [cosd(rot_deg) sind(rot_deg) 0;-sind(rot_deg) cosd(rot_deg) 0;0 0 1];
       
     % Affine object
    tform = affine2d(R);
    
% Image translation
elseif strcmp(action,'trans')
    
    % Translation
    displ = varargin{1};
    
    % Translation matrix
    T = [1 0 0;0 1 0;displ(1) displ(2) 1];
    
     % Affine object
    tform = affine2d(T);
    
% Rotate, then translate the image
elseif strcmp(action,'rot, then trans')
    
    % Rotation in degrees
    rot_deg = varargin{1};
    
    % Rotation matrix
    R = [cosd(rot_deg) sind(rot_deg) 0;-sind(rot_deg) cosd(rot_deg) 0;0 0 1];
    
    % Translation
    displ = varargin{2};
    
    % Translation matrix
    T = [1 0 0;0 1 0;displ(1) displ(2) 1];
    
     % Affine object
    tform = affine2d(R*T);
    
% Translate, then rotate the image
elseif strcmp(action,'trans, then rot')
    
     % Rotation in degrees
    rot_deg = varargin{1};
    
    % Rotation matrix
    R = [cosd(rot_deg) sind(rot_deg) 0;-sind(rot_deg) cosd(rot_deg) 0;0 0 1];
    
    % Translation
    displ = varargin{2};
    
    % Translation matrix
    T = [1 0 0;0 1 0;displ(1) displ(2) 1];
    
    % Affine object
    tform = affine2d(T*R);

else
    error(['Do not understand ' action]);
end










