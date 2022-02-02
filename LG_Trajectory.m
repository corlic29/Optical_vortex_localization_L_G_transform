clear all
clc

%% Singular point localization based on Laguerre - Gaussian transform
%{
Scritps detects the position of a single singular point based on the images
representing intensity distribution

Authors: Mateusz Szatkowski, Emilia Burnecka, Hanna Dyla
Contact: mateusz.szatkowski@pwr.edu.pl

Script additionally uses functions:

intersections
Douglas Schwarz (2022). Fast and Robust Curve Intersections (https://www.mathworks.com/matlabcentral/fileexchange/11837-fast-and-robust-curve-intersections), MATLAB Central File Exchange. Retrieved January 20, 2022.
%}

%Choose files
[FileName,PathName] = uigetfile({'*.jpg; *.bmp; *.png; *.tif'},'Select file','Multiselect','on');
[~,num_of_files]=size(FileName);

%Read files and save in a structure
for z=1:1:num_of_files
a=imread([PathName strjoin(FileName(z))]);

%Convert to grayscale
A{z}=mat2gray(a);
end

%% Bandwith determination
% Manual
core_diameter=20; %[px]
w=core_diameter; 


%% Localization algorithm initialization

% Crop (optional, reduces computation time)
crop_start=100;
crop_end=180;


% Polar coordination system
pX=round(2*(crop_end-crop_start));
pY=round(2*(crop_end-crop_start));

vector_x=-1:2/(pX-1):1;
vector_y= -1:2/(pY-1):1;
[X, Y]=meshgrid(vector_x, vector_y);

[alpha, r]=cart2pol(X, Y); %Polarne

% LG function definition
LG = (1i*pi.^2*w.^4).*(r.*exp(-pi.^.2*r.^2*w.^2).*exp(1i.*alpha));


% Pre-definition of variable
positions = [];

%% Detection
for i=1:1:num_of_files
%Load the image
    im=A{i};
    im=(im(:,:,1));
    %im=flipud(im);
    I=im;
    im=imcrop(im,[crop_start, crop_start, crop_end, crop_end]);
    
    %Improve the contrast (optional)
    im2=imadjust(im,[0.1 0.5]);

    %Pseudo-complex amplitude  
    im_complex=conv2(im2,LG,'same');
    im_complex_II=im_complex.*conj(im_complex);
   
    real_im=real(im_complex);
    imag_im=imag(im_complex);
    
    %Surface approximation
    c=contourc(real_im,[0,0]); %Zero level contour
    c(:,~c(1,:)) = NaN;
    
    c2=contourc(imag_im,[0,0]); %Zero level contour
    c2(:,~c2(1,:)) = NaN;
    
    %Crossection of two contours
    [x,y,~,~]=intersections(c(1,2:end),c(2,2:end),c2(1,2:end),c2(2,2:end),true);
    
    %Criterion for correct vortex posistion (lowest intensity)
    x_round=round(x); 
    y_round=round(y);
    
    [length_,~]=size(x_round);
    
    %Lowest intensity criterion
    for j=1:1:length_
    points(j)=im(y_round(j),x_round(j));
    end
    [~, idx]=min(points);
    
    %Exchange of the point
    x_round=x_round(idx);
    y_round=y_round(idx);
    points=[];
    
    %positions=[positions; x_round, y_round];
    
    x=x(idx);
    y=y(idx);

    positions=[positions; x+crop_start, y+crop_start];
    
 
end

%% Display trajectory
figure()
set(gcf,'color','w');

plot(positions(:,1),positions(:,2),'.k','MarkerSize',5)

axis([0 length(A{1}) 0 length(A{1})])
axis([113 163, 112 162])
xlabel('Position x [px]','FontSize',10)
ylabel('Position y [px]','FontSize',10)
title('Trajectory')

%% Save (optional)
%{
xlswrite('Trajectory.xls',positions,1);
%}