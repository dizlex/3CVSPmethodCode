clear all;close all; clc;
load('FinalCoh_def.mat')
angle=0:180;
V=1000;            % initial SCANNING VELOCITY in m/s
LastVelocityToScan=5000; %the final Velocity to be scanned in m/s
V_step=25;
velocities=V:V_step:LastVelocityToScan; 
[X,Y] =meshgrid(angle,velocities);
    %contour(X,Y,FinalCoh(:,:,mastercount),'--')
    Z=FinalCoh(:,:,25)';
    Z=smoothdata(Z,'gaussian',6); %you need the curveFittingToolBox
    zmin = (min(Z(:))); 
    zmax = (max(Z(:)));
    zinc = (zmax - zmin) / 30; %choose the number of contours
    zlevs = zmin:zinc:zmax;
    %contour(X,Y,FinalCoh(:,:,mastercount)',zlevs)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %uncomment line below for contour plot
    contour(X,Y,Z,zlevs)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %uncomment line below for surf plot
    %surf(X,Y,FinalCoh(:,:,mastercount)')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    shading interp
    colormap turbo
    colorbar
    title(sprintf('Vcoherency VS alpha top depth of:2391 and bottom of 2541'));
    xlabel('Angle (deg)')
    xticks(0:10:180)
    ylabel('Apparent Velocity (m/s)')
    ylim([1750 2750])
    view(2)