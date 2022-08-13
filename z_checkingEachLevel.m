% this script allows the user to select a specific depth interval
% plot it and save the result in a .txt file

clear all;close all; clc;
load('FinalCoh_def.mat')
clear all;close all; clc;
load('seisXVfrom2600to3200Vstep20AngStep1WinSize=006Num_Rec=11DepthInterval=2940to3140FinalCoh.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Select plotting bounds
minVelocityToPlot=2600;
maxVelocityToPlot=3100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [minDistance, indexOfMin] = min(abs(minVelocityToPlot-velocities));
    [minDistance,indexOfMax]=min(abs(maxVelocityToPlot-velocities));
    [X,Y] =meshgrid(angle,velocities(indexOfMin:indexOfMax));
    %contour(X,Y,FinalCoh(:,:,mastercount),'--')
    Z=FinalCoh(:,indexOfMin:indexOfMax,1)'; %index of Min velocity, last dimension is slide number
    %Z=circshift(Z,-10,2);
    Z=smoothdata(Z,'gaussian',15); %you need the curveFittingToolBox
    
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
    %surf(X,Y,Z(:,:,1))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    shading interp
    colormap turbo % selecting the color map
    %colorbar
    hCB=colorbar;
    %hCB.Title.String='Title'
    [n1,n2]=max(Z); %obtaining the [valueVector,PositionVector]
    maxColorbar=max(n1); % obtaining the max value
    set(hCB.XLabel,{'String','Rotation','Position'},{'Coherency',90,[-1 maxColorbar*0.50]}) %setting the colorbar label 
    %title(sprintf('Vcoherency VS alpha top depth of:2300 and bottom of 2500'));
    title('(b)','FontSize',14)
    xlabel('Angle (deg)')
    xticks(0:10:180)
    ylabel('Apparent Velocity (m/s)')
    ylim([minVelocityToPlot maxVelocityToPlot])
    view(2)
    %%%%%%%%%%%%%%%%%%%%%%%%

%% getting the slow event coordinates in angle and velocity
    [maxFinCoh, AngIndexOfFinCoh] = max(FinalCoh);
            %velIndex
 [TrueMax, indTrueMax]=max(maxFinCoh); %obtaining the [valueVector,PositionVector]
 maxCohVel_SlowEvent=velocities(indTrueMax) %the velocity coordinate of the max Coherency Value for the slow velocity
 maxCohAng_SlowEvent=AngIndexOfFinCoh(indTrueMax) %the angle coordinate of the max Coherency value for the slow event
 SlowEventCoh=FinalCoh(AngIndexOfFinCoh(indTrueMax),indTrueMax); %theValueAt (maxCohAng_slow,maxCohVel_slow)
%% getting the fast event coordinates in angle and velocity
 [maxFinCoh, AngIndexOfFinCoh] = max(FinalCoh(:,floor(length(velocities)/2):end)); %obtaining the [valueVector,PositionVector]
            %velIndex
 [TrueMax, indTrueMax]=max(maxFinCoh); %obtaining the [valueVector,PositionVector]

 maxCohVel_FastEvent=velocities(floor(length(velocities)/2)+indTrueMax) %the velocity coordinate of the max Coherency Value for the slow velocity
 maxCohAng_FastEvent=AngIndexOfFinCoh(indTrueMax) %the angle coordinate of the max Coherency value for the slow event
 FastEventCoh=FinalCoh(AngIndexOfFinCoh(indTrueMax),indTrueMax); %theValueAt (maxCohAng_slow,maxCohVel_slow)
 writematrix(FinalCoh,'FigureFileName.txt') %saving the result in a txt file
