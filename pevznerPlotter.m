%plotter of the data
close all;
%Pevzner Plot v2.0

 %% Plot Fig 1.a (easiest way)
 clear all;
 clc;
 load seisX
 load seisY
 hist=1;
 if hist ==1
 figure('Name','Histogram of Hx');
 histogram(seisX);
 figure('Name','Histogram of Hy');
 histogram(seisZ);
 figure('Name','Seismograms');
 end
 r=size(seisX);
 rr=size(seisZ);
%INPUT FROM USER
 dt=0.002;
 depthRecord=1145;
 dz=50;

 [X,Y] =meshgrid(1:r(1) , 0:dt:dt*(r(2)-1)); %depthRecord(1)+dz:dz:depthRecord(1)+dz*70
 title('another way to see the original seismograms')
 subplot(1,3,1)
 surf(X,Y,seisX')
 shading interp
 colormap bone
 colorbar
 title('Hx')
 xlabel('receiver')
 ylabel('time[s]')
 %caxis([-25557 18579])
 view(2) %to view it as 2D
 view([0 -90])

 subplot(1,3,2)
 [X,Y] =meshgrid(1:rr(1) , 0:dt:dt*(r(2)-1));
 surf(X,Y,seisZ')
 shading interp
 colormap gray
 colorbar
 title('Hy')
 xlabel('receiver')
 ylabel('time[s]')
 %caxis([-25557 32500])
 view(2) %to view it as 2D
 view([0 -90])
%% Plot Fig 1.b
seisRot = pevznerRotation(seisX,seisZ,45);
subplot(1,3,3)
surf(X,Y,seisRot')
 shading interp
 colormap gray
 colorbar
 title('Rotated at 45 deg')
 %xticks(1:10:70)
 %yticks(0:.1:3)
 xlabel('receiver')
 ylabel('time[s]')
 %caxis([-25557 32500])
 view(2) %to view it as 2D
 view([0 -90])