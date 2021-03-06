clear all;clc;close all;
%CHECK THAT YOU HAVE THE seisX.mat and seisZ.mat saved at same directory
%                        you need the curveFittingToolBox installed to plot
%                        you need the exportToPPTX at same directory
%This script is the MATLAB coding of the Pevzner et al.2011 method:
%Estimation of azimuthal anisotropy from VSP data using Multicomponent S-wave velocity analysis
%------------------------------------------------------------------------
Num_Rec=14;    % total number of receivers in window
plotSeis=0;    %if you want to plot seismogram of choice plotSeis=1
%------------------------------------------------------------------------
%% Loading the data
tStart = tic; %start chronometer to know who long the program takes to run
%put the datafiles here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
filenameX='vspSurvey_xComponent_stack_new.mat'; % matlab does not differentiate capital letters
filenameZ='vspSurvey_yComponent_stack_new.mat';
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
load(filenameX);
load(filenameZ);
seisZ=seisY; %this is because the code was written originally for a Z which didnt represent depth
%checking that both files contain the same number of traces
if size(seisZ) ~= size(seisX)
    return
end
filenameX=split(filenameX,'.');%getting the part before the .mat
filenameX=filenameX{1}; %saving before the .mat
diary(sprintf('Out_%s.txt',filenameX))%saving the command window to a txt file


%% INPUTS FOR PEVZNER METHOD-----------------------------------------------------------------
dt=0.001;          %T(1,length(T))/length(T); %getting the time step in time domain 
dz=15;             %space between receivers at depth (z-axis) in meters
z_0=4105.8;        %initialDepth of first receiver in feet
z_0=round(z_0.*(1/3.281));     %converting to meters
angle=0;           %the first angle in degrees to rotate the seismogram
angleStep=1;       %the angle step used by Pevzner (no decimal points allowed)
lastAngle=180;     %last angle to be scanned
V=2250;            % initial SCANNING VELOCITY in m/s
LastVelocityToScan=4000; %the final Velocity to be scanned in m/s
V_step=25;         %note that it may be that V_steps < 10 produce time differences less than dt sampling interval. They may not affect the coherency)
winSizeInSec=0.03; %half the window size in seconds has to be smaller than or equal to t_0

%% variables for pevzner method--------------------------------------
[TotalNumOfReceivers,nt] = size(seisX); %must be equal to seisZ
T=zeros(Num_Rec,nt);                    %space for the time array
totalTime=dt*nt;                        %dt is the time step, nt is the number of time samples
TotalDepth=z_0+dz*TotalNumOfReceivers;  %Total depth of the survey Num_Rec is the number of traces
timeRecord=0:dt:totalTime-dt;           %The row vector containing all the [s] points
depthRecord=z_0:dz:TotalDepth-dz;       %The row vector  containing all the [z depth] points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_0Step=winSizeInSec/2;   %The reference time Step so that each window overlaps by half
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
velocities=V:V_step:LastVelocityToScan; %The row vector containing all the velocities to be scanned
NumOfScannedV=length(velocities);       %the number of the Scanned Velocities
angle=angle:angleStep:lastAngle;        %theAngleVectorContaining all the angles to be scanned
NumOfAnglesScanned=length(angle);       %theNumberOfScannedAngles
samplesInWindow=round(winSizeInSec/dt); %250 for half second if dt=0.002
N=samplesInWindow;M=Num_Rec;            %variables Names for easier use later
maxNumOfScannedt_0=(totalTime-( (Num_Rec-1)*dz*1/LastVelocityToScan ))*1/t_0Step; %corresponds to the totalNumberOfWindowsFittingInTimeDomain
CohArr=zeros(NumOfScannedV,floor(maxNumOfScannedt_0));
CohArrStack=zeros(NumOfAnglesScanned,NumOfScannedV); %the space for the stacking
REcc=1:Num_Rec-1:TotalNumOfReceivers-Num_Rec;% for the starting receivers needed later
%REcc is a vector steping every Num_Rec-1 from the first receiverTotheEnd
%REcc=50:58;                            %If you want to analyze by trace number just change here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IF YOU WANT TO SKIP TO A CERTAIN RECEIVER DEPTH, THEN ACTIVATE THE FOLLOWING LINE
REcc=[sum(depthRecord(1,:) <= 3036) 2]; %vector containing receiver numbers which are at the top of the distance window
%IF YOU WANT TO ANALYZE A CERTAIN DEPTH INTERVAL, THEN ACTIVATE THE FOLLOWING LINE
%REcc=[sum(depthRecord(1,:) <= 3000):sum(depthRecord(1,:) <= 3100)+1] 
%the last one wont be read  SELECT THE TOP RECEIVERS OF ANALYSIS WINDOW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D=zeros(Num_Rec,samplesInWindow,length(REcc)-1); %Creating the space for the D matrix
%checking that the selected receiver number is not too deep
for i=1:length(REcc)-1
    if REcc(i)+Num_Rec > TotalNumOfReceivers
        fprintf('ERROR \n');
        fprintf('The selected top of analysis window, ReceiverNumber=%d is too deep \n',REcc(i));
        fprintf('The max receiver should be: %d \n',TotalNumOfReceivers-Num_Rec);
        fprintf('In line 67 Change REcc to [%d : %d] \n',REcc(1),TotalNumOfReceivers-Num_Rec+1);
        return
    end
end
%checking that the velocity intervals are ok
if V >= LastVelocityToScan
    fprintf('ERROR \n');
    fprintf('The selected lastVelocity=%d is greater than the InitialVelocity=%d \n',LastVelocityToScan,V);
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FinalCoh=zeros(NumOfAnglesScanned,NumOfScannedV,length(REcc)-1); %creating the space for the finalCoherency as a function of V,alpha,z
initialAngle=angle(1); %saving this one for later print
initialV=V;%saving for later print
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%loop for filling TheTimeRecord space
for count=1:Num_Rec
     T(count,:)=dt:dt:dt*nt; %filling TheTimeRecord space
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%THE MAIN LOOP of the program that goes through depth
for mastercount=1:length(REcc)-1 %THE MAIN LOOP that goes through depth
    %selecting only the receivers we are interested in
 Vx=seisX(REcc(mastercount):REcc(mastercount)+Num_Rec-1,:); %horizontal1 direction
 Vz=seisZ(REcc(mastercount):REcc(mastercount)+Num_Rec-1,:); %horizontal2 direction
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% WIGGLE PLOTTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     if plotSeis == 1 %&& mastercount == floor(TotalNumOfReceivers/Num_Rec) %So that it only prints the last analysis window
     %loop for obtaining the value of the displacement for the vertical axis plot in matlab
             for count=1:Num_Rec
                 [m1(count),p]=max(abs(Vx(count,:))); %locating the maximum amplitude of Vz and saving it in vector m1 containing all max amplitudes of all traces, p is the position of the maxAmplitude in the Vx vector
                 [m3(count),p]=max(abs(Vz(count,:))); %locating the maximum amplitude of Vz and saving it in vector m3 containing all max amplitudes of all traces, p is the position in the Vx vector
             end

             [mx,~]=max(m1);[mz,p]=max(m3); %getting the maximum amplitude from all the vx and vz traces, p is the index in the m vector

             %loop for adding the displacement for each trace
             for count2=1:Num_Rec
                 for k=1:nt
                     Vx(count2,k)=Vx(count2,k)+count2*mx; %adding the maximum amplitude at each trace i*mx is the separation between traces
                     Vz(count2,k)=Vz(count2,k)+count2*mz;
                 end
             end

         figure;% Vx
         for i=1:Num_Rec
         plot(T(i,:),Vx(i,:),'-k','LineWidth',1.5);hold on;
         end
         %SAVING THE FIGURE VX
         set(gca,'YTickLabel','','YTick',[]);
         g=sprintf('analysis window of Vx for interval %d top trace is %d at %dm',mastercount,REcc(mastercount),depthRecord(REcc(mastercount)));
         title(g);
         xlabel('Time(s)');ylabel('Receiver');
         set(gca, 'YDir','reverse')
         filename=sprintf('analysis window of Vx for interval %d top trace is %d at %dm',mastercount,REcc(mastercount),depthRecord(REcc(mastercount)));
         saveas(gcf,filename,'png')
         hold off;set(gcf,'Color','w');clear filename;

         figure;% Vz
         for i=1:Num_Rec
         plot(T(i,:),Vz(i,:),'-k','LineWidth',1.5);hold on;
         end
         %SAVING THE FIGURE VZ
         g=sprintf('analysis window of Vz for interval %d top trace is %d at %dm',mastercount,REcc(mastercount),depthRecord(REcc(mastercount)));
         title(g);
         set(gca,'YTickLabel','','YTick',[]);
         xlabel('Time(s)');ylabel('Receiver');
         set(gca, 'YDir','reverse')
         filename=sprintf('analysis window of Vz for interval %d top trace is %d at %dm',mastercount,REcc(mastercount),depthRecord(REcc(mastercount)));
         saveas(gcf,filename,'png')
         hold off;set(gcf,'Color','w');clear filename;
         %Recovering the amplitudes
         load(filenameX);
         load(filenameZ);
         Vx=seisX(REcc(mastercount):REcc(mastercount)+Num_Rec-1,:); %horizontal1 direction
         Vz=seisZ(REcc(mastercount):REcc(mastercount)+Num_Rec-1,:); %horizontal2 direction
         close all;
           
     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
     %END OF WIGGLE PLOTTER     %after running that previous loop
                                %the True registered amplitude of the traces is TrueAmp=Vx(i,k)=Vx(1,k)-i*mx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 %printing some values to know the parameters
  if mastercount==1
  format shortg
  disp(fix(clock))
  fprintf('The file being used is %s \n',filenameX);
  fprintf('the time step is: %d \n', t_0Step);
  fprintf('the velocity step is: %d \n', V_step);
  fprintf('The number of receivers in the window is: %d \n', Num_Rec);
  fprintf('The length of the depth window is %d \n',(Num_Rec-1)*dz);
  fprintf('The initial scanning velocity was %d \n',V);
  fprintf('The initial scanning angle was %d \n',initialAngle);
  fprintf('The total number of receivers is: %d \n',TotalNumOfReceivers);
  fprintf('The window Size is %d \n', winSizeInSec);
  fprintf('If you are saving in pptx file check that powerpoint is closed');
  end
  %% STEP 1: rotate at current angle----------------------------------------------
 %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 for counterAlpha=1:NumOfAnglesScanned %THE ANGLE LOOP------------------------------------------------------------------
 %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
         alpha=angle(counterAlpha);
         seisRot = methodRotation(Vx,Vz,alpha); %Creating the rotated seismogram (see function for details)
 %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 for counterV = 1:NumOfScannedV %THE VELOCITY LOOP-------------------------------------------------------------------------------
 %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     %% STEP 2: picking a velocity
                        velocity=velocities(counterV); %the current scanning velocity
                         %VECTOR CONTAINING ALL THE DEPTHS OF CHOSEN RECEIVERS
                         delta_z= depthRecord(REcc(mastercount)):dz:depthRecord(REcc(mastercount))+(Num_Rec-1)*dz;
                         initialTime=delta_z(1)/velocity; %because we start by reading the first receiver check Dij eq in pevzner and notice how only t_0 is left
                         indexOfReferenceTime= sum(T(1,:) <= initialTime); %getting the index of the initial time t_0
                            
                         
                        %% STEP 3: picking a time and getting the coherency
                         %this loop scans at many t_0, at single v and then saves the
                         %results in CohArr
                         NumOfScannedt_0= (totalTime-( (Num_Rec-1)*dz*1/velocity ) - initialTime )*1/t_0Step; %the number of windows
 %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
                    for countert_0=1:floor(NumOfScannedt_0) %THE refTime LOOP how many refTimes Or windows you want to scan-------------------------------------------------------------------
 %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        %getting the coherency after the velocity loop
                        t_0=initialTime; %current reference time
                        t_0=t_0+t_0Step; %changing the reference time to be evaluated | t_0Step=winSizeInSec/2;
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        Coh = methodCoherency(t_0,delta_z,velocity,seisRot,Num_Rec,D,timeRecord,nt,N,M); %getting the coherency with this function
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                %fprintf('The Coherency at %d[m/s] alpha= %d and t_0= %d is %d \n', velocities(counterV),angle(counterAlpha),refTimes(countert_0),Coh); %---------------------------------%%TO CHECK THE RESULTS REDUCE THE NUMBER OF ANGLES            
                            CohArr(counterV,countert_0) = Coh; %saving the current coherency of the scanning t_0 at a fixed V  (it is also at a fixed alpha but this is 2D son not important)
                            %the row are for V the columns are refTime
                    end %OfRefTimeLoop
                 end%OfVelocityLoop
    %% loop for stacking at reference times (means adding the elements of a column to produce one value)
    for stackedV=1:counterV %is the same as NumOfScannedV
        CohArrStack(counterAlpha,stackedV)=sum(CohArr(stackedV,:)); 
        %each ROW ELEMENT of this array is the coherency of one ScannedVelocity stacked at multiple travel times
        %the total number of rows is the number of angles scanned
    end %of alphaVS Velocities matrix
end %end of alpha(azimuth) loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FinalCoh(1:NumOfAnglesScanned,1:NumOfScannedV,mastercount)=CohArrStack(:,:);%matrix that takes the angles,Velocities,depths in each dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTTING THE FIGURES
figH=figure(mastercount); %The figure handle
    %minx = round(min(FinalCoh(:,:,1)),1);
    %maxx = round(max(FinalCoh(:,:,1)),1);
    %levels =  minx:0.1:maxx;
    [X,Y] =meshgrid(angle,velocities);
    %contour(X,Y,FinalCoh(:,:,mastercount),'--')
    Z=FinalCoh(:,:,mastercount)';
    Z=smoothdata(Z,'gaussian',5); %you need the curveFittingToolBox
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
    title(sprintf('WinSize=%d Vcoherency VS alpha top depth of:%d and bottom of %d ',winSizeInSec,depthRecord(REcc(mastercount)),depthRecord(REcc(mastercount))+(Num_Rec-1)*dz));
    xlabel('Angle (deg)')
    xticks(0:10:180)
    ylabel('Apparent Velocity (m/s)')
    view(2)
    %% Saving the figure
    filename=strcat('Vfrom',num2str(V),'to',num2str(velocity),'Vstep',num2str(V_step),'AngStep',num2str(angleStep),'WinSize=',erase(num2str(winSizeInSec),'.'),'Num_Rec=',num2str((Num_Rec)),'DepthInterval=',num2str( depthRecord(REcc(mastercount)) ),'to',num2str( depthRecord(  REcc(mastercount) ) +(Num_Rec-1)*dz));
    filename=strcat(filenameX,filename);
    saveas(gcf,filename,'png')
    saveas(gcf,filename,'tif')

    %% Start new presentation
if mastercount == 1
    pptx    = exportToPPTX('', ...
    'Dimensions',[12 6], ...
    'Title','Example Presentation', ...
    'Author','MATLAB', ...
    'Subject','Automatically generated PPTX file', ...
    'Comments','This file has been automatically generated by exportToPPTX');
end

% Another way of setting some properties
pptx.title          = 'Demonstration Presentation';
pptx.author         = 'exportToPPTX Example';
pptx.subject        = 'Demonstration of various exportToPPTX commands';
pptx.description    = 'Description goes in here';

% Additionally background color for all slides can be set as follows:
% exportToPPTX('new','BackgroundColor',[0.5 0.5 0.5]);


%% Add some slides
%figH = figure('Renderer','zbuffer'); mesh(peaks); view(0,0);


    slideId = pptx.addSlide();
    fprintf('Added slide %d\n',slideId);
    pptx.addPicture(figH);
    pptx.addTextbox(sprintf('Slide Number %d',slideId));
    pptx.addNote(sprintf('Notes data: slide number %d',slideId));
    
    % Rotate mesh on each slide
    %view(18*islide,18*islide);

close(figH);


%% Check current presentation
%fprintf('Presentation size: %f x %f\n',pptx.dimensions);
%fprintf('Number of slides: %d\n',pptx.numSlides);


%% Save presentation  -- overwrite file if it already exists
% Filename automatically checked for proper extension
% Filename automatically checked for proper extension
[~,name,~]=fileparts(pwd); %getting the working folder name
newFile = pptx.save(strcat(name,'_numRec=',num2str(Num_Rec)));
%newFile     = pptx.save(strcat('results','numRec=',num2str(Num_Rec)));

end % end of master depth loop

clear pptx %close presentation
%% 
    disp('the velocity loop in the program ran this many times:');
    disp(counterV);
    disp('the angle loop in the program ran this many times:');
    disp(counterAlpha);
    disp('the referenceTime loop in the program ran this many times:');
    disp(countert_0);
    fprintf('The last velocity was: %d \n',velocities(counterV));
    fprintf('The last time is: %d \n', t_0+( delta_z(Num_Rec)-delta_z(1) )/ velocities(counterV));
    fprintf('The last angle was: %d \n',angle(counterAlpha));
    fprintf('The last reference time was: %d \n', t_0Step*NumOfScannedt_0 );
    
    figure
    
    plot(velocities,CohArrStack(1,:)) %the CohArrStack(1,:) corresponds to the velocity stacking at one angle angle(counterAlpha=1) for the last layer
    title(sprintf('Coherency at fixed [t_0 and alpha= %d ] C(t_0,v,alpha)',angle(1)));
    xlabel('velocities')
    xlim([velocities(1) velocities(end)])
    filename=sprintf('Coherency at fixed [t_0 and alpha= %d ] C(t_0,v,alpha)',angle(1));
    filename=strcat(filenameX,filename);
    saveas(gcf,filename,'png')
    %depth = z_0:dz*(Num_Rec-1):z_0+(dz-1)*TotalNumOfReceivers;
 %delete(gcp('nocreate'));%shutdown the parallel pool
 tEnd = toc(tStart);
 fprintf('the program took: %d [s] \n', tEnd);
 fprintf('the program took: %d minutes \n', tEnd/60);
%% END OF PROGRAM
 diary off
save FinalCoh.mat FinalCoh
save ElapsedTime.mat tEnd
