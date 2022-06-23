%% this script draws a quick blocky column with the depths of strata
% and after shows the way the BLOCKS OF depths that were analyazed
%THIS PROGRAM REQUIRES 'syms'
%syms package of matlab must be installed
close all; clear all; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FIRST part that needs some input from seismic 
%% INPUTS-----------------------------------------------------------------

filenameX='Well_x_stack_new.mat'; %filename handle
load(filenameX); %%Loading the data from seismic
dz=15;             %space between receivers at depth (z-axis)
z_0=1251;          %initialDepth of first receiver
Num_Rec=34;    % total number of receivers in window
save=0; % if you want to save EVERY analysis window as an image save=1
LastMeassuredDepth=3413; %last MD depth in meters from formation tops file
%%variables for pevzner method REQUIRED FOR THIS SCRIPT--------------------------------------
[TotalNumOfReceivers,nt] = size(seisX); %must be equal to seisZ
T=zeros(Num_Rec,nt);%space for the time array
TotalDepth=z_0+dz*TotalNumOfReceivers;          %Total depth of the survey Num_Rec is the number of traces
depthRecord=z_0:dz:TotalDepth-dz;   %The row vector  containing all the [z depth] points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Second part of the program that draws the column

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms x %using the syms package of matlab as exercise
f1= 0.00001*x;% so we draw the horizontal or vertical lines
%% INPUT HERE MUST BE A FILE WITH FORMATION TOPS ON THE FIRST COLUMN AND DEPTHS ON THE RIGHT COLUMN
fileName='FormationTops_Well.xlsx';
colors = {'white'};% { 'black';'red'; 'blue'; 'green'; 'yellow' };

%below the formation names
    formationNames={'Dammam';'Rus';'Umm Er Radhuma';'Simsima';'Fiqa';'Halul';'Laffan';'Mishrif';'Shilaif';'nahr-umr';'thamama1';'thamama2';'thamama3';'thamama4';'thamama5';'thamama6';'Hith';'ArabA';'ArabB';'ArabC';'ArabD';'Diyab';'end'};
%below the lithologies
    Litho={}; %still under progress so ignore
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data=xlsread(fileName); %reads only the numbers in the excel file
size(data);
depthToUnit=data(:,10); %the two is the second column representing meters
fprintf('The number of receivers in the window is: %d \n', Num_Rec);
  fprintf('The length of the depth window is %d \n',(Num_Rec-1)*dz);
%this condition below is in case there is no first zero value representing
%the beginning of the measures
if depthToUnit(1) > 0
depthToUnit=[0;depthToUnit];
fprintf('Added a zero at the beginning of MD to represent surface \n');
else
end
% appending adding a beginning to formation names
formationNames=['surface';formationNames];

horizontalLim=depthToUnit(end)/2; % the limit of the x axis
unitThickness=zeros(size(data,1),1); % the space of the thickness of units

%checking that the data has the last measured depth so it can plot all tops
depthToUnit=[depthToUnit;LastMeassuredDepth];


% loop for getting the unit thicknesses
for i=1:size(data,1)-1
    unitThickness_temp=depthToUnit(i+1)-depthToUnit(i);
    unitThickness(i)=unitThickness_temp;
end

%% section for drawing the column
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i=1:length(depthToUnit)-1 %you have to add a base
    figure(1)
    set(gcf, 'Position', get(0, 'Screensize')); %making it full screen
ezplot(f1+depthToUnit(i)) %plotting the horizontal lines representing each formation top
hold on
axis equal
ylim([-100 depthToUnit(end)]+100)%100s just to give some edges while plotting
xlim([0 horizontalLim])
title('Well formations column')
ylabel('Depth To Top Of Unit [m]')
set(gca, 'YDir','reverse') %plotting with the ydirection down


% selecting the lims of the area between depths
upperLim=depthToUnit(i);
lowerLim=depthToUnit(i+1);

%vectors containing coordinates clockwise dir from bottom to left corner
x = [-horizontalLim horizontalLim horizontalLim -horizontalLim];
y = [upperLim upperLim lowerLim lowerLim];


% this if else exists because the color labels may be less than the number
% of units
if i < size(colors,1)
patch(x,y,char(colors(i)))%this line draws the surface

else
    patch(x,y,char(colors(end)))
end
text(x(2),y(2),strcat( 'Top ',{' '},char( formationNames(i) ) )) %writing the formation names
%text(x(2)/2-100,y(2)+10,char(formationNames(i)))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%
%IF YOU WANT A SPECIFIC WINDOW UNLOCK THIS BLOCK

wantedDepth=2301; %check depth record so it can be the exact same depth
scanwin_x= [-horizontalLim horizontalLim horizontalLim -horizontalLim];
scanwin_y= [depthRecord(depthRecord==wantedDepth) depthRecord(depthRecord==wantedDepth) depthRecord(depthRecord==wantedDepth+(Num_Rec-1)*dz) depthRecord(depthRecord==wantedDepth+(Num_Rec-1)*dz)];
patch(scanwin_x,scanwin_y,'red')
alpha(0.3)
%}
%run til here (SELECT ALL TEXT ABOVE AND PRESS F9 )
% or press run section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% loop for drawing the intervals of depth analysis
figure(1)
pause(5)
for k=1:length(depthRecord)-Num_Rec
%drawing the depth interval analysis box
scanwin_x= [-horizontalLim horizontalLim horizontalLim -horizontalLim];
scanwin_y= [depthRecord(k) depthRecord(k) depthRecord(k+Num_Rec-1) depthRecord(k+Num_Rec-1)];
patch(scanwin_x,scanwin_y,'red') %drwaing the analysis window in red

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%% Saving the figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %filename=strcat(filenameX,filename);
    %if you want to save
    if save==1
    filenameFig=strcat('Well_DepthInterval=',num2str( depthRecord(k) ),'to',num2str( depthRecord(  k+Num_Rec-1 ) +(Num_Rec-1)*dz));
    saveas(gcf,filenameFig,'png') %gcf is get current figure
    else
        if k==1
            pause(3)
        else
        end
        pause(0.2)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%% third part putting everything into a pptx

figH=figure(1); %The figure handle this figure will be put in every slide
    %% Start new presentation
    pptx    = exportToPPTX('', ...
    'Dimensions',[12 6], ...
    'Title','Example Presentation', ...
    'Author','MATLAB', ...
    'Subject','Automatically generated PPTX file', ...
    'Comments','This file has been automatically generated by exportToPPTX');

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

%close(figH);


%% Check current presentation
%fprintf('Presentation size: %f x %f\n',pptx.dimensions);
%fprintf('Number of slides: %d\n',pptx.numSlides);


%% Save presentation  -- overwrite file if it already exists
% Filename automatically checked for proper extension
newFile     = pptx.save(strcat('ColumOfDepthWindow','numRec=',num2str(Num_Rec)));
end
close all;
clear pptx %close presentation
%}





