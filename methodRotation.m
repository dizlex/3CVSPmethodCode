%{
this is the first function needed
the one that rotates the seismograms
this will only rotate the Vz first, then I will expand
BEFORE running the function do in command window:
seisX=load('inputX.txt'); seisY=load('inputY.txt');
%}

%{
the INPUTS are:
depthRecord
V_H1 the seismogram measurement in H1 direction
V_H2 the seismogram measuremnet in H2 direction
alpha the angle of azimuth
the OUTPUTS are: 
seisRot = pevznerRotation
%}

%coder.inline('never');
%{
 so that this function wont be inlined in the generated code
 meaning that a separate pevznerRotation.c will be created
%}

%good practice to call the function the same name as the file
function rotated = methodRotation(V_H1,V_H2,alpha)
    rotated = V_H1*cos(alpha*3.1415926535/180) + V_H2*sin(alpha*3.1415926535/180); %from deg to rad you do pi/180
end
%{
    seisXrot=V_H1*cos(alpha*3.1415926535/180);
    seisYrot=V_H2*sin(alpha*3.1415926535/180);
    seisRot = seisXrot + seisYrot;
 %}


