clear all;clc;
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%INPUT: THE FILES TO BE READ in .txt sorted in increasing order 
%PROCESS: remove the selected header and put every trace in a row of a matrix.
%the sorted means the minimum value is in the first row, and the maximum
%value is in the last row
%stack at same keyword value available
%OUTPUT: THE sarb-09 variables
%GENERATE THE sarb-09 variables
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%INPUT CHANGE HERE
    fileName = 'shah_x.txt'; %the file containg the x component H1
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    %fileName= 'testForStacking.txt';
    FID = fopen(fileName);
    data = textscan(FID,'%s');
    fclose(FID);
    stringData = string(data{:});
    %finding the maximum integer so you can get the number of samples
%-------------------------------------------------------------------------------
    %USR INPUT
    ns=2500; %number of elements in each matrix row
    %searches for all the gelev and keeps it, else gives 0
    keyword = "scalel"; %if you want to remove all the headers just put  keyword ="=";
%-----------------------------------------------------------------------------------
    %CreatingTheZeroMatrix
    newStringData = stringData;
    for i=1:length(stringData)
        newStringData(i) = 0; 
    end

    
    TF = contains(stringData,keyword);

    %make all the values that you do not care for 0
    for i=1:length(stringData)
        if TF(i)==1
        newStringData(i) = stringData(i); %copying it to the new matrix
        end
    end

    %% searches for all the e+ or e- and keeps it
    %this are the values of the amplitudes
    pat= "e+";
    TF = contains(stringData,pat);
    %make all the values that you do not care for logical 0
    for i=1:length(stringData)
        if TF(i)==1
        newStringData(i) = stringData(i); %copying it to the new matrix
        end
    end
    pat= "e-";
    TF = contains(stringData,pat);
    %make all the values that you do not care for logical 0
    for i=1:length(stringData)
        if TF(i)==1
        newStringData(i) = stringData(i); %copying it to the new matrix
        end
    end
    %after this part the amplitude numbers were picked

     %% If non zero the values that you are interested in get written to
    %finalStringData
     counter=1;
    for i=1:length(stringData)

        if newStringData(i) ~= "0"
        %this is into a single column
        finalStringData(counter)      = newStringData(i); %this is your data scrapped
        counter=counter+1;
        end
    end
    %The prvious loop got the data with the selected keyword in the the row
    %header

    %% this is the data in matrix form, where each row is a trace
    %I have it all into one column of string
    finalStringData=finalStringData';
    %turning into rows
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    MatrixWithARowForEachTrace  = reshape(finalStringData, ns+1, []).';
    %EX ns+1=2501 = numberOfSamples + numberOfKeywords

    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    clear finalStringData;clear pat; clear i;clear TF; clear counter; clear newStringData;
    [n1,n2]=size(MatrixWithARowForEachTrace);
%% now you extract the number associated with the keyword so the matrix has the number associated to keyword
 
 % this loop removes the keyword, leaving only the value at the first row
 for i=1:n1
     MatrixWithARowForEachTrace(i,1)=str2num(extractAfter(MatrixWithARowForEachTrace(i,1),"="));
 end
    seisX=zeros(n1,n2);
    % this loop creates the numerical matrix that can be handled by matlab
    for i=1:n1
        for j=1:n2
        seisX(i,j)=str2num(MatrixWithARowForEachTrace(i,j));
        end
    end

    %save keyword as a text file
    fid = fopen(strcat(keyword,'.txt'),'w');
    fprintf(fid,'%s\n',MatrixWithARowForEachTrace(:,1));
    fclose(fid);

    fid = fopen('CorrectGelevComparison.txt','w');
    fprintf(fid,'%s\n',mod);
    fclose(fid);

    stem(x,CorrectGelev)
    axis tight
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RUN UNTIL HERE TO CHECK THAT THE KEYWORD HAS BEEN EXTRACTED
% SELECT EVERYTHING ABOVE AND PRESS F9 TO RUN UNTIL HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% now if the keyword difference between each  trace is between 26 and 51
    %%you sum the traces into a single trace
   fid = fopen('GelevDifference.txt','w');
  for i=1:length(CorrectGelev)-1%because you are taking pairs to move the pointer
      if abs( CorrectGelev(i,1)-CorrectGelev(i+1,1) ) <= 55 && abs( CorrectGelev(i,1)-CorrectGelev(i+1,1) ) >= 20 %if the stacking is between more than two traces the program does not work
          fprintf(fid,'%0.5f \n', abs( CorrectGelev(i,1)-CorrectGelev(i+1,1) ));
      else
      end
  end
  fclose(fid);
  %{
          for j=1:n2-1
          temp=seisX(i,j+1)+seisX(i+1,j+1);
          seisX(i,j+1)=temp; %the traces are merged into one
          end
          seisX(i+1,2:n2)= 0; %setting to zero the second\repeated row
          i+1;
      end
      end
      %}
  %% removing the zero rows
  % To remove all zeros rows from seisX
  seisX=seisX(any(seisX(:,2:n2),2),:);
  %size(seisX,2);
  %saving the header
  gelev=seisX(:,1);
  %removing the header column
  seisX(:,1)=[];
  %%%%%%%
  % visual inspection traces to be removed, because they are anomalous
  %%%%%%%
  anom_tracx=[6 11 19 46 76 132];
  %for that removes anomalous tracks
  for i=1:length(anom_tracx)
      seisX(anom_tracx(i)-i-1,:)=[]; 
  end

  seisY=seisX;
  save sarb09_x_stack_new.mat seisX

    





