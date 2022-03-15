clear all;clc;
%read the data.txt generated from data.su using seismic unix
%load the ascii file
%INPUT: THE FILES TO BE READ in .txt
%PROCESS: removes all the headers and put every trace in a row of a final matrix.
%OUTPUT: THE seisx and seisZ variables
%GENERATE THE seisx and seisZ variables
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%INPUT
    fileName = 'sarb09_x_stack_new.txt'; %the file containg the x component H1
    fileName1 = 'sarb09_y_stack_new.txt'; %the file containg the Y component H2
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FID = fopen(fileName);
    data = textscan(FID,'%s');
    fclose(FID);
    stringData = string(data{:});

    %% ReadingTheDataX FIRST FILE
    %CreatingTheZeroMatrix for storing
    newStringData = stringData;
    for i=1:length(stringData)
        newStringData(i) = 0; 
    end
%% REMOVING ALL THE HEADERS
    %searches for all the e and keeps it, else gives 0
    pat = "e";
    TF = contains(stringData,pat);
    %make all the values that you do not care for 0
    for i=1:length(stringData)
        if TF(i)==1
        newStringData(i) = stringData(i);
        end
    end

    %searches for all equals and replaces them with zero
    pat = "=";
    TF = contains(newStringData,pat);
    %make all the values that you do not care for 0
    for i=1:length(stringData)
        if TF(i)==1
        newStringData(i) = 0;
        end
    end

    %If non zero the values that you are interested in get written to
    %finalStringData
     counter=1;
    for i=1:length(stringData)

        if newStringData(i) ~= "0"
        finalStringData(counter) = newStringData(i);
        counter=counter+1;
        end

    end
%% END
    %I have it all into one column of string
    finalStringData=finalStringData';
    %turning into rows
    %%%%%%%%%%%%% check how many number of samples are there ns=?
    %%% the example is Rows_of_2000 because there were ns=2000
    Rows_of_2000  = reshape(finalStringData, 6000, []).'; 
    clear finalStringData;clear pat; clear i;clear TF; clear counter; clear newStringData;
    [n1,n2]=size(Rows_of_2000);
    seisX=zeros(n1,n2);

    for i=1:n1
        for j=1:n2
        seisX(i,j)=str2num(Rows_of_2000(i,j));
        end
    end
    clear n1;clear n2; clear Rows_of_2000;
    
    
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    %% ReadingTheDataZ THE SECOND FILE
    
    FID = fopen(fileName1);
    data = textscan(FID,'%s');
    fclose(FID);
    stringData = string(data{:});

    
    %CreatingTheZeroMatrix
    newStringData = stringData;
    for i=1:length(stringData)
        newStringData(i) = 0; 
    end

%% removing the headers
    %searches for all the e and keeps it, else gives 0
    pat = "e";
    TF = contains(stringData,pat);
    %make all the values that you do not care for 0
    for i=1:length(stringData)
        if TF(i)==1
        newStringData(i) = stringData(i);
        end
    end

    %searches for all equals and replaces them with zero
    pat = "=";
    TF = contains(newStringData,pat);
    %make all the values that you do not care for 0
    for i=1:length(stringData)
        if TF(i)==1
        newStringData(i) = 0;
        end
    end

    %If non zero then replaces the values
     counter=1;
    for i=1:length(stringData)

        if newStringData(i) ~= "0"
        finalStringData(counter) = newStringData(i);
        counter=counter+1;
        end

    end

    %% reshaping the string of numbers only into rows
    %I have it all into one column of string
    finalStringData=finalStringData';
    %turning into rows
    Rows_of_2000  = reshape(finalStringData, 6000, []).';
    clear finalStringData;clear pat; clear i;clear TF; clear counter; clear newStringData;
    [n1,n2]=size(Rows_of_2000);
    seisY=zeros(n1,n2);

    for i=1:n1
        for j=1:n2
        seisY(i,j)=str2num(Rows_of_2000(i,j));
        end
    end
    clear n1;clear n2; clear Rows_of_2000; clear stringData;

save seisX.mat seisX
save seisY.mat seisY




