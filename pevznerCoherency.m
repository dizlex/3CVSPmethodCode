%this returns the Coherency at a single reference time but with multiple
%velocities that were scanned
function Coh = pevznerCoherency(t_0,delta_z,velocity,seisRot,Num_Rec,D,timeRecord,nt,N,M)
    Num_interior_sumOfAllTracesAtFixSample=zeros(1,nt); %preallocating    
    %FIRST a loop that goes through every trace and picks the initial time at each trace!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          for i=1:Num_Rec %so we go through every trace
              
               initialTime=t_0+(delta_z(i)-delta_z(1))/velocity; %first point in time domain at current scanning trace i
               index=sum(timeRecord(1,:) <= initialTime); %getting the index of the first point at the ith trace
               
              
               %fprintf('The initial time is %d \n',timeRecord(index));  
          %SECOND a loop that writes D_ij | that is to say %writing all the j values at the ith trace 
                   for j=0:N-1 %so we go every sample check the slides for the -1
                       %uncomment for debugging error caveat mistake
                       %{
                       if j+index-N/2 > length(timeRecord) %why would this case happen?
                           disp('the pointer is larger')
                           disp('Is z_0 in meters?')
                           disp(j+index-N/2)
                           
                       end
                       if j+index-N/2 <= 0
                           disp(j+index-N/2)
                           disp(j+1)
                       end
                       %}
                   D(i,j+1)=seisRot(i,j+ index -round(N/2)); %writing the values of the seisRot In the matrix
                   end
           
           end
          %THIRD calculating coherency at current V // NOTE THAT alpha and t_0 remain unchanged
                  %coherency is a function C(t_0,V,a), remember it returns something lower than one
                  %fprintf('The values of the seismogram are being read at this reference time: %d \n',T(indexOfReferenceTime));
                  %CALCULATING THE DENOMINATOR
                              %inside the parenthesis says: IN TIME the 4th power of the sum of all the traces at a fixed time
                                                          % IN INDEX the 4th power of the sum of all the i at a fixed j
                              %outside the parenthesis says: sum the inside of the parenthesis for all the times
       
                              %loop for outside the parenthesis check the slides
                              for j=1:N %going from the first sample in the window until the N sample in the window
                                  %inside the parenthesis
                                  %% Sum of the (M traces=:) at this time index and saving the results at this row vector
                                  Num_interior_sumOfAllTracesAtFixSample(1,j) = sum( D(:,j) ); %this vector at each entry has the sum of all traces at a fixed time k
                              end
                          
                              Num_interior_sumOfAllTracesAtFixSample = Num_interior_sumOfAllTracesAtFixSample.^4; % rising every element of the row vector containing the sum of all the traces at a common time to the 4th power
                              Num = sum ( Num_interior_sumOfAllTracesAtFixSample(1,:) ); %sum of all the elements of the vector
                  %the denominator of coherency
                      %it says: multiply M times; the sum of every element of the matrix squared
                              %M = Num_Rec; %setting M as in the paper,M is the number of traces in the window
                              squareTheElements = D.^2; %squaring all the elements of the matrix
                              %total=0;
           %loop for summing all the elements of D
                              %for quefijado=1:N*M %from 1 unitl the number of elements in the matrix
                              %    total= total + squareTheElements(quefijado); %going through every element of the matrix in horizontal
                              %end
                              sumTheElements = sum(squareTheElements,'all'); %summing all the elements of the matrix
                  Den = M * sumTheElements; % multiplying by M
                  
                  
                  Coh= Num / Den;
                  %dbstop if naninf
                 %fprintf('The Coherency at %d[m/s] alpha= %d and t_0= %d is %d \n', velocities(counterV),angle(counterAlpha),refTimes(countert_0),Coh); %---------------------------------%%TO CHECK THE RESULTS REDUCE THE NUMBER OF ANGLES
                          
              %CohArr(countert_0,counterV) = Coh; %saving the current coherency of the scanning V at a fixed t_o  (it is also at a fixed alpha but this is 2D son not important)
              %the row are for refTimes the columns are for V
                              %c is the velocity counter
                          
          
end
      