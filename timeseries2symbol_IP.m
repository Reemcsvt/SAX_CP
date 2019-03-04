function [symbolic_data, pointers] =  sax_modified(data, N, n, alphabet_size, union_arr)
NR_opt = 2;
%win_size = floor(N/n);                              % win_size is the number of data points on the raw time series that will be mapped to a single symbol
pointers         = [];                                                  % Initialize pointers,
symbolic_data = zeros(1,n);                                             % Initialize symbolic_data with a void string, it will be removed later.
all_string = zeros(length(data)-N+1,n);
% Scan accross the time series extract sub sequences, and converting them to strings.
for i = 1 : length(data) - (N -1)                                       
    if mod(i, 1000) == 0
        disp(num2str(i));
    end
    sub_section = data(i:i + N -1);            

    for i_s_ip=1:length(union_arr)
        clear sss;
        if i_s_ip == 1
           sss=data(1:union_arr(i_s_ip));
           PPA(i_s_ip)=mean(sss);            
        else
           sss=data((union_arr(i_s_ip-1)+1):union_arr(i_s_ip));
           PPA(i_s_ip)=mean(sss);
        end
    end
    if union_arr(i_s_ip)<length(data)
       sss=data(union_arr(i_s_ip)+1:length(data));
       PPA(i_s_ip+1)=mean(sss);
    end;

    current_string = map_to_string(PPA,alphabet_size);          % Convert the PAA to a string. 

    if NR_opt == 2     
        if ~all(current_string == symbolic_data(end,:))             % If the string differs from its leftmost neighbor...
            symbolic_data    = [symbolic_data; current_string];     % ... add it to the set...
            pointers         = [pointers ; i];                      % ... and add a new pointer.
        end;
    end;    
    
end;

% Delete the first element, it was just used to initialize the data structure
symbolic_data(1,:) = [];                                               

%--------------------------------------------------------------------------------------------------------------------------------------------------------
%----------------Local Functions----------------------Local Functions----------------Local Functions----------------------Local Functions----------------
%--------------------------------------------------------------------------------------------------------------------------------------------------------

function string = map_to_string(PAA,alphabet_size)

string = zeros(1,length(PAA));

switch alphabet_size
        case 2, cut_points  = [-inf 0];
        case 3, cut_points  = [-inf -0.43 0.43];
        case 4, cut_points  = [-inf -0.67 0 0.67];
        case 5, cut_points  = [-inf -0.84 -0.25 0.25 0.84];
        case 6, cut_points  = [-inf -0.97 -0.43 0 0.43 0.97];
        case 7, cut_points  = [-inf -1.07 -0.57 -0.18 0.18 0.57 1.07];
        case 8, cut_points  = [-inf -1.15 -0.67 -0.32 0 0.32 0.67 1.15];
        case 9, cut_points  = [-inf -1.22 -0.76 -0.43 -0.14 0.14 0.43 0.76 1.22];
        case 10, cut_points = [-inf -1.28 -0.84 -0.52 -0.25 0. 0.25 0.52 0.84 1.28];
        case 11, cut_points = [-inf -1.34 -0.91 -0.6 -0.35 -0.11 0.11 0.35 0.6 0.91 1.34];
        case 12, cut_points = [-inf -1.38 -0.97 -0.67 -0.43 -0.21 0 0.21 0.43 0.67 0.97 1.38];
        case 13, cut_points = [-inf -1.43 -1.02 -0.74 -0.5 -0.29 -0.1 0.1 0.29 0.5 0.74 1.02 1.43];
        case 14, cut_points = [-inf -1.47 -1.07 -0.79 -0.57 -0.37 -0.18 0 0.18 0.37 0.57 0.79 1.07 1.47];
        case 15, cut_points = [-inf -1.5 -1.11 -0.84 -0.62 -0.43 -0.25 -0.08 0.08 0.25 0.43 0.62 0.84 1.11 1.5];
        case 16, cut_points = [-inf -1.53 -1.15 -0.89 -0.67 -0.49 -0.32 -0.16 0 0.16 0.32 0.49 0.67 0.89 1.15 1.53];
        case 17, cut_points = [-inf -1.56 -1.19 -0.93 -0.72 -0.54 -0.38 -0.22 -0.07 0.07 0.22 0.38 0.54 0.72 0.93 1.19 1.56];
        case 18, cut_points = [-inf -1.59 -1.22 -0.97 -0.76 -0.59 -0.43 -0.28 -0.14 0 0.14 0.28 0.43 0.59 0.76 0.97 1.22 1.59];
        case 19, cut_points = [-inf -1.62 -1.25 -1 -0.8 -0.63 -0.48 -0.34 -0.2 -0.07 0.07 0.2 0.34 0.48 0.63 0.8 1 1.25 1.62];
        case 20, cut_points = [-inf -1.64 -1.28 -1.04 -0.84 -0.67 -0.52 -0.39 -0.25 -0.13 0 0.13 0.25 0.39 0.52 0.67 0.84 1.04 1.28 1.64];
        otherwise disp('Error! alphabet_size is too big');           
end;
        
for i = 1 : length(PAA)    
    string(i) = sum( (cut_points <= PAA(i)), 2 );         % order is now: a = 1, b = 2, c = 3..
end; 