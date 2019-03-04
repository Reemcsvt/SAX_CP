function [sax_stringSD] = sax_demo_NPlot(data,nseg,alphabet_size, union_arr)
   
    data_len      = length(data);  
    %Normalization
    data = (data - mean(data))/std(data);
    %%%%%%%%%%%%%%%%%%%%SD%%%%%%%%%%%%%%%%%%%%%%%
    for i_s_ip=1:length(union_arr)
        clear sss;
        if i_s_ip == 1
           sss=data(1:union_arr(i_s_ip));
           PAA_SD(i_s_ip)=std(sss,1);            
        else
           sss=data((union_arr(i_s_ip-1)+1):union_arr(i_s_ip));
           PAA_SD(i_s_ip)=std(sss,1);
        end
    end
    if union_arr(i_s_ip)<length(data)
       sss=data(union_arr(i_s_ip)+1:length(data));
       PAA_SD(i_s_ip+1)=std(sss,1);
    end;
    %%%%%%%%%%%%%%%%%%Sum of(x-mean(x))^2%%%%%%%%%%%%%%%%%%%%
    for i_s_ip=1:length(union_arr)
        clear sss a_s2;
        if i_s_ip == 1
           sss= data(1:union_arr(i_s_ip));
           M_Seq(i_s_ip)=(mean(sss));
           a_s2=sss-M_Seq(i_s_ip); 
           PPA_s2(i_s_ip)=sum(a_s2.^2);
        else
           sss= data((union_arr(i_s_ip-1)+1):union_arr(i_s_ip));
           M_Seq(i_s_ip)=(mean(sss));
           a_s2=sss-M_Seq(i_s_ip); 
           PPA_s2(i_s_ip)=sum(a_s2.^2);
        end
    end
    clear sss a_s2;
    if union_arr(i_s_ip)<length(data)
        sss=data(union_arr(i_s_ip)+1:length(data));
        M_Seq(i_s_ip+1)=(mean(sss));
        a_s2=sss-M_Seq(i_s_ip+1); 
        PPA_s2(i_s_ip+1)=sum(a_s2.^2); 
    end;     
    %%%%%%%%%%%%%%%%%%%%%PAA_sign%%%%%%%%%%%%%%%%%%%%%%
    sidy=1; sid_r=1; IP_Pos=union_arr(sid_r);
    for sidx=1:data_len 
        if sidx>IP_Pos 
           sidy=sidy+1; 
           sid_r=sid_r+1;
           if sid_r>length(union_arr)
              IP_Pos=data_len;
           else
           IP_Pos=union_arr(sid_r);
           end;
        end    
            if data(sidx)-M_Seq(sidy)>=0
               PAA_sign(sidx)=1;
            else
               PAA_sign(sidx)=-1;
            end;        
    end;  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    str = timeseries2symbol_IP(data, data_len, nseg, alphabet_size, union_arr);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get the breakpoints
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
    
    symbols = {'1', '2', '3', '4', '5', '6', '7', '8', '9','10', '11', '12', '13', '14', '15', '16', '17', '18', '19','20'};
    sax_string = symbols(str);
    sax_stringSD.a= sax_string; %list
    sax_stringSD.s=PPA_s2;      %number
    sax_stringSD.n=PAA_sign;    %list
    sax_stringSD.sd=PAA_SD;

