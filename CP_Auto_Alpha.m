%------------------------------------------------------------------------------------%
%Copyright and terms of use (DO NOT REMOVE THE HEADER):
%  
% This file is part of SAX_CP.
% SAX_CP is a free project
%
% "A Novel Trend based SAX Reduction Technique for Time Series" 
% Authors: Hamdi Yahyaoui and Reem Al-Daihani.
%  
% SAX_CP can not be copied and/or distributed without the express
% permission of the authors
%
% Copyright (C) 2019 Written by Reem Aldaihani  All rights reserved.
%
%------------------------------------------------------------------------------------%
function alpha85App = CP_Auto_Alpha(nnseq)
global  TrainFile 
data1 = xlsread(TrainFile);  
dataT=data1(:,2:end);
%dataT
[size_s,Slength]=size(dataT(:,:));
for i=1:size_s
dataT(i,:) = zscore(dataT(i,:));
end
union_arr = findchangepts(dataT,'Statistic','rms','MaxNumChanges',nnseq-1);
%%%%%%%%%%%   Mean  %%%%%%%%%%%%
for TT=1:size_s
    data=dataT(TT,:);
    for i_s_ip=1:length(union_arr)
            clear sss;
            if i_s_ip == 1
               sss=data(1:union_arr(i_s_ip));
               PPA(i_s_ip)=mean(sss);            
            else
               sss=data((union_arr(i_s_ip-1)+1):union_arr(i_s_ip));
               PPA(i_s_ip)=mean(sss);
            end;
    end;
    if union_arr(i_s_ip)<length(data)
       sss=data(union_arr(i_s_ip)+1:length(data));
       PPA(i_s_ip+1)=mean(sss);
     end;
     Smean(TT,:)=PPA;
end;    
L_S=length(Smean(1,:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for seg_no=1:size_s
    for obs=1:L_S
        for alpha=2:20
            exists_a=Try_a((Smean(seg_no,obs)),alpha);
            if exists_a==1
               Res_85New(seg_no,obs)= alpha;
               break
            end
        end
    end
end  
alpha85App=ceil(max(median(Res_85New)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function exs=Try_a(segno,Ea)
switch Ea
        case 2, if segno >= -Inf && segno<=0
                exs=1; 
                else
                exs=0;     
                end    
        case 3, if segno >= -Inf && segno<=0.43
                exs=1; 
                else
                exs=0;     
                end
        case 4, if segno >= -Inf && segno<= 0.67  
                exs=1; 
                else
                exs=0;     
                end
         case 5, if segno >= -Inf && segno<=0.84
                exs=1; 
                else
                exs=0;     
                end                
        case 6, if segno >= -Inf && segno<=0.97  
                exs=1; 
                else
                exs=0;     
                end
        case 7, if segno >= -Inf && segno<=1.07  
                    exs=1; 
                else
                exs=0;     
                end
        case 8, if segno >= -Inf && segno<=1.15  
                    exs=1; 
                else
                exs=0;     
                end
        case 9, if segno >= -Inf && segno<=1.22  
                    exs=1; 
                else
                exs=0;     
                end
        case 10,  if segno >= -Inf && segno<=1.28  
                    exs=1; 
                else
                exs=0;     
                end
        case 11,  if segno >= -Inf && segno<=1.34  
                     exs=1; 
                  else
                  exs=0;     
                  end
        case 12,  if segno >= -Inf && segno<=1.38  
                    exs=1;
                  else
                  exs=0;     
                  end
        case 13,  if segno >= -Inf && segno<=1.43  
                    exs=1;
                  else
                  exs=0;     
                  end
        case 14,  if segno >= -Inf && segno<=1.47  
                    exs=1; 
                  else
                  exs=0;     
                  end
        case 15,  if segno >= -Inf && segno<=1.5  
                   exs=1; 
                  else
                  exs=0;     
                  end
        case 16,  if segno >= -Inf && segno<=1.53  
                    exs=1; 
                  else
                  exs=0;     
                  end
        case 17,  if segno >= -Inf && segno<=1.56  
                    exs=1; 
                  else
                  exs=0;     
                  end
        case 18,  if segno >= -Inf && segno<=1.59  
                    exs=1; 
                  else
                  exs=0;     
                  end
        case 19,  if segno >= -Inf && segno<=1.62  
                     exs=1; 
                  else
                  exs=0;     
                  end
        case 20,  if segno >= -Inf  && segno<=1.64  
                     exs=1; 
                  else
                  exs=1;     
                  end
    otherwise, exs=1;          
end;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cutlines = FBP(alphabet_size)
switch alphabet_size
        case 2, cutlines  = [-inf 0];
        case 3, cutlines  = [-inf -0.43 0.43];
        case 4, cutlines  = [-inf -0.67 0 0.67];
        case 5, cutlines  = [-inf -0.84 -0.25 0.25 0.84];
        case 6, cutlines  = [-inf -0.97 -0.43 0 0.43 0.97];
        case 7, cutlines  = [-inf -1.07 -0.57 -0.18 0.18 0.57 1.07];
        case 8, cutlines  = [-inf -1.15 -0.67 -0.32 0 0.32 0.67 1.15];
        case 9, cutlines  = [-inf -1.22 -0.76 -0.43 -0.14 0.14 0.43 0.76 1.22];
        case 10, cutlines = [-inf -1.28 -0.84 -0.52 -0.25 0. 0.25 0.52 0.84 1.28];
        case 11, cutlines = [-inf -1.34 -0.91 -0.6 -0.35 -0.11 0.11 0.35 0.6 0.91 1.34];
        case 12, cutlines = [-inf -1.38 -0.97 -0.67 -0.43 -0.21 0 0.21 0.43 0.67 0.97 1.38];
        case 13, cutlines = [-inf -1.43 -1.02 -0.74 -0.5 -0.29 -0.1 0.1 0.29 0.5 0.74 1.02 1.43];
        case 14, cutlines = [-inf -1.47 -1.07 -0.79 -0.57 -0.37 -0.18 0 0.18 0.37 0.57 0.79 1.07 1.47];
        case 15, cutlines = [-inf -1.5 -1.11 -0.84 -0.62 -0.43 -0.25 -0.08 0.08 0.25 0.43 0.62 0.84 1.11 1.5];
        case 16, cutlines = [-inf -1.53 -1.15 -0.89 -0.67 -0.49 -0.32 -0.16 0 0.16 0.32 0.49 0.67 0.89 1.15 1.53];
        case 17, cutlines = [-inf -1.56 -1.19 -0.93 -0.72 -0.54 -0.38 -0.22 -0.07 0.07 0.22 0.38 0.54 0.72 0.93 1.19 1.56];
        case 18, cutlines = [-inf -1.59 -1.22 -0.97 -0.76 -0.59 -0.43 -0.28 -0.14 0 0.14 0.28 0.43 0.59 0.76 0.97 1.22 1.59];
        case 19, cutlines = [-inf -1.62 -1.25 -1 -0.8 -0.63 -0.48 -0.34 -0.2 -0.07 0.07 0.2 0.34 0.48 0.63 0.8 1 1.25 1.62];
        case 20, cutlines = [-inf -1.64 -1.28 -1.04 -0.84 -0.67 -0.52 -0.39 -0.25 -0.13 0 0.13 0.25 0.39 0.52 0.67 0.84 1.04 1.28 1.64];
        otherwise disp('Error! alphabet_size is too big');           
end;
end