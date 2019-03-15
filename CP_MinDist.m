% Copyright and terms of use (DO NOT REMOVE):
% The code is made freely available for non-commercial uses only, provided that the copyright 
% header in each file not be removed, and suitable citation(s) (see below) be made for papers 
% published based on the code.
%
% The code is not optimized for speed, and we are not responsible for any errors that might
% occur in the code.
%
% The copyright of the code is retained by the authors.  By downloading/using this code you
% agree to all the terms stated above.
%
%   [1]"A Novel Trend based SAX Reduction Technique for Time Series" 
%       Authors: Hamdi Yahyaoui and Reem Al-Daihani, 2019
%   
%   [2] Lin, J., Keogh, E., Lonardi, S. & Chiu, B. 
%   "A Symbolic Representation of Time Series, with Implications for Streaming Algorithms." 
%   In proceedings of the 8th ACM SIGMOD Workshop on Research Issues in Data Mining and 
%   Knowledge Discovery. San Diego, CA. June 13, 2003. 
%
%   [3] Lin, J., Keogh, E., Patel, P. & Lonardi, S. 
%   "Finding Motifs in Time Series". In proceedings of the 2nd Workshop on Temporal Data Mining, 
%   at the 8th ACM SIGKDD International Conference on Knowledge Discovery and Data Mining. 
%   Edmonton, Alberta, Canada. July 23-26, 2002
%
%
% This function computes the minimum (lower-bounding) distance between two strings.  The strings
% should have equal length.
%   Input:
%       str1: first string
%       str2: second string
%       alphabet_size: alphabet size used to construct the strings
%       compression_ratio: original_data_len / symbolic_len
%   Output:
%       dist: lower-bounding distance
%
%   usage: dist = min_dist(str1, str2, alphabet_size, compression_ratio)
%
% This distance measure is not the best measure to use for comparing strings,
% if you are NOT going to follow up with access to the original data. This is
% because it cannot discriminate between two strings that differ only in the
% ith place, by consecutive symbols.  For example the min_dist between 'abba'
% and  'abbb' is zero.
% However, in practice, the min_dist function works very well for
% classification and clustering, even when you do not follow up with access to
% the original data. See [1].
%
%
%
% Copyright (c) 2003, Eamonn Keogh, Jessica Lin, Stefano Lonardi, Pranav Patel, Li Wei.  All rights reserved.
%------------------------------------------------------------------------------------%  
%------------------------------------------------------------------------------------%  
% This file is part of SAX_CP.
% SAX_CP is a free project
%
% "A Novel Trend based SAX Reduction Technique for Time Series" 
% Authors: Hamdi Yahyaoui and Reem Al-Daihani.
%  
% Copyright (C) 2019 Written by Reem Aldaihani All rights reserved.
%------------------------------------------------------------------------------------%
function dist_arr = CP_MinDist(str11, str_list1)
data_length=length(str11);
global SDT str111 dataO IPsax alphabet_size nnseq No_SeqOT No_SeqOS TrainFile TestFile
IPsax_n=nnseq-1;
if length(IPsax) ~= IPsax_n
    for isaxi=1:10
    [SSSAAA,err] = findchangepts(str_list1,'Statistic','rms','MaxNumChanges',IPsax_n);
        if length(SSSAAA) == IPsax_n
           IPsax= SSSAAA;
           isaxi=10;
        end;
    end;
end
win_size(1)=IPsax(1);
for i=2:length(IPsax)
    win_size(i)= abs(IPsax(i)-IPsax(i-1));
end
win_size(nnseq)=abs(IPsax(length(IPsax))-data_length); 

if exist('SAXT85.xlsx', 'file')==0 
    dataTTT = xlsread(TrainFile); 
    dataTT=dataTTT(:,2:end);
    for SEQO=1:No_SeqOT
        X=CP_SAX_TimeSeries_Computation(dataTT(SEQO,:),nnseq,alphabet_size,IPsax);
        AT(SEQO,:)=X.a;  AAT(SEQO,:)=X.s;  AAAT(SEQO,:)=X.n;  AAAAT(SEQO,:)=X.sd;
    end;
    xlswrite('SAXT85.xlsx',AT);   xlswrite('SEQT85.xlsx',AAT);  xlswrite('BETAT85.xlsx',AAAT);  xlswrite('SDT85.xlsx',AAAAT);
    SDT.a= xlsread('SAXT85.xlsx');
    SDT.s= xlsread('SEQT85.xlsx');
    SDT.n= xlsread('BETAT85.xlsx');
    SDT.sd= xlsread('SDT85.xlsx');
else
end;

if exist('SAXS85.xlsx', 'file')==0
    dataOO = xlsread(TestFile); 
    dataO=dataOO(:,2:end);
    for SEQO=1:No_SeqOS
    Y=CP_SAX_TimeSeries_Computation(dataO(SEQO,:),nnseq,alphabet_size,IPsax);
    AS(SEQO,:)=Y.a;    AAS(SEQO,:)=Y.s;   AAAS(SEQO,:)=Y.n;   AAAAS(SEQO,:)=Y.sd;
    end
    xlswrite('SAXS85.xlsx',AS);  xlswrite('SEQS85.xlsx',AAS); xlswrite('BETAS85.xlsx',AAAS); xlswrite('SDS85.xlsx',AAAAS);
    str111.a= xlsread('SAXS85.xlsx');
    str111.s= xlsread('SEQS85.xlsx');
    str111.n= xlsread('BETAS85.xlsx');
    str111.sd= xlsread('SDS85.xlsx');  
end
for SEQO=1:No_SeqOS
    if dataO(SEQO,:) == str11
     str1 = str111.a(SEQO,:);
     strs=  str111.s(SEQO,:);
     strn=  str111.n(SEQO,:);
     strsd= str111.sd(SEQO,:);
    end
end
[s_row,s_col]=size(str_list1);
% for each sequence 
for sax_i=1:s_row
      str_list(sax_i,:)=SDT.a(sax_i,:);
      str_S(sax_i,:)=SDT.s(sax_i,:);
      str_N(sax_i,:)=SDT.n(sax_i,:);
      str_SD(sax_i,:)=SDT.sd(sax_i,:);   
    if (length(str1) ~= length(str_list(sax_i,:)))
        display('error: the strings must have equal length!');
        return;
    end    
    if (any(str1 > alphabet_size) | any(str_list(sax_i,:) > alphabet_size))
        display('error: some symbol(s) in the string(s) exceed(s) the alphabet size!');
        return;
    end   
    dist_matrix = build_dist_table(alphabet_size);  
    dist = 0.0; 
    distt =  diag(dist_matrix(str1,str_list(sax_i,:)));
    for ws=1:nnseq
        dist=dist+(win_size(ws)*distt(ws));
    end
    
    dist2=0.0; 
    IPP=1;  IP_Pos=IPsax(IPP);  sqn=0; sum_n=0.0; nnn=1;
    for bb=1:data_length 
        if bb>IP_Pos  
           sbb(nnn)=1+(double(sum_n)/double(sqn));
           IPP=IPP+1;
           sum_n=0; sqn=0; nnn=nnn+1; 
           if IPP>length(IPsax)
              IP_Pos=data_length;
           else
           IP_Pos=IPsax(IPP);
           end;
        end;
        sum_n=sum_n+(strn(bb)* str_N(sax_i,bb));
        sqn=sqn+1;
    end
    sbb(nnn)=1+(double(sum_n)/double(sqn));
    
    sumBeta=0.0;
    for sbeta=1:nnseq
        sumBeta=sumBeta+((2-sbb(sbeta))*(sqrt(strs(sbeta)*str_S(sax_i,sbeta))));
    end;

    SD_dist=0.0;
    for fd=1:nnseq
        SD_dist = SD_dist + (win_size(fd)*((strsd(fd))-str_SD(sax_i,fd))^2);
    end   
    dist=dist+SD_dist;
    dist2=dist+sumBeta;
    dist_arr(sax_i,1)= dist2;  
end;
%------------------------------------------------------------------------------------------------------
%                               build_dist_table(alphabet_size)
%------------------------------------------------------------------------------------------------------

function dist_matrix = build_dist_table(alphabet_size)
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

    dist_matrix=zeros(alphabet_size,alphabet_size);
    for i = 1 : alphabet_size
        for j = i+2 : alphabet_size            
            dist_matrix(i,j)=(cutlines(i)-cutlines(j-1))^2;           
            dist_matrix(j,i)= dist_matrix(i,j);
        end;
    end;  