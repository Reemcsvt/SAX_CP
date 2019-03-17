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
% Copyright (C) 2019 Written by Reem Aldaihani reema@vt.edu  All rights reserved.
%
%------------------------------Global variables--------------------------------------%
global alphabet_size nnseq No_SeqOT No_SeqOS TrainFile TestFile
%------------------------------------Input-------------------------------------------%
TrainFile='TRAIN85';
TestFile='TEST85';
nnumbre_of_Segs=4;
%-----------------------------------------------------------------------------------%
data1 = xlsread(TrainFile); 
data2 = xlsread(TestFile); 
dataT=data1(:,2:end);   
dataS=data2(:,2:end);
class=data1(:,1); 
truelabels=data2(:,1);
No_SeqOT=length(data1);
No_SeqOS=length(data2);
%-----------------------------------------------------------------------------------%
disp('Running...');
for  nnseq= nnumbre_of_Segs
    if check_ip(dataT,nnseq-1) == 1
        afaf=CP_Auto_Alpha(nnseq);
        alphabet_size= afaf;
        classout = knnclassify_Mod(dataS, dataT, class,1,2); 
        CP = classperf(truelabels, classout);
        fprintf('Error rate  %f\n',CP.LastErrorRate);
    end
end 
delete('SAXT85.xlsx'); delete('SEQT85.xlsx'); delete('BETAT85.xlsx'); delete('SDT85.xlsx');
delete('SAXS85.xlsx'); delete('SEQS85.xlsx'); delete('BETAS85.xlsx'); delete('SDS85.xlsx');
%-------------------------------------------------------------------------------------%
%-----------------------------Local Function: check_ip--------------------------------%
%-------------------------------------------------------------------------------------%
function ip_exists= check_ip(data_ip,ip_no)
    Flag_ip=0;
    for isaxiii=1:10
        [res_ip,erra] = findchangepts(data_ip,'Statistic','rms','MaxNumChanges',ip_no);
            if length(res_ip) == ip_no
                ip_exists=1;
                isaxiii=10;
                Flag_ip=1;
            end
    end;
            if Flag_ip~=1    
                ip_exists=0;
            end;       
end