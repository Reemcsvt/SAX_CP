global alphabet_size nnseq
data1 = xlsread('TRAIN85');    
dataT=data1(:,2:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%
data2 = xlsread('TEST85');     
dataS=data2(:,2:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%
nnumbre_of_Seqs=4;
%%%%%%%%%%%%%%%%%%%%%%%%%%
class= xlsread('TRAIN85','A:A');
i_a=1; i_b=1;
clear IPsax
for  nnseq= nnumbre_of_Seqs
    i_b=1;
    if check_ip(dataT,nnseq-1) == 1
        afaf=Alpha85NewApp_Aa85(nnseq);
        for alphabet_size= afaf
        classout = knnclassify_Mod(dataS, dataT, class,1,2); 
        truelabels= xlsread('TEST85','A:A');
        CP = classperf(truelabels, classout);
        a_result_all85(i_a,1)= CP.LastErrorRate;
        a_result_all85(i_a,2)=nnseq;
        a_result_all85(i_a,3)=alphabet_size;
        delete('SAXT85.xlsx'); delete('SEQT85.xlsx'); delete('BETAT85.xlsx'); delete('SDT85.xlsx');
        delete('SAXS85.xlsx'); delete('SEQS85.xlsx'); delete('BETAS85.xlsx'); delete('SDS85.xlsx');
        i_b=i_b+1;
        end
        clear IPsax
    end
    i_a=i_a+1;
end
%a_result_all85
xlswrite('ResultWow.xlsx',a_result_all85); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ip_exists= check_ip(data_ip,ip_no)
    tttt=0;
    for isaxiii=1:10
        [res_ip,erra] = findchangepts(data_ip,'Statistic','rms','MaxNumChanges',ip_no);
            if length(res_ip) == ip_no
                ip_exists=1;
                isaxiii=10;
                tttt=1;
            end
    end;
            if tttt~=1    
                ip_exists=0;
            end;       
    clear res_ip
end