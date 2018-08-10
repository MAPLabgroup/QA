function [] = getslicextpsignal_batch%(subj)
% getslicextpsignal_batch
% batches getslicextpsignal to run on multiple subs using multiple masks
%throwaway = subj;
%addpath(genpath(['/Users/vacarr/Desktop/MCI/mciExpt/scripts/kendrick/']));
%Study name
S.exp_name = 'Circmaze'; %change this to flexibly redirect the script to different studies

%groups = {'mci','old'};

%dataFormats = {'older','newer'};

masklist={'anti','brain','white'};


subjlist={'cm100'};%'cm07', 'cm08', 'cm09', 'cm10', 'cm12', 'cm13',


for s=1:length(subjlist)
    for m=1:length(masklist)
        %getslicextpsignal(subType,dataType,subjlist{s},masklist{m});
        getslicextpsignal(subjlist{s},masklist{m}, S.exp_name);
    end
    
end

end

