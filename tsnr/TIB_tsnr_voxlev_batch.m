function [] = TIB_tsnr_voxlev_batch
% batches getslicextpsignal to run on multiple subs using multiple masks
% throwaway = subj;
% addpath(genpath(['/Users/vacarr/Desktop/MCI/mciExpt/scripts/kendrick/']));
% Study name

S.exp_name = 'Circmaze'; %change this to flexibly redirect the script to different studies

masklist={'rca1_trimmed'};%{'anti','brain','white'};

subjlist={'cm12', 'cm16', 'cm17', 'cm21', 'cm22', 'cm23', 'cm24', 'cm25'};%'cm07', 'cm08', 'cm09', 'cm10', 'cm12', 'cm13','cm12'};%,   
runsvecs={[195	177	168	161	161	171	156	183	156	160	165	172	152	161	175	161], [174	167	166	163	162	164	161	165	159	162	161	161	160	161	160	160], [165	165	171	166	169	155	156	161	155	158	159	151	163	160], [158	170	163	162	160	157	159	154	161	156	157	159	162	155	153	155], [168	169	172	163	187	166	162	171	155	169	159	161	169	163], [157	162	161	159	159	159	155	164	160	164	160	161	153	162	165	159], [166	163	174	168	164	171	159	159	162	166	165	161	160	157	158	157], [161	162	151	153	168	160	162	154	156	166	156	161	152	156	162	150]};%
%subjlist={'cm12', 'cm16', 'cm17', 'cm21', 'cm22', 'cm24'};%'cm07', 'cm08', 'cm09', 'cm10', 'cm12', 'cm13','cm12'};%,   
%runsvecs={[195	177	168	161	161	171	156	183	156	160	165	172	152	161	175	161], [174	167	166	163	162	164	161	165	159	162	161	161	160	161	160	160], [165	165	171	166	169	155	156	161	155	158	159	151	163	160], [158	170	163	162	160	157	159	154	161	156	157	159	162	155	153	155], [168	169	172	163	187	166	162	171	155	169	159	161	169	163], [166	163	174	168	164	171	159	159	162	166	165	161	160	157	158	157]};%

for s=1:length(subjlist)
    for m=1:length(masklist)
        %getslicextpsignal(subType,dataType,subjlist{s},masklist{m});
        TIB_tsnr_voxlev(subjlist{s}, runsvecs{s}, masklist{m}, S.exp_name);
    end
    
end

end

