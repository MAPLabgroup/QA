
function [] = tsnr_resamp_agg()
% aggregate tsnr across subjects
% karen larocque, may 2014, modified june 2014

%% general params

% experiment & subjects
%exp=[]; % assigned later based on subject
sublist={'PSER01','PSER02','PSER03','PSER04','PSER05','PSER06','PSER07','PSER08','PSER09','PSER10','PSER11','PSER12','PSER13','PSER14','PSER15','EB01','EB02','EB03','EB04','EB05','EB06'};

% paths
outputpath=fullfile('/Users/anthonywagner/Dropbox/spikeqa/'); % this is where your textfiles will print
respath=fullfile('/Users/anthonywagner/Desktop'); % beginning of path to where your slicemean.mat files are stored -- *will need to make final edit in loop below*

% what masks to use
masklist={'anti','white','brain'};

% what threshold
threshlist=[4,6];

%% prep text files

fid=fopen(fullfile(outputpath,'tsnr_group.csv'),'w');
fprintf(fid,'subj,mask,thresh,run,field,sample,value\n');


%% get data & output tsnr for each mask type & subject

for sbj=1:length(sublist)
    
    subj=sublist{sbj};
    if strfind(subj,'PSER'), exp = 'PSERscan';
    elseif strfind(subj,'EB'), exp = 'EBscan';
    else
        fprintf('Unidentified experiment for %s, skipping',subj);
        continue
    end
    resdir=fullfile(respath,exp,subj,'spikedetection'); % NEED TO CHANGE THIS TO APPROPRIATE DIRECTORY FOR YOUR EXPERIMENT
    fprintf('Outputting data for subject %s\n',subj);
    
    for msk=1:length(masklist)
        
        masktype=masklist{msk};
        
        for thresh=threshlist

            d=load(fullfile(resdir,sprintf('tsnr_%s_t%d.mat',masktype,thresh)));
            
            fields=fieldnames(d);
            
            for f = 1:length(fields)
                
                tmp=d.(fields{f});
                
                for b=1:size(tmp,1)
                    
                    for r=1:size(tmp,2)

                        fprintf(fid,'%s,%s,%d,%d,%s,%d,%f\n',subj,masktype,thresh,b,fields{f},r,tmp(b,r));
                        
                    end

                end
            end
        
        end
        
    end
end