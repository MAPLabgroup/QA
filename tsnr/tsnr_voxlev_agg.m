
function [] = tsnr_voxlev_agg()
% aggregate tsnr across subjects
% karen larocque, may 2014

%% general params

% experiment & subjects

exp='PSERscan';
switch exp
    case {'PSERscan'}
        sublist={'PSER01','PSER02','PSER03','PSER04','PSER05','PSER06','PSER07','PSER08','PSER09','PSER10','PSER11','PSER12','PSER13','PSER14','PSER15'};
    case {'EBscan'}
        sublist={'EB01','EB02','EB03','EB04','EB05','EB06'};
end

% paths
outputpath=fullfile('/Users/anthonywagner/Desktop/',exp); % this is where your textfiles will print
respath=fullfile('/Users/anthonywagner/Desktop',exp); % beginning of path to where your slicemean.mat files are stored -- *will need to make final edit in loop below*

% what masks to use
masklist={'anti','white','brain'};

% what threshold
thresh=3;

%% prep text files

fid=fopen(fullfile(outputpath,sprintf('tsnr_%s.csv',exp)),'w');
fprintf(fid,'subj,mask,run,tsnr\n');


%% get data & output tsnr for each mask type & subject

for sbj=1:length(sublist)
    
    subj=sublist{sbj};
    resdir=fullfile(respath,subj,'spikedetection'); % NEED TO CHANGE THIS TO APPROPRIATE DIRECTORY FOR YOUR EXPERIMENT
    fprintf('Outputting data for subject %s\n',subj);
    
    for msk=1:length(masklist)
        
        masktype=masklist{msk};

        load(fullfile(resdir,sprintf('tsnr_%s_t%d.mat',masktype,thresh)));
        
        for b=1:length(tsnr)
           
            fprintf(fid,'%s,%s,%d,%f\n',subj,masktype,b,tsnr(b));
            
        end
        
    end
end