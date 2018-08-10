
function [] = TIB_tsnr_voxlev(subj,runsvec,masktype, exp_name)

% masks out spikes, computes voxel-wise tsnr, averages across mask, stores.
%
% NOTE: assumes getslicextpsignal has already been run (which implies brain masks
% already made)

% subj: subject
% masktype: which brain mask to use: 'anti' (create liberal brain mask, use inverse), 'white' (white matter), or 'brain' (brain)

% dependencies: spm, par_savefields.m, matchfiles.m (optional to identify scans)

% Written by Karen larocque, may 2014
%
% modified by TIB, Oct 2014. Adds "runsvec" input (a vector of run
% durations in TR) to account for bad or undesired scans at end of runs

%% params

% directory & file paths

%if strfind(subj,'PSER')

    %par.subdir=fullfile('/Volumes/LaCie/thackery_brown/', exp_name,'/fermi_nr_4D', subj);%fullfile('/Users/anthonywagner/Desktop/PSERscan',subj,'spikedetection'); % where your results will end up
    par.subdir = fullfile('/Users/thackery/Work/', exp_name, '/QA/fermi_nr_4D', subj);%fullfile('/Users/thackery/Work/Circmaze/bobspiketest/', subj)%fullfile('/mnt/wagner/thackery/', exp_name, '/subs_preprocessed/', subj);

    resdir=[par.subdir '/spikedetection'];
    %scans=matchfiles(['/Volumes/LaCie2/PSERraw/' subj '/*EPI_fMRI*/*.nii*'],'tr'); % cell array of paths to functional runs
    %brainmask=fullfile(resdir,sprintf('tsmask_%s.nii',masktype));
    par.funcdir =[par.subdir '/bolds/'] %subdirectory of subject folder that contains bold folders
    maskdir = fullfile('/mnt/wagner/thackery/', exp_name, '/subs_preprocessed_fermi_nospikerepair/', subj, '/Masks'); %if we are specifying our own mask(s)
    brainmask = fullfile(maskdir, 'rca1_trimmed.nii');
% elseif strfind(subj,'EB')
% 
%     resdir=fullfile('/Users/anthonywagner/Desktop/EBscan',subj,'spikedetection'); % where your results will end up
%     scans=matchfiles(['/Volumes/LaCie2/EBraw/' subj '/raw/*EPI_fMRI*/*.nii*'],'tr'); % cell array of paths to 4D functional runs
%     brainmask=fullfile(resdir,sprintf('tsmask_%s.nii',masktype));
%     
% else
%    
%     fprintf('cannot determine matching experiment directory for %s, exiting\n',subj);
%     return
%     
% end

runfolds = dir(fullfile(par.funcdir, 'run_*'));
for idxr = 1:length(runfolds)
    allrawfilenames{idxr,1} = dir(fullfile(par.funcdir, runfolds(idxr).name, '/*_1.nii*'));%'/run*_4d.nii'));
    for idxf = 1:length(allrawfilenames{idxr})
        allrawfilepaths{idxr,1}{idxf,1} = runfolds(idxr).name;
    end
end
allrawfilenames = vertcat(allrawfilenames{:});
allrawfilepaths = vertcat(allrawfilepaths{:});
for idx = 1:length(allrawfilenames)%-1%note, we are filling in the beta file names based on how many betas OF INTEREST we have (length(betaidx)). We don't care about the error reg betas for this analysis
    raw_filenames{idx,1} = [par.funcdir char(allrawfilepaths(idx)) '/' allrawfilenames(idx).name]; %create the analog to "raw_filenames.mat" - i.e. a list of all filenames including the path
end
for idx = 1:length(raw_filenames)
    if length(raw_filenames{idx,1}) == 77%83%longer names (runs 10+)
        raw_filenames{idx,2} = str2double(raw_filenames{idx,1}(length(raw_filenames{idx,1})-7:length(raw_filenames{idx,1})-6));
    else
        raw_filenames{idx,2} = str2double(raw_filenames{idx,1}(length(raw_filenames{idx,1})-6:length(raw_filenames{idx,1})-6));
    end
end
a = sortrows(raw_filenames, 2);
scans = a(:,1);


% drop how many scans
ndrop=6; % drop the first 'ndrop' volumes before calculating M/SD

% spike detection params
thresh=3;
iters=2;

%% load slice x timepoint signal, identify spike tps

tmp=load(fullfile(resdir,sprintf('slcmeans_%s.mat',masktype)));
slcmeans=tmp.slcmeans.globalcent; % use global centered data to detect spikes
clear('tmp');

%% identify spike tps, calculats tsnr, mean tsnr for each scan

Mh=spm_vol(brainmask);
M=spm_read_vols(Mh);

for s=1:length(scans)
    rundur = runsvec(s)%pull up current run's length, add dropped volumes count so we read the 4D files to the right length
    fprintf('Starting scan %d\n',s);
    
    % identify spikes
    slc=slcmeans{s}; 
    cleanslc=slc; % spikes removed from cleanslc; cleanslc used to calculate mean / sd
    spk=zeros(size(slc));

    for i=1:iters

        cleanslc(spk~=0)=NaN;

        mm=repmat(nanmean(cleanslc,2),1,size(cleanslc,2));
        sd=repmat(nanstd(cleanslc,[],2),1,size(cleanslc,2));

        zscored=((slc-mm)./sd);
        spk(abs(zscored)>=thresh)=zscored(abs(zscored)>=thresh);

    end

    spk=spk>0;
    
    % calculate tsnr
    scn=scans{s};

    if any(strfind(scans{s},'.gz'))
        system(sprintf('gunzip %s',scn));
        fprintf('unzipping scan %d\n',s);
        scn=strrep(scn,'.gz','');
    end

    Vh=spm_vol(scn);
    V=spm_read_vols(Vh); % inplane x inplane x slice x time

    assert(size(M,1)==size(V,1) & size(M,2)==size(M,2) & size(V,3)==size(V,3),'mismatch between brain mask & scan dimensions!')
    
    % drop leadin, take mean across inplane 
    V=V(:,:,:,ndrop+1:rundur);%end);bounding by number of good TRs in run (rundur)
    assert(size(slc,1)==size(V,3) & size(slc,2)==size(V,4),'mismatch between stored slice means & current scan dimensions!')
    
    % mask out spike timepoints
    for j=1:size(V,3)
        V(:,:,j,spk(j,:)~=0)=NaN; %mask out spike timepoints
    end
    
    % mask out voxels outside of the mask
    for j=1:size(V,4)
        tmp=V(:,:,:,j);
        tmp(M==0)=NaN; % ignore voxels that aren't in mask
        V(:,:,:,j)=tmp;
    end
    
    voxmean=nanmean(V,4);
    voxstd=nanstd(V,[],4);
    voxtsnr=voxmean./voxstd;
    voxtsnr(voxtsnr==Inf)=NaN;
    
    tsnr(s)=nanmean(voxtsnr(:)); %#ok<AGROW,NASGU>

    clear('V','Vh');

    fprintf('Done with scan %d\n',s);

end

save(fullfile(resdir,sprintf('tsnr_%s_t%d.mat',masktype,thresh)),'tsnr');


