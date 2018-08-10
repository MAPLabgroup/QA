
function [] = tsnr_voxlev(subj,masktype)

% masks out spikes, computes voxel-wise tsnr, averages across mask, stores.
% assumes getslicextpsignal has already been run (which implies brain masks
% already made)

% subj: subject
% masktype: which brain mask to use: 'anti' (create liberal brain mask, use inverse), 'white' (white matter), or 'brain' (brain)

% dependencies: spm, par_savefields.m, matchfiles.m (optional to identify scans)

% karen larocque, may 2014

%% params

% directory & file paths

if strfind(subj,'PSER')

    resdir=fullfile('/Users/anthonywagner/Desktop/PSERscan',subj,'spikedetection'); % where your results will end up
    scans=matchfiles(['/Volumes/LaCie2/PSERraw/' subj '/*EPI_fMRI*/*.nii*'],'tr'); % cell array of paths to functional runs
    brainmask=fullfile(resdir,sprintf('tsmask_%s.nii',masktype));
    
elseif strfind(subj,'EB')

    resdir=fullfile('/Users/anthonywagner/Desktop/EBscan',subj,'spikedetection'); % where your results will end up
    scans=matchfiles(['/Volumes/LaCie2/EBraw/' subj '/raw/*EPI_fMRI*/*.nii*'],'tr'); % cell array of paths to 4D functional runs
    brainmask=fullfile(resdir,sprintf('tsmask_%s.nii',masktype));
    
else
   
    fprintf('cannot determine matching experiment directory for %s, exiting\n',subj);
    return
    
end

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
    V=V(:,:,:,ndrop+1:end);
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


