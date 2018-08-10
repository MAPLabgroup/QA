
function [] = TIB_getslicextpsignal(subj, runsvec, masktype, exp_name)%getslicextpsignal(subType,dataType,subj,masktype)

% store mean signal per slice at each time point, masked by various brain
% masks

% subj: subject
% masktype: which brain mask to use: 'anti' (create liberal brain mask, use inverse), 'white' (white matter), or 'brain' (brain)

% dependencies: spm, par_savefields.m, matchfiles.m (optional to identify scans)

%% params

% directory & file paths
%par.subdir = fullfile('/Volumes/LaCie/thackery_brown/', exp_name, '/fermi_nr_4D', subj);%fullfile('/Users/thackery/Work/Circmaze/bobspiketest/', subj)%fullfile('/mnt/wagner/thackery/', exp_name, '/subs_preprocessed/', subj);
par.subdir = fullfile('/Users/thackery/Work/', exp_name, '/QA/fermi_nr_4D', subj);%fullfile('/Users/thackery/Work/Circmaze/bobspiketest/', subj)%fullfile('/mnt/wagner/thackery/', exp_name, '/subs_preprocessed/', subj);

resdir=[par.subdir '/spikedetection']; %subdirectory of subject folder that will contain output of the scripts for that sub
par.funcdir =[par.subdir '/bolds/'] %subdirectory of subject folder that contains bold folders
maskdir = fullfile('/mnt/wagner/thackery/', exp_name, '/subs_preprocessed_fermi_nospikerepair/', subj, '/Masks'); %if we are specifying our own mask(s)

%resdir=fullfile('/Users/vacarr/Desktop/fMRI_MCI',subType,subj,'spikedetection'); % directory where your results will end up
%if strcmp(dataType,'newer')
%    scans=matchfiles(['/Volumes/LaCie6/MCI_raw/' subType '/' subj '/raw/niftis/*EPI_fMRI*/*_*.nii*'],'tr'); % cell array of paths to functional runs
%elseif strcmp(dataType,'older')
%    scans=matchfiles(['/Volumes/LaCie6/MCI_raw/' subType '/' subj '/raw/niftis/*EPI_fMRI*.nii*'],'tr'); % cell array of paths to functional runs
%end

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

%hack to get rid of weird duplicate lines with ._run in them...
% for l = 2:2:length(scans)
% test{l,1} = scans{l}
% end
% scans = test(~cellfun('isempty',test))

% make our own brainmasks? mean inplane scan -- used to generate relevant brain masks
if strcmp(masktype, {'anti','brain','white'})
    inplane=[par.funcdir 'run_01/meanurun1_0006.nii']; % mean inplane scan -- used to generate relevant brain masks
    %brainmask=fullfile(resdir,sprintf('tsmask_%s.nii',masktype)); % if your
    %brain mask already exists, where it lives. otherwise use [];
    brainmask=[];
else
    brainmask = fullfile(maskdir, 'rca3_trimmed.nii');
end

if~exist(resdir,'dir')
    system(sprintf('mkdir %s',resdir));
end

% drop how many scans
ndrop=6; % drop the first 'ndrop' volumes before calculating M/SD - if using raw bolds this might be first 6

%% create brain masks

if ~exist(brainmask,'file')
    
    fprintf('No existing brain mask found; creating brain mask using %s\n',inplane);
    
    if any(strfind(inplane,'.gz'))
        system(sprintf('gunzip %s',inplane));
        fprintf('unzipping inplane\n');
        inplane=strrep(inplane,'.gz','');
    end
    
    opts.biasfwhm=60;
    segopts = struct('biascor',1,'GM',[0 0 1],'WM',[0 0 1],'CSF',[0 0 0],'cleanup',0);
    
    res=spm_preproc(inplane,opts);
    [sn,~]=spm_prep2sn(res);
    [pth,nam,ext]  = spm_fileparts(inplane);
    par_savefields(fullfile(pth,[nam '_seg_sn.mat']),sn);
    spm_preproc_write(sn, segopts);
    
    grsegs(1,:) = fullfile(pth, sprintf('c1%s%s',nam,ext));
    grsegs(2,:) = fullfile(pth, sprintf('c2%s%s',nam,ext));
    
    maskimg = fullfile(pth,'mask.nii');
    smaskimg = fullfile(pth,'smask.nii');
    tsmaskimg = fullfile(pth,'tsmask.nii');
    
    switch masktype
        
        case {'brain'}
            spm_imcalc_ui(grsegs,maskimg,'i1+i2');
            spm_smooth(maskimg,smaskimg,[20 20 20]);
            spm_imcalc_ui(smaskimg,tsmaskimg,'i1>.2');
        case {'anti'}
            spm_imcalc_ui(grsegs,maskimg,'i1+i2');
            spm_smooth(maskimg,smaskimg,[35 35 35]);
            spm_imcalc_ui(smaskimg,tsmaskimg,'i1<.01');
        case {'white'}
            spm_imcalc_ui(sprintf('%s/c2%s%s',pth,nam,ext),tsmaskimg,'i1>.95')
        otherwise
            fprintf('Brain mask does not already exist and specified type is not recognized.\n');
            return
    end
    
    % save PDF of mask
    tmskh=spm_vol(sprintf('%s/tsmask.nii',pth));
    tmsk=spm_read_vols(tmskh);
    imwrite(uint8(255*makeimagestack(tmsk)),hot(256),sprintf('%s/tsmask_%s.png',resdir,masktype));
    clear('tmskh','tmsk');
    
    % cleanup
    system(sprintf('rm %s/c1%s%s',pth,nam,ext));
    system(sprintf('rm %s/c2%s%s',pth,nam,ext));
    system(sprintf('rm %s/m%s%s',pth,nam,ext));
    system(sprintf('rm %s/%s_seg_sn.mat',pth,nam));
    system(sprintf('rm %s/mask.nii',pth));
    system(sprintf('mv %s/smask.nii %s',pth,resdir));
    system(sprintf('mv %s/tsmask.nii %s/tsmask_%s.nii',pth,resdir,masktype));
    brainmask=fullfile(resdir,sprintf('tsmask_%s.nii',masktype));
    
    fprintf('Brain mask %s created\n',masktype);
    
end
%% save slice means at each timepoint for each scan

Mh=spm_vol(brainmask);
M=spm_read_vols(Mh);

for s=1:length(scans)
    
    rundur = runsvec(s);%+ndrop;%pull up current run's length, add dropped volumes count so we read the 4D files to the right length
    fprintf('Starting scan %d\n',s);
    
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
    Vc=V;
    
    % mask voxels outside of brain & take out global sig
    % loop for now :( fix later
    for tp=1:size(V,4)
        tmp=V(:,:,:,tp);
        tmp(M==0)=NaN; % ignore voxels that aren't in mask
        V(:,:,:,tp)=tmp;
        tmp_vec = tmp(:);
        tmp_vec = tmp_vec(~isnan(tmp_vec));
        tmp=tmp-mean(tmp_vec); %take out global signal
        Vc(:,:,:,tp)=tmp;
        
    end
    
    slcmeans.raw{s}=squeeze(nanmean(nanmean(V,1),2));
    slcmeans.globalcent{s}=squeeze(nanmean(nanmean(Vc,1),2));
    
    clear('V','Vh','Vc');
    
    fprintf('Done with scan %d\n',s);
    
end

save(fullfile(resdir,sprintf('slcmeans_%s.mat',masktype)),'slcmeans');


