% computes SNR in MTL subfields for comparing EPI to mux data
%
% note: if using motion-corrected or smoothed data, will first need to
% combine 3D files into a 4D image: system('fslmerge -t rmux rmux0*')

addpath(genpath('/Users/vacarr/Desktop/MCI/mciExpt/scripts/kendrick/'));

subs={'mt01'};

fType = 'mux'; %epi or mux

mask_list = {'rcL-PRC.nii','rcR-PRC.nii','rcL-PHC.nii','rcR-PHC.nii',...
    'rcL-ERC.nii','rcR-ERC.nii','rcL-AntHipp.nii','rcR-AntHipp.nii',...
    'rcL-PostHipp.nii','rcR-PostHipp.nii'};
mask_names = {'LPRC','RPRC','LPHC','RPHC','LERC','RERC','LAHIPP','RAHIPP','LPHIPP','RPHIPP'};
fp=fullfile('/Users/vacarr/Desktop/mux_data/muxtest/');

for s=1:length(subs)
    
    sub=subs{s};
    par = mt_params_noknk(sub);
    
    subdir = [fp sub];
%     
%     % original data
%     if strcmp(fType,'epi')
%         func = {[subdir '/epi/functional/epi.nii']};
%     elseif strcmp(fType,'mux')
%         func = {[subdir '/mux/functional/mux.nii']};
%     end
    
    % motion corrected data
    if strcmp(fType,'epi')
        func = {[subdir '/epi/functional/repi.nii']};
    elseif strcmp(fType,'mux')
        func = {[subdir '/mux/functional/rmux.nii']};
    end
    
%     % motion corrected data + smoothed mux
%     if strcmp(fType,'epi')
%         func = {[subdir '/epi/functional/repi.nii']};
%     elseif strcmp(fType,'mux')
%         func = {[subdir '/mux/functional/srmux.nii']};
%     end
    
    run=load_untouch_nii(func{1});
    run=run.img;
    
    run=double(run);
    
    mean_run=mean(run(:,:,:,8:end),4);
    std_run=std(run(:,:,:,8:end),0,4);
    
    tsnr=mean_run./std_run;
    
    % write out volume & figures
    %     V=spm_vol(epi);
    %     V=V(1);
    %     V.fname=sprintf('/Users/anthonywagner/Desktop/pilotfm/SNR/SNR_%s.nii',sub);
    %     spm_write_vol(V,temporalsnr);
    %     imwrite(uint8(255*makeimagestack(tsnr,[0 20])),hot(256),sprintf('%stemporalsnr_%s.png',fp,sub));
    
    % agggregate ROI values
    if strcmp(fType,'epi')
        rp=fullfile(subdir,'/anatomy/rois_regEPI/');
    elseif strcmp(fType,'mux')
        rp=fullfile(subdir,'/anatomy/rois_regMux/');
    end
    
    for r=1:length(mask_list)
        
        roi=load_untouch_nii(fullfile(rp, mask_list{r}));
        roi=double(roi.img);
        
        val.(mask_names{r})(s) = mean(tsnr(roi==1));
        
    end
    
    
end