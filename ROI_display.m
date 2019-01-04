AAL_map = spm_read_vols(spm_vol(fullfile(pwd,'aal_labels.nii')));
show_maps(AAL_map,1)
set(gca,'clim',[0 max(AAL_map(:))])
colormap colorcube


Df = spm_read_vols(spm_vol(fullfile(pwd,'roi_fmri.nii')));
fMRI_comp_ind = csvread('comp_fMRI.csv',1,0);
fMRI_ci = find(fMRI_comp_ind == 24); 
show_maps(Df,fMRI_ci)


Ds = spm_read_vols(spm_vol(fullfile(pwd,'roi_mri.nii')));
sMRI_comp_ind = csvread('comp_MRI.csv',1,0);
sMRI_ci = find(sMRI_comp_ind == 10); 
show_maps(Ds,sMRI_ci)


sz = size(Df);
Dftab = reshape(Df,prod(sz(1:3)),sz(4));
sz = size(Ds);
Dstab = reshape(Ds,prod(sz(1:3)),sz(4));


msk = sum(abs(Dftab),2) ~= 0;


cor_fs = corr(Dftab(msk,:),Dstab(msk,:)); 
figure
imagesc(cor_fs,[-.5 .5])
set(gca,'fontsize',5)
set(gca,'Xtick',1:length(sMRI_comp_ind),'Ytick',1:length(fMRI_comp_ind))
set(gca,'XtickLabel',sMRI_comp_ind,'YtickLabel',fMRI_comp_ind)
axis equal tight
colorbar