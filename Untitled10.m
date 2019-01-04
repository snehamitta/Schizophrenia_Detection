addpath('C:\Users\Admin\Desktop\spm8')
addpath('C:\Users\Admin\Desktop\libtiff_4')

datadir = 'C:\Users\Admin\Desktop';
fname = 'sub-01_task-letter0backtask_bold.nii';

Data2Read = fullfile(datadir,fname);

HeaderInfo = spm_vol(Data2Read);

Data = spm_read_vols(HeaderInfo);

signal_filt(Data,1,0.5,0.5,1);

mri1 = uint8(Data);
mri2 = squeeze(mri1);
mri3 = mri2(:,:,3);

csvwrite(test2,mri3);