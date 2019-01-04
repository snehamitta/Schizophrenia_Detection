function show_maps(D,cc)

switch length(size(D))
    case 4
        
        D = permute(D,[2 1 3 4]);
        t = squeeze(D(:,:,:,cc));
    case 3
        
        t = permute(D,[2 1 3]);
end

[val, idx] = max(abs(t(:)));
[y x z] = ind2sub(size(t),idx);


figure
subplot(1,3,1)  
I = transpose(squeeze(t(:,x,:)));
imagesc(I,max(abs(I(:)))*[-1 1]), axis xy equal tight
subplot(1,3,2)  
I = fliplr(transpose(squeeze(t(y,:,:))));
imagesc(I,max(abs(I(:)))*[-1 1]), axis xy equal tight
subplot(1,3,3)  
I = fliplr(squeeze(t(:,:,z)));
imagesc(I,max(abs(I(:)))*[-1 1]), axis xy equal tight