function []=showPatches(patches)
    figure
    it=1;
    p_size=sqrt(size(patches,1));
    for i=1:size(patches,2)
        subplot(p_size,p_size,it),imshow(reshape(patches(:,i),p_size,p_size),[min(patches(:,i)),max(patches(:,i))]);
        it=it+1;
    end
end