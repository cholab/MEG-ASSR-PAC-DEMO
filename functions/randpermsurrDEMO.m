function [surrdata,surrtemplate] = randpermsurrDEMO(data,surrs,tr,method,wind)
xV = zeros(1,length(wind)*tr);
% go through filtered data and shuffle
surrdata = zeros(size(data,1),length(wind),tr,surrs);
temp = data(:,wind,:);
[xx,yy,zz] = size(temp);
temp = temp(:,:);
% get samples
randtime = randsample(round(length(xV)*.8),surrs)+round(length(xV)*.1);
parfor ii = 1:surrs
    temp2 = [temp(:,randtime(ii):end) temp(:,1:randtime(ii)-1)];
    temp2 = reshape(temp2,xx,yy,zz);
    if method == 1
        temp2 = abs(hilbert(permute(temp2,[2,1,3]))).^2;
        temp2 = permute(temp2,[2,1,3]);
    elseif method == 2
        temp2 = angle(hilbert(permute(temp2,[2,1,3])));
        temp2 = permute(temp2,[2,1,3]);
    end
    surrdata(:,:,:,ii) = temp2;
    %     clear temp*
end
surrtemplate = zeros(size(data,1),size(data,2), size(data,3),size(surrdata,4));
surrtemplate(:,wind,:,:) = surrdata;
end