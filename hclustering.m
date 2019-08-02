

function CM = hclustering(pixel_vector, Xd)

[ylen, xlen] = size(Xd);

% feature vectors are divided into three categories by using FCM

options = [2.0; 100; 1e-5; 0];

fprintf('... ... 1st round clustering ... ...\n');
[center,U,obj_fcn] = fcm(pixel_vector,2, options);

maxU = max(U);
index{1} = find(U(1,:) == maxU);
index{2} = find(U(2,:) == maxU);  
if numel(index{1})<numel(index{2})
    ttr = numel(index{1})/(ylen*xlen)*1.25;
    ttl = numel(index{1})/(ylen*xlen)/1.10;
else
    ttr = numel(index{2})/(ylen*xlen)*1.25;
    ttl = numel(index{2})/(ylen*xlen)/1.10;
end

c_num = 5;
fprintf('... ... 2nd round clustering ... ...\n');
[center,U,obj_fcn] = fcm(pixel_vector,c_num, options);

Xdk =  zeros(ylen*xlen, 1);
CMk0 = zeros(ylen*xlen, 1);

Xdk = reshape(Xd, ylen*xlen, 1);

maxU = max(U);

for i = 1:c_num
    index{i} = find(U(i,:) == maxU);    
end

for i = 1:c_num
    idx_mean(i) = mean(Xdk(index{i}));
end

[idx_mean, idx] = sort(idx_mean);

for i = 1:c_num
    idx_num(i) = numel(index{idx(i)});
end

CMk0(index{idx(c_num)}) = 0.0;
c = idx_num(c_num);
mid_lab = 0;

for i = 1:c_num-1
    c = c+idx_num(c_num-i);
    if c / (ylen*xlen) < ttl
       CMk0(index{idx(c_num-i)}) = 0.0;
    elseif c / (ylen*xlen) >= ttl && c / (ylen*xlen) < ttr
        CMk0(index{idx(c_num-i)}) = 0.5;
        mid_lab = 1;
    else
        if mid_lab == 0
            CMk0(index{idx(c_num-i)}) = 0.5;
            mid_lab = 1;
        else
            CMk0(index{idx(c_num-i)}) = 1;
        end
    end
end

CM = reshape(CMk0, ylen, xlen);















