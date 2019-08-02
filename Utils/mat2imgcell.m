function Img = mat2imgcell(D,height,width,ImgFormat)

N = size(D,2);
Img = cell(N,1);
if strcmp(ImgFormat,'gray')
    for i = 1:N
        Img{i} = reshape(D(:,i),height,width);
    end
elseif strcmp(ImgFormat,'color')
    for i = 1:N
        Img{i} = reshape(D(:,i),height,width,3);
    end
end


