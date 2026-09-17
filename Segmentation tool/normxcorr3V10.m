function C = normxcorr3V10(T, A, shape)
% C = normxcorr3(TEMPLATE, IMAGE, SHAPE)
%
%       TEMPLATE - type double, ndims==3, size <= size of image
%       IMAGE    - type double, ndims==3
%       SHAPE    - one of: 'valid', 'same', 'full'. same as conv2 shape parameter
%                  'full' by default
%
%       C        - values in [-1,1]. size depends on SHAPE
%
% the syntax of this function is identical to Matlab's
% normxcorr2, except that it's been extended to 3D matrices,
% and, the SHAPE parameter has been introduced as a convenience
%
% the SHAPE parameter has the same effect as it does for the CONVN function.
% see the documentation for CONVN for a more detailed explanation
%
% caveat emptor: this function does not perform the argument checking that
% normxcorr2 does. for example, it doesn't ensure that std(T(:))~=0
%
% daniel eaton, 2005, danieljameseaton@gmail.com




% 3D Gaussian filter reconstruction
% filter_size = [5, 5, 5]; stdd = 1;
% step_size = 100; pixel_size = 30;
% filter_size_in = (filter_size-1)/2;  %The default filter size is 2*ceil(2*std)+1.
% [xx,yy,zz] = meshgrid(-filter_size_in(2):filter_size_in(2),-filter_size_in(1):filter_size_in(1),-filter_size_in(3):filter_size_in(3));
% aspect_ratio = step_size/pixel_size;
% filter_Gau = exp(-(xx.^2+yy.^2+(aspect_ratio*zz).^2)/(2*stdd^2));
% filter_Gau(filter_Gau < eps*max(filter_Gau(:))) = 0;
% filter_Gau = filter_Gau/sum(filter_Gau(:));
%
%
%
% Gx_oprt = [-1; 0; 1];
% Gy_oprt = [-1; 0; 1];
% Gz_oprt = [-1; 0; 1];


%median -> Gauss -> Prewitt
img_medianA = zeros(size(A));
for i = 1:size(A,3)
    img_medianA(:,:,i) = medfilt2(A(:,:,i));
end
% img_filter = imfilter(img_median,filter_Gau);
% A = double(div_img(img_filter,Gx_oprt,Gy_oprt,Gz_oprt));
A = img_medianA;


if nargin<3
    shape = 'full';
end

if ndims(A)~=3 || ndims(T)~=3
    error('A and T must be 3 dimensional matrices');
end

szT = size(T);
szA = size(A);

if any(szT>szA)
    error('template must be smaller than image');
end

pSzT = prod(szT);

% make the running-sum/integral-images of A and A^2, which are
% used to speed up the computation of the NCC denominator
intImgA = integralImage(A,szT);
intImgA2 = integralImage(A.*A,szT);

szOut = size(intImgA);

% compute the numerator of the NCC
% emulate 3D correlation by rotating templates dimensions
% in 3D frequency-domain correlation is MUCH faster than the spatial-domain
% variety
rotT = flipdim(flipdim(flipdim(T,1),2),3); % this is rot90 in 3d
fftRotT = fftn(rotT,szOut);
fftA = fftn(A,szOut);
corrTA = real(ifftn(fftA.*fftRotT));
num = (corrTA - intImgA*sum(T(:))/pSzT ) / (pSzT-1);

% compute the denominator of the NCC
denomA = sqrt( ( intImgA2 - (intImgA.^2)/pSzT ) / (pSzT-1) );
denomT = std(T(:));
denom = denomT*denomA;

% compute the NCC
s = warning('off', 'MATLAB:divideByZero');
C = num ./ denom;
s = warning('on', 'MATLAB:divideByZero');

% replace the NaN (if any) with 0's
zeroInd = find(denomA==0);
C(zeroInd) = 0;

switch( lower(shape) )
    case 'full'
    case 'same'
        szTp = fix((szT-1)/2);
        C = C( szTp(1)+1:szTp(1)+szA(1), szTp(2)+1:szTp(2)+szA(2), szTp(3)+1:szTp(3)+szA(3) );
    case 'valid'
        C = C(szT(1):end-szT(1)+1,szT(1):end-szT(2)+1,szT(3):end-szT(3)+1);
    otherwise
        error(sprintf('unknown SHAPE %s, assuming FULL by default', shape));
end

    function integralImageA = integralImage(A,szT)
        % this is adapted from Matlab's normxcorr2
        
        szA = size(A);
        
        B = zeros( szA+2*szT-1);
        B( szT(1)+1:szT(1)+szA(1), szT(2)+1:szT(2)+szA(2), szT(3)+1:szT(3)+szA(3) ) = A;
        
        s = cumsum(B,1);
        c = s(1+szT(1):end,:,:)-s(1:end-szT(1),:,:);
        s = cumsum(c,2);
        c = s(:,1+szT(2):end,:)-s(:,1:end-szT(2),:);
        s = cumsum(c,3);
        integralImageA = s(:,:,1+szT(3):end)-s(:,:,1:end-szT(3));
    end

CrossSize = size(C);
FinalX = []; FinalY = []; FinalZ = [];
for z = 1: CrossSize(3)
    data = C(:,:,z);
    [y,x] = find(data >= 0.85*max(data(:)));
    
    for i= 1:length(x)-1
        if (x(i+1)-x(i)<7)&&(y(i+1)-y(i)<7)
            x(i) = 0;
            y(i) = 0;
        end
    end
    
    x(x==0)=[];
    y(y==0)=[];
    z = ones(numel(x),1)*z;
    
    FinalX = [x;FinalX];
    FinalY = [y;FinalY];
    FinalZ = [z;FinalZ];
    
end

Final = [FinalZ, FinalX, FinalY];
Final = sortrows(Final);

figure(5)
imshow(max(A,[],3),[0 200]);
axis equal;
colormap('gray');


for i= 1:length(FinalZ)-1
    if (Final(i+1,1)-Final(i,1)<7)&&(Final(i+1,2)-Final(i,2)<7)&&(Final(i+1,3)-Final(i,3)<7)
        if C(Final(i),Final(i),Final(i))< C(Final(i+1),Final(i+1),Final(i+1))
            Final(i,1) = 0; Final(i,2) = 0; Final(i,3) = 0;
        else
            Final(i+1,1) = Final(i,1); Final(i+1,2) = Final(i,2); Final(i+1,3) = Final(i,3);
            Final(i,1) = 0; Final(i,2) = 0; Final(i,3) = 0;
        end
    end
end

zFinal = Final(:,1);
xFinal = Final(:,2);
yFinal = Final(:,3);


xFinal(xFinal==0)=[];
yFinal(yFinal==0)=[];
zFinal(zFinal==0)=[];

Final2 = [xFinal, yFinal, zFinal];
[~,ia,~] = unique(Final2(:,1:2),'rows');
Final2 = Final2(ia,:);


Final2 = sortrows(Final2);


for i= 1:length(Final2)-1
    if (Final2(i+1,1)-Final2(i,1)<7)&&(Final2(i+1,2)-Final2(i,2)<7)&&(Final2(i+1,3)-Final2(i,3)<7)
        Final2(i,1) = 0;
        Final2(i,2) = 0;
        Final2(i,3) = 0;
    end
end

xFinal = Final2(:,1); yFinal = Final2(:,2); zFinal = Final2(:,3);

xFinal(xFinal==0)=[]; yFinal(yFinal==0)=[]; zFinal(zFinal==0)=[];


Final2 = [xFinal, yFinal, zFinal];
Final2 = sortrows(Final2,2);

for i= 1:length(Final2)-1
    if (Final2(i+1,1)-Final2(i,1)<7)&&(Final2(i+1,2)-Final2(i,2)<7)&&(Final2(i+1,3)-Final2(i,3)<7)
        Final2(i,1) = 0;
        Final2(i,2) = 0;
        Final2(i,3) = 0;
    end
end

xFinal = Final2(:,1); yFinal = Final2(:,2); zFinal = Final2(:,3);

xFinal(xFinal==0)=[]; yFinal(yFinal==0)=[]; zFinal(zFinal==0)=[];



%% removing edge detections
for j = 1:length(Final2)
    if ((Final2(j,1)< size(T,2) || Final2(j,1)> size(A,1)-size(T,2))...
            || (Final2(j,2)<size(T,1) || Final2(j,2)> size(A,2)-size(T,1))...
            || (Final2(j,3)<size(T,3) || Final2(j,3)>size(A,3)-size(T,3)))
        
        Final2(j,:) = 0;
        
    end
end
Final2 = Final2(any(Final2,2),:);



%FineLocalization
numofbeads = size(Final2,1);
fineLoc = zeros(numofbeads,3);

for n = 1:numofbeads
    bead = A(Final2(n,1)-size(T,1)/1.1 : Final2(n,1)+size(T,1)/1.1,...
        Final2(n,2)-size(T,2)/1.1 : Final2(n,2)+size(T,2)/1.1,...
        Final2(n,3)-size(T,3)/1.1 : Final2(n,3)+size(T,3)/1.1);
    bead = bead/mean(bead(:));

    CM = centerOfMass(bead);

        
    newrow = Final2(n,1)-(CM(2)-size(T,1)/1.1);
    newcol = Final2(n,2)-(CM(1)-size(T,2)/1.1);
    newAxi = Final2(n,3)-(CM(3)-size(T,3)/1.1);
    
    fineLoc(n,1) = newrow;
    fineLoc(n,2) = newcol;
    fineLoc(n,3) = newAxi;
  
end

ynew = fineLoc(:,2) - size(T,2)/2;
xnew = fineLoc(:,1) - size(T,1)/2;

for i = 1:length(xnew)
    rectangle('Position',[xnew(i), ynew(i), size(T,2), size(T,1)],'EdgeColor','r', 'LineWidth',0.5);
    hold on;
end




fwhm_list = zeros(size(fineLoc)); %final position array
PSFdata = struct('PosSet',[],'FWHMset',[]);



for n = 1:numofbeads
    FwhmX = 0; FwhmY = 0; FwhmZ = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%
    bead1 = A(round(fineLoc(n,1)-size(T,1)/2) : round(fineLoc(n,1)+size(T,1)/2),...
        round(fineLoc(n,2)-1-size(T,2)/2) : round(fineLoc(n,2)-1+size(T,2)/2),...
        round(fineLoc(n,3)-size(T,3)/2) : round(fineLoc(n,3)+size(T,3)/2));
    bead1 = bead1/mean(bead1(:));
    %%%%%%%%%%%%%%%%%%%%%%%%%
    bead2 = A(round(fineLoc(n,1)-size(T,1)/2) : round(fineLoc(n,1)+size(T,1)/2),...
        round(fineLoc(n,2)+1-size(T,2)/2) : round(fineLoc(n,2)+1+size(T,2)/2),...
        round(fineLoc(n,3)-size(T,3)/2) : round(fineLoc(n,3)+size(T,3)/2));
    bead2 = bead2/mean(bead2(:));
    %%%%%%%%%%%%%%%%%%%%%%%%%
    bead3 = A(round(fineLoc(n,1)-1-size(T,1)/2) : round(fineLoc(n,1)-1+size(T,1)/2),...
        round(fineLoc(n,2)-size(T,2)/2) : round(fineLoc(n,2)+size(T,2)/2),...
        round(fineLoc(n,3)-size(T,3)/2) : round(fineLoc(n,3)+size(T,3)/2));
    bead3 = bead3/mean(bead3(:));
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    bead4 = A(round(fineLoc(n,1)+1-size(T,1)/2) : round(fineLoc(n,1)+1+size(T,1)/2),...
        round(fineLoc(n,2)-size(T,2)/2) : round(fineLoc(n,2)+size(T,2)/2),...
        round(fineLoc(n,3)-size(T,3)/2) : round(fineLoc(n,3)+size(T,3)/2));
    bead4 = bead4/mean(bead4(:));
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    bead5 = A(round(fineLoc(n,1)-size(T,1)/2) : round(fineLoc(n,1)+size(T,1)/2),...
        round(fineLoc(n,2)-size(T,2)/2) : round(fineLoc(n,2)+size(T,2)/2),...
        round(fineLoc(n,3)-1-size(T,3)/2) : round(fineLoc(n,3)-1+size(T,3)/2));
    bead5 = bead5/mean(bead5(:));
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    bead6 = A(round(fineLoc(n,1)-size(T,1)/2) : round(fineLoc(n,1)+size(T,1)/2),...
        round(fineLoc(n,2)-size(T,2)/2) : round(fineLoc(n,2)+size(T,2)/2),...
        round(fineLoc(n,3)+1-size(T,3)/2) : round(fineLoc(n,3)+1+size(T,3)/2));
    bead6 = bead6/mean(bead6(:));
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    beadstruc = struct('matrix',{bead1,bead2,bead3,bead4,bead5,bead6});
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:6
        bead = beadstruc(i).matrix;
        
        xline = 1:size(bead,1);
        yline = bead(:,ceil(size(bead,2)/2),ceil(size(bead,3)/2));
        
        
        try
            f = fit(xline.',yline,'gauss1');
            %     figure()
            %     plot(f,xline,yline)
            coeffvals = coeffvalues(f); %standard deviation of Normal distribution = coeffvals(3);
            FwhmX(i) = 2.3548*coeffvals(3)*806;
        catch
        end
        FwhmX = FwhmX(i) + FwhmX;
        
        
        xline = 1:size(bead,2);
        yline = bead(ceil(size(bead,1)/2),:,ceil(size(bead,3)/2));
        
        try
            f = fit(xline.',yline.','gauss1');
            %     figure()
            %     plot(f,xline,yline)
            coeffvals = coeffvalues(f);
            FwhmY(i) = 2.3548*coeffvals(3)*806;
        catch
        end
        FwhmY = FwhmY(i) + FwhmY;
        
        xline = 1:size(bead,3);
        yline = squeeze(bead(ceil(size(bead,1)/2),ceil(size(bead,2)/2),:));
        try
            f = fit(xline.',yline,'gauss1');
            %     figure()
            %     plot(f,xline,yline)
            coeffvals = coeffvalues(f);
            FwhmZ(i) = 2.3548*coeffvals(3)*806;
        catch
        end
        FwhmZ = FwhmZ(i) + FwhmZ;
    end
    
    fwhm_list(n,1) = FwhmX./6;
    fwhm_list(n,2) = FwhmY./6;
    fwhm_list(n,3) = FwhmZ./6;
end


for n = 1:numofbeads
    if ((fwhm_list(n,1)< 1000 || fwhm_list(n,1)> 3500)...
            || (fwhm_list(n,2)<1000 || fwhm_list(n,2)> 3500)...
            || (fwhm_list(n,3)<1000 || fwhm_list(n,3)>10000))
        fwhm_list(n,:) = 0;
        fineLoc(n,:) = 0;
    end
    
end

fwhm_list = fwhm_list(any(fwhm_list,2),:);
fineLoc = fineLoc(any(fineLoc,2),:);


corrected_fwhm_list = sqrt((fwhm_list).^2-(80).^2);
PSFdata.FWHMset = round(corrected_fwhm_list);
PSFdata.PosSet = fineLoc;


end
