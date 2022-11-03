img = imread('~/Pictures/totom-totoro-final.png');
 img = double(mean(img(:,:,:,1),3));
 img=img(1:723,:);
 i0=img;
a=deg2rad(40);
 s=size(img);
 R = imrotmat(s,a);
 Rt=imrotmat(s,-a);
 img = R*Rt*flat(img);
 img=reshape(img,s);
 clf;
 imagesc(img-i0)
