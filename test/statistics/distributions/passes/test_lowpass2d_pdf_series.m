% 2023-03-03 21:02:29.829064145 +0100
% the series expansion has only a very small roc, so not useful
a=1;
 order=1;
 L=1/(2*pi);
 n=100;
 fr = linspace(0,L,n)';
 S = lowpass2d_pdf(fr,a,order);
 S(:,2) = lowpass2d_pdf_b(fr,a,order,10)/(2*pi);
 plot(S)
