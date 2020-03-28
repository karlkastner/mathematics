% Sun 12 Jan 19:40:40 +08 2020
a = 0.25;
n = 1e7;
sk = 0.99;
x = skewrnd_walsh(5,2,sk,n);
[mean(x), std(x), skewness(x)]

