function [row_shift,col_shift] = fftreg(buf1ft,buf2ft)
% calculates the registration shift via discrete Fourier transform.
%
% [row_shift,col_shift] = fftreg(buf1ft,buf2ft) calculates the row shift
% and the column shift for two Fourier transformed arrays.

[m,n]=size(buf1ft);
CC = ifft2(buf1ft.*conj(buf2ft));
[max1,loc1] = max(CC);
[~,loc2] = max(max1);
rloc=loc1(loc2);
cloc=loc2;

md2 = fix(m/2);
nd2 = fix(n/2);

if rloc > md2
    row_shift = rloc - m - 1;
else
    row_shift = rloc - 1;
end

if cloc > nd2
    col_shift = cloc - n - 1;
else
    col_shift = cloc - 1;
end
