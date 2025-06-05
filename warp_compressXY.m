%dum=imread('4d_noChiron_S01_frame462.tif');
%[i1,i2]=size(dum);
%subplot(1,2,1), imagesc(dum)
%dum2=reshape(dum,i1*2,i2/2);
%subplot(1,2,2), imagesc(dum2(1:i1,1:i2/2))


%addpath 'C:\Users\endrejm\OneDrive - Universitetet i Oslo\Documents\Richard'

V=double(tiffreadVolume('4d_noChiron_S01_711_decon_stitched_aligned-c0_8bit.tif'));
[i1,i2,i3]=size(V);
X=1:i2;
Y=1:i1;
Z=1:i3;
xstep = 0.4228; %0.234216/0.554
ystep=xstep;
%zstep=2.36;
%aspect ratio = 2.36
%[Xq,Yq,Zq] = meshgrid(X,Y,1:1/zstep:i3);
[Xq,Yq,Zq] = meshgrid(1:1/xstep:i2,1:1/ystep:i1,Z);
Vq = uint8(interp3(X,Y,Z,V,Xq,Yq,Zq));
[numRows,numCols,numSlices]=size(Vq);
% for slice = 1:numSlices
%     sliceData = Vq(:, :, slice);
%     imwrite(sliceData,['4d_noChiron_S01_interpol',int2str(slice),'.tif'])
% end

t = Tiff('4d_noChiron_S01_711_decon_stitched_aligned-c0_8bit_warped_xy.tif', 'w');
%t = Tiff('4d_noChiron_S01_interpol.tif', 'w');
for slice = 1:numSlices
    sliceData = Vq(:, :, slice);
    
    % Write the current slice to the multi-page TIFF stack
    t.setTag('Photometric', Tiff.Photometric.MinIsBlack);
    t.setTag('Compression', Tiff.Compression.None);
    t.setTag('BitsPerSample', 8); % Set BitsPerSample to 32 for single precision data
    % t.setTag('SamplesPerPixel', 1);
    % t.setTag('SampleFormat', Tiff.SampleFormat.IEEEFP);
    t.setTag('ImageLength', numRows);
    t.setTag('ImageWidth', numCols);
    t.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
    t.write(sliceData);
    
    % Add a new page to the multi-page TIFF stack
    if slice ~= numSlices
        t.writeDirectory();
    end

end
t.close();