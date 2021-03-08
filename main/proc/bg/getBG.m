function bg = getBG(octa_ffname, scanObj, wb)
%getBG Returns a background vector for subtraction based on the mean signal
%of the whole volume

nFrames = scanObj.B*scanObj.xB;
averaging_mat = single(read_octa_frames(octa_ffname, scanObj, 1));
for ii=2:nFrames
    averaging_mat = averaging_mat + ...
        single(read_octa_frames(octa_ffname, scanObj, ii));
    
    if mod(ii,10)==0 && exist('wb', 'var') && ~isempty(wb)
        waitbar(ii/nFrames, wb, 'Calculating background');
    end
end
averaging_mat = averaging_mat./nFrames;
bg = mean(averaging_mat, 2);

end