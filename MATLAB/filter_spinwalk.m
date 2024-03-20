function m1_t = filter_spinwalk(m1, xyz1, tissue_type, mask, fov)
%
% m1    : 3 * spins
% xyz1  : 3 * spins
% tissue_type : array of tissue types
% mask  : 
% fov   : in mm
%

m1   = squeeze(m1);
xyz1 = squeeze(xyz1);

if size(m1,1) ~= 3
    error('m1 must be of size 3*n')
end

if size(m1) ~= size(xyz1)
    error('size mismatch between m1 and xyz1')
end

sz      = size(mask);
scale   = sz(:) ./ fov(:);
I123    = floor(xyz1 .* scale) + 1;

ind     = sub2ind(sz, I123(1,:,:,:), I123(2,:,:,:), I123(3,:,:,:));
% ind_t   = mask(ind) == tissue_type;
ind_t   = ismember(mask(ind), tissue_type);
m1_t    = m1(:,ind_t);
