function xyz = write_xyz0(fov, n_spin, filename)

grid_size = round(n_spin ^ (1/3));

[x, y, z] = ndgrid(linspace(0.01*fov(1), 0.99*fov(1), grid_size), ...
                   linspace(0.01*fov(2), 0.99*fov(2), grid_size), ...
                   linspace(0.01*fov(3), 0.99*fov(3), grid_size));

xyz = transpose([x(:), y(:), z(:)]);
xyz = single(xyz);

fileID = fopen(filename, 'w');
fwrite(fileID, xyz, 'single'); % x1,y1,z1, x2,y2,z2, ..., xn,yn,zn
fclose(fileID);
