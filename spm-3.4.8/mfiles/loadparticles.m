

    %load data
    fid = fopen(['./output/particles',num2str(fnr),'.dat']);
    npa = fread(fid,1,'long');
    xp = fread(fid,npa,'double');
    yp = fread(fid,npa,'double');
    zp = fread(fid,npa,'double');
    bf = fread(fid,npa,'double');
    sedi = fread(fid,npa,'double');
    N10 = fread(fid,npa,'double');
    bx = fread(fid,npa,'double');
    by = fread(fid,npa,'double');
    birthday = fread(fid,npa,'double');
    age = fread(fid,npa,'double');
    dl = fread(fid,npa,'double');
    vx = fread(fid,npa,'double');
    vy = fread(fid,npa,'double');
    vz = fread(fid,npa,'double');
    erate = fread(fid,npa,'double');
    fclose(fid);

