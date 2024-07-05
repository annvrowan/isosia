

    %load data
    fid = fopen(['./output/output',num2str(fnr),'.dat']);
    xc = fread(fid,[ny,nx],'double');
    yc = fread(fid,[ny,nx],'double');
    bed = fread(fid,[ny,nx],'double');
    ice = fread(fid,[ny,nx],'double');
    bslope = fread(fid,[ny,nx],'double');
    tslope = fread(fid,[ny,nx],'double');
    te2 = fread(fid,[ny,nx],'double');
    frost = fread(fid,[ny,nx],'double');
    sediment = fread(fid,[ny,nx],'double');
    tn = fread(fid,[ny,nx],'double');
    te = fread(fid,[ny,nx],'double');
    ts = fread(fid,[ny,nx],'double');
    sliding = fread(fid,[ny,nx],'double');
    deformation = fread(fid,[ny,nx],'double');
    pw = fread(fid,[ny,nx],'double');
    bmelt = fread(fid,[ny,nx],'double');
    quarrying = fread(fid,[ny,nx],'double');
    accrate = fread(fid,[ny,nx],'double');
    meltrate = fread(fid,[ny,nx],'double');
    frostrate = fread(fid,[ny,nx],'double');
    fluvial = fread(fid,[ny,nx],'double');
    landslide = fread(fid,[ny,nx],'double');
    abrasion = fread(fid,[ny,nx],'double');
    isostasy = fread(fid,[ny,nx],'double');
    lee = fread(fid,[ny,nx],'double');
    Ts = fread(fid,[ny,nx],'double');
    hillslope = fread(fid,[ny,nx],'double');
    Ta = fread(fid,[ny,nx],'double');
    Tb = fread(fid,[ny,nx],'double');
    sfac = fread(fid,[ny,nx],'double');
    margin = fread(fid,[ny,nx],'int');
    for i=1:20,
        Vs{i} = fread(fid,[ny,nx],'double');
    end
    fclose(fid);
