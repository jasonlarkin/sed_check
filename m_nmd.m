clear
%m_nmd_load.m

NMD.str.main

NMD.alat = 1.5560;
NMD.ucell = [
0.0  0.0  0.0
0.5  0.5  0.0
0.5  0.0  0.5
0.0  0.5  0.5
];
NMD.Nx = 4; NMD.Ny = NMD.Nx; NMD.Nz = NMD.Nx;
NMD.x0 = load('./lj_copy/x0.data');
NMD.NUM_ATOMS = size(NMD.x0,1);
NMD.NUM_MODES = size(NMD.ucell,1)*3;
NMD.NUM_ATOMS_UCELL = 4;
NMD.NUM_UCELL_COPIES = NMD.NUM_ATOMS / NMD.NUM_ATOMS_UCELL;

NMD.mass = NMD.x0(:,2);

NMD.kptlist = [2 0 0]; %0.5 0 0 in GULP
NMD.NUM_KPTS = size(NMD.kptlist,1);
NMD.kpt_index = 1;

NMD.NUM_SEEDS=1;

NMD.NUM_TSTEPS = 2^17 / 2^4;
NMD.NUM_FFTS = 1;
NMD.NUM_OMEGAS = NMD.NUM_TSTEPS/2;
NMD.w_max = 2*pi / (2^4*2);

NMD.eigvec = dlmread('./lj_copy/eigvec.dat');

save('./lj_copy/nmd.mat');
