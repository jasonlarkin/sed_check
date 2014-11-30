load('./nmd.mat');

[tmp,str.main]=system('pwd');

%--------------------------------------------------------------------------
    iseed = 1;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
    ikslice = 1;
%-------------------------------------------------------------------------- 

SED.SED(1:NMD.NUM_KPTS,1:(NMD.NUM_TSTEPS/2)) = 0.0;

%--------------------------------------------------------------------------
tic  
%--------------------------------------------------------------------------    
for ifft = 1:NMD.NUM_FFTS  
%VElOCITIES
    str_read=...
        strcat(...
        str.main ,'./dump_',int2str(iseed),'_',int2str(ifft),'.vel')
    fid=fopen('./dump_1_1.vel','r')
    dummy = textscan(fid,'%f%f%f','delimiter',' ','commentstyle', '--'); 
    fclose(fid);
%Store velocity data of all atoms: subtract off the last time step
    velx = zeros(NMD.NUM_ATOMS,NMD.NUM_TSTEPS);
    vely = zeros(NMD.NUM_ATOMS,NMD.NUM_TSTEPS);
    velz = zeros(NMD.NUM_ATOMS,NMD.NUM_TSTEPS);
%--------------------------------------------------------------------------
toc
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
tic  
%--------------------------------------------------------------------------
    for iatom = 1:NMD.NUM_ATOMS  
        velx(iatom,1:NMD.NUM_TSTEPS) =...
            dummy{1}...
            (iatom:NMD.NUM_ATOMS:(length(dummy{1}(:))-NMD.NUM_ATOMS));
        vely(iatom,1:NMD.NUM_TSTEPS) =...
            dummy{2}...
            (iatom:NMD.NUM_ATOMS:(length(dummy{1}(:))-NMD.NUM_ATOMS));
        velz(iatom,1:NMD.NUM_TSTEPS) =...
            dummy{3}...
            (iatom:NMD.NUM_ATOMS:(length(dummy{1}(:))-NMD.NUM_ATOMS));
    end
%--------------------------------------------------------------------------
toc
%--------------------------------------------------------------------------
%Remove dummy
    clear dummy  
%Set mass array
%     m = repmat(NMD.mass(:,1),1,NMD.NUM_TSTEPS);     
    m = NMD.mass(:,1);
%EIGENVECTORS
    eigenvec = NMD.eigvec;
%FREQUENCIES
%    freq = NMD.freq;              
%Zero main SED FP: this gets averaged as you loop over the NUM_FFTS      
    Q = zeros(1,NMD.NUM_TSTEPS);
    QDOT = zeros(1,NMD.NUM_TSTEPS);
    
SED.SED(...
    size(NMD.kptlist(:,1:3,ikslice),1),...
    1:(NMD.NUM_TSTEPS/2)) = 0.0;

%--------------------------------------------------------------------------
tic  
%--------------------------------------------------------------------------
    for ikpt = 1:size(NMD.kptlist(:,1:3,ikslice),1)
        for ib = 1:NMD.NUM_ATOMS_UCELL   
            spatial = 2*pi*1i*(...
            NMD.x0(ib:NMD.NUM_ATOMS_UCELL:end,3)*...
            ( (NMD.kptlist(ikpt,1,ikslice))/(NMD.alat*NMD.Nx) ) +...
            NMD.x0(ib:NMD.NUM_ATOMS_UCELL:end,4)*...
            ( (NMD.kptlist(ikpt,2,ikslice))/(NMD.alat*NMD.Ny) ) +...
            NMD.x0(ib:NMD.NUM_ATOMS_UCELL:end,5)*...
            ( (NMD.kptlist(ikpt,3,ikslice))/(NMD.alat*NMD.Nz) ) );
    
            kindex = NMD.kpt_index(ikpt,ikslice);

            QDOTX =...
                bsxfun(@times,...
                velx(ib:NMD.NUM_ATOMS_UCELL:end,:) , exp(spatial) );
            QDOTY =...
                bsxfun(@times,...
                vely(ib:NMD.NUM_ATOMS_UCELL:end,:) , exp(spatial) );
            QDOTZ =...
                bsxfun(@times,...
                velz(ib:NMD.NUM_ATOMS_UCELL:end,:) , exp(spatial) );

%            print("size QDOTX");
            size(QDOTX)

            QDOTX_FFT =sum(fft(QDOTX,[],2),1);
            QDOTY_FFT =sum(fft(QDOTY,[],2),1);
            QDOTZ_FFT =sum(fft(QDOTZ,[],2),1);

%            print("size QDOTX_FFT");
            size(QDOT)

%sum over x,y,z, and also b with SED.SED = SED.SED...
%            print("size SED.SED");
            size(SED.SED) 
            size(m(ib:NMD.NUM_ATOMS_UCELL:end))
            SED.SED =...
                SED.SED +...
                sqrt(m(ib)/NMD.NUM_UCELL_COPIES).*...
                (...
                real(QDOTX_FFT(1:(NMD.NUM_TSTEPS/2))).^2 +...
                imag(QDOTX_FFT(1:(NMD.NUM_TSTEPS/2))).^2 +...
                real(QDOTY_FFT(1:(NMD.NUM_TSTEPS/2))).^2 +...
                imag(QDOTY_FFT(1:(NMD.NUM_TSTEPS/2))).^2 +...
                real(QDOTZ_FFT(1:(NMD.NUM_TSTEPS/2))).^2 +...
                imag(QDOTZ_FFT(1:(NMD.NUM_TSTEPS/2))).^2 ...
                );

        end 
    end 
%--------------------------------------------------------------------------
toc 
%--------------------------------------------------------------------------
end %END ifft
    
%Average over FFTS
    SED.SED = SED.SED/NMD.NUM_FFTS;
%Define frequencies
    omega = (1:NMD.NUM_OMEGAS)*(NMD.w_max/NMD.NUM_OMEGAS);
%Output SED
    for ikpt = 1:size(NMD.kptlist(:,1:3,ikslice),1)
        str_write_single=...
            strcat('./SED_Phip_',...
            num2str(NMD.kptlist(ikpt,1,ikslice)),...
            num2str(NMD.kptlist(ikpt,2,ikslice)),...
            num2str(NMD.kptlist(ikpt,3,ikslice)),...
            '_',int2str(iseed),'.txt');
        output(1:length(omega),1) = omega;
        output(1:length(omega),2) = SED.SED(ikpt,:);
        dlmwrite(str_write_single,output,'delimiter',' ');
        clear output
    end %END ikpt    
%end %END iseed


