load('./nmd.mat');

[tmp,str.main]=system('pwd');

%--------------------------------------------------------------------------
    iseed = 1;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
    ikslice = 1;
%-------------------------------------------------------------------------- 

SED.SED(1:NMD.NUM_KPTS,1:(NMD.NUM_TSTEPS/2),1:NMD.NUM_MODES) = 0.0;

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
    1:(NMD.NUM_TSTEPS/2),1:NMD.NUM_MODES) = 0.0;

%--------------------------------------------------------------------------
tic  
%--------------------------------------------------------------------------
    for ikpt = 1:size(NMD.kptlist(:,1:3,ikslice),1)
        for imode = 1:NMD.NUM_MODES
    
            spatial = 2*pi*1i*(...
    NMD.x0(:,3)*( (NMD.kptlist(ikpt,1,ikslice))/(NMD.alat*NMD.Nx) ) +...
    NMD.x0(:,4)*( (NMD.kptlist(ikpt,2,ikslice))/(NMD.alat*NMD.Ny) ) +...
    NMD.x0(:,5)*( (NMD.kptlist(ikpt,3,ikslice))/(NMD.alat*NMD.Nz) ) );
    
    kindex = NMD.kpt_index(ikpt,ikslice);
            
            eigx = repmat(...
                conj(...
                eigenvec(...
                ((NMD.NUM_ATOMS_UCELL*3)*(kindex-1)+1)... 
                :3:...
                ((NMD.NUM_ATOMS_UCELL*3)*kindex),imode...
                )...
                ),NMD.NUM_UCELL_COPIES,1);
            
            eigy = repmat(... 
                conj(...
                eigenvec(...
                ((NMD.NUM_ATOMS_UCELL*3)*(kindex-1)+2)... 
                :3:...
                ((NMD.NUM_ATOMS_UCELL*3)*kindex),imode...
                )...
                ),NMD.NUM_UCELL_COPIES,1);
            
            eigz = repmat(...
                conj(...
                eigenvec(...
                ((NMD.NUM_ATOMS_UCELL*3)*(kindex-1)+3)... 
                :3:...
                ((NMD.NUM_ATOMS_UCELL*3)*kindex),imode...
                )...
                ),NMD.NUM_UCELL_COPIES,1);

            QDOT = sum(...
                bsxfun(@times,...
                bsxfun(@times, velx, eigx) + ...
                bsxfun(@times, vely, eigy) + ...
                bsxfun(@times, velz, eigz) ...
                , exp(spatial).*(sqrt(m/NMD.NUM_UCELL_COPIES)) )...
                , 1 );

            KEFFT = real(fft(QDOT)).^2 + imag(fft(QDOT)).^2;
            
        SED.SED(ikpt,:,imode) =...
            SED.SED(ikpt,:,imode)+KEFFT(1:(NMD.NUM_TSTEPS/2)) ;
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
            strcat('./SED_Phi_',...
            num2str(NMD.kptlist(ikpt,1,ikslice)),...
            num2str(NMD.kptlist(ikpt,2,ikslice)),...
            num2str(NMD.kptlist(ikpt,3,ikslice)),...
            '_',int2str(iseed),'.txt');
        output(1:length(omega),1) = omega;
        output(1:length(omega),2:(NMD.NUM_MODES+1)) = SED.SED(ikpt,:,:);
        dlmwrite(str_write_single,output,'delimiter',' ');
        clear output
    end %END ikpt    
%end %END iseed

