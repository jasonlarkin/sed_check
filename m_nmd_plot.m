Phi = load('./lj_copy/SED_Phi_200_1.txt');
Phip = load('./lj_copy/SED_Phip_200_1.txt');

semilogy(...
Phi(1:2048,1),sum(Phi(1:2048,:),2),...
Phip(1:2048,1),Phip(1:2048,2)/(4*pi)...
)

xlabel ("\omega (LJ units)");
ylabel ("PS (LJ units)");
legend ("\Phi", "\Phi'");

print -deps -color sed_check.eps

%sum(...
%Phi(1:2048,1),sum(Phi(1:2048,:),2) - 
%Phip(1:2048,1),Phip(1:2048,2)/(4*pi)...
%)

