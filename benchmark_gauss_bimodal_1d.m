function fx=benchmark_gauss_bimodal_1d(xx,mu,sigma,delta) %multimodal

%the function is one-dimensional (length(x)=1)
%sigma and mu have length equal to the number of peaks
%MU: position of the peak
%DELTA: set the heights of the peaks
%SIGMA: standard deviation (sqrt(variance))

mu=[-14;-10.5;-7;-3;2;4;12;21;26;32]; %modality X Nx
sigma=[1.8;1.8;1.8;0.5;1;1;3;2;1;1]; % modality X Nx
delta=[1.1;1.1;1.1;1.7;1;1;1.5;1.55;1.6;1.2];

fx=0;
for i=1:length(mu) %modality
    x=(xx-mu(i,:))^2./sigma(i,:)^2;
    fx=fx+delta(i)*exp(-1/2*(x));
end
