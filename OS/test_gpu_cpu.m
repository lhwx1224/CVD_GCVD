clear
clc
% Configure System Information
disp(['========= CPU Device Information =========='])
system('wmic cpu get caption');
system('wmic cpu get name');
system('wmic cpu get numberofcores');
system('wmic cpu get maxclockspeed');
gpuinfo = gpuDevice();
disp(['========= GPU Device Information =========='])
disp(['GPU Name:',gpuinfo.Name])
disp(['GPU Available Memory:',num2str(gpuinfo.AvailableMemory/1024/1024),'MB'])
disp(['GPU NO. of Cores:',num2str(gpuinfo.MultiprocessorCount)])
disp(['GPU Clock Rate:',num2str(gpuinfo.ClockRateKHz/1000),'MHz'])
msize = [10 100 500 1000 5000 10000];
if (gpuinfo.AvailableMemory)/1024/1024 > 6000
    M = 6;
    disp(['Maximum Matrix Size: ',num2str(msize(M))])
else
    M = 4;
    disp(['Maximum Matrix Size: ',num2str(msize(M))])
end

%%
disp('======== Running GPU-CPU Computation =========')
N = 10;
disp('======== Example 1: Matrix Power =========')
mcput = zeros(M,1);
mgput = zeros(M,1);
for j = 1:M
    cput = zeros(N,1);
    gput = zeros(N,1);
    for i = 1:N
        A = rand(msize(j),msize(j));
        tic
        C = A^3;
        cput(i) = toc;
    end
    
    for i = 1:N
        A = gpuArray(rand(msize(j),msize(j)));
        tic
        C = A^3;
        gput(i) = toc;
    end

mcput(j) = mean(cput(2:N));
mgput(j) = mean(gput(2:N));
goc(j) = mcput(j)/mgput(j);
disp('======== Example 1: Matrix Power Done !=========')
disp(['Matrix A Size: ',num2str(size(A))])
disp(['CPU time: ',num2str(mcput(j)),' Seconds'])
disp(['GPU time: ',num2str(mgput(j)),' Seconds'])
disp(['GPU speedup over CPU: ',num2str(goc(j))])
end

figure(1),clf
x = 1:M;
yyaxis left
plot(x,mcput,'o'), hold on, plot(x,mgput,'^')
set(gca,'Yscale','log')
yyaxis right
plot(x,goc,'+')
hold on
plot([1-0.1 M+0.1],[1 1])
set(gca,'Yscale','log')
xlabel('Square Matrix Size')
xticks([x])
xticklabels(split(num2str(msize)))
xlim([1-0.1 M+0.1])
legend('CPU time','GPU time','R = GPU/CPU','R = 1','Location','Northwest')
pbaspect([1.618 1 1])
set(gcf,'papersize',[6 6*0.618])
set(gcf,'paperposition',[0 0 6 6*0.618])
%%
N = 10;
disp('======== Example 2: Matrix Inverse =========')
msize = [10 100 500 1000 5000 10000];
mcput = zeros(M,1);
mgput = zeros(M,1);
for j = 1:M
    cput = zeros(N,1);
    gput = zeros(N,1);
    for i = 1:N
        A = rand(msize(j),msize(j));
        tic
        C = inv(A);
        cput(i) = toc;
    end
    for i = 1:N
        A = gpuArray(rand(msize(j),msize(j)));
        tic
        C = inv(A);
        gput(i) = toc;
    end
    mcput(j) = mean(cput(2:N));
    mgput(j) = mean(gput(2:N));
    goc(j) = mcput(j)/mgput(j);
    disp('======== Example 2: Matrix Inverse Done !=========')
    disp(['Matrix A Size: ',num2str(size(A))])
    disp(['CPU time: ',num2str(mcput(j)),' Seconds'])
    disp(['GPU time: ',num2str(mgput(j)),' Seconds'])
    disp(['GPU speedup over CPU: ',num2str(goc(j))])
end

figure(2),clf
x = 1:M;
yyaxis left
plot(x,mcput,'o'), hold on, plot(x,mgput,'^')
set(gca,'Yscale','log')
yyaxis right
plot(x,goc,'+')
hold on
plot([1-0.1 M+0.1],[1 1])
set(gca,'Yscale','log')
xlabel('Square Matrix Size')
xticks([x])
xticklabels(split(num2str(msize)))
xlim([1-0.1 M+0.1])
legend('CPU time','GPU time','R = GPU/CPU','R = 1','Location','Northwest')
pbaspect([1.618 1 1])
%%
disp('======== Example 3: Matrix SVD =========')
msize = [10 100 500 1000 5000 10000];
mcput = zeros(M,1);
mgput = zeros(M,1);
for j = 1:M
for i = 1:N
    A = rand(msize(j),msize(j));
    tic
    C = svd(A);
    cput(i) = toc;
end
for i = 1:N
    A = gpuArray(rand(msize(j),msize(j)));
    tic
    C = svd(A);
    gput(i) = toc;
end
    mcput(j) = mean(cput(2:N));
    mgput(j) = mean(gput(2:N));
    goc(j) = mcput(j)/mgput(j);
    disp('======== Example 3: Matrix SVD Done !=========')
    disp(['Matrix A Size: ',num2str(size(A))])
    disp(['CPU time: ',num2str(mcput(j)),' Seconds'])
    disp(['GPU time: ',num2str(mgput(j)),' Seconds'])
    disp(['GPU speedup over CPU: ',num2str(goc(j))])
end

figure(3),clf
x = 1:M;
yyaxis left
plot(x,mcput,'o'), hold on, plot(x,mgput,'^')
set(gca,'Yscale','log')
yyaxis right
plot(x,goc,'+')
hold on
plot([1-0.1 M+0.1],[1 1])
set(gca,'Yscale','log')
xlabel('Square Matrix Size')
xticks([x])
xticklabels(split(num2str(msize)))
xlim([1-0.1 M+0.1])
legend('CPU time','GPU time','R = GPU/CPU','R = 1','Location','Northwest')
pbaspect([1.618 1 1])
%%
disp('======== Example 4: 2D-Convolution =========')
msize = [10 50 100 500 1000 1500 2000];
mcput = zeros(M,1);
mgput = zeros(M,1);
for j = 1:M
    for i = 1:N
        A = rand(msize(j),msize(j));
        B = rand(msize(j),msize(j));
        
        tic
        C = conv2(A,B);
        cput(i) = toc;
    end
    
    for i = 1:N
        A = gpuArray(rand(msize(j),msize(j)));
        B = gpuArray(rand(msize(j),msize(j)));
        tic
        C = conv2(A,B);
        gput(i) = toc;
    end
    mcput(j) = mean(cput(2:N));
    mgput(j) = mean(gput(2:N));
    goc(j) = mcput(j)/mgput(j);
    disp('======== Example 4: 2D-Convolution Done !=========')
    disp(['Matrix A Size: ',num2str(size(A))])
    disp(['CPU time: ',num2str(mcput(j)),' Seconds'])
    disp(['GPU time: ',num2str(mgput(j)),' Seconds'])
    disp(['GPU speedup over CPU: ',num2str(goc(j))])
end

figure(4),clf
x = 1:M;
yyaxis left
plot(x,mcput,'o'), hold on, plot(x,mgput,'^')
set(gca,'Yscale','log')
yyaxis right
plot(x,goc,'+')
hold on
plot([1-0.1 M+0.1],[1 1])
set(gca,'Yscale','log')
xlabel('Square Matrix Size')
xticks([x])
xticklabels(split(num2str(msize)))
xlim([1-0.1 M+0.1])
legend('CPU time','GPU time','R = GPU/CPU','R = 1','Location','Northwest')
pbaspect([1.618 1 1])
%%
disp('======== Example 5: 2D-FFT =========')
msize = [10 100 500 1000 5000 10000];
mcput = zeros(M,1);
mgput = zeros(M,1);
for j = 1:M
for i = 1:N
    A = rand(msize(j),msize(j));
    tic
    C = fft2(A);
    cput(i) = toc;
end
for i = 1:N
    A = gpuArray(rand(msize(j),msize(j)));
    tic
    C = fft2(A);
    gput(i) = toc;
end
    disp('======== Example 5: 2D-FFT Done !=========')
    disp(['Matrix A Size: ',num2str(size(A))])
    disp(['CPU time: ',num2str(mcput(j)),' Seconds'])
    disp(['GPU time: ',num2str(mgput(j)),' Seconds'])
    disp(['GPU speedup over CPU: ',num2str(goc(j))])
end

figure(5),clf
x = 1:M;
yyaxis left
plot(x,mcput,'o'), hold on, plot(x,mgput,'^')
set(gca,'Yscale','log')
yyaxis right
plot(x,goc,'+')
hold on
plot([1-0.1 M+0.1],[1 1])
set(gca,'Yscale','log')
xlabel('Square Matrix Size')
xticks([x])
xticklabels(split(num2str(msize)))
xlim([1-0.1 M+0.1])
legend('CPU time','GPU time','R = GPU/CPU','R = 1','Location','Northwest')
pbaspect([1.618 1 1])