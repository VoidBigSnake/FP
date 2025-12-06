% CLONALG - Clonal Selection Algorithm for Optimization Problems
% For 3D MINIMIZATION Problems

% OUTPUTS:

% x,y,fx	-> fx = f(x,y) is the minimal value of the function
% vfx		-> Best fitness of population through generations

% Reference:
% de Castro, L. N., and Von Zuben, F. J. Learning and Optimization Using
% the Clonal Selection Principle. IEEE Transactions on Evolutionary
% Computation. v. 6(3). 2002. DOI: 10.1109/TEVC.2002.1011539

clear
close
clc

rng('shuffle') %default时会输出一样的结果

N    = 6;                  % Population size
Ab   = cadeia(N,44);         % Antibody population
gen  = 5;                   % Number of generations
pm   = 0.2;                  % Mutation probability
d    = 0.25;                  % Population to suffer random reshuffle %保证N*d>1
beta = 0.5;                  % Proportion of clones %保证N*beta>1 否则不足一个克隆，会报错

% Function to optimization
f = @myfunc;

varMin = 20; % Lower bound
varMax =  40; % Upper bound

fbest = 0; % Global f best (minimum)

% x = meshgrid(linspace(varMin, varMax, 61));
% y = meshgrid(linspace(varMin, varMax, 61))';
% 
% vxp = x;
% vyp = y;
% 
% vzp = f([x(:),y(:)]);
% vzp = reshape(vzp,size(x));

% 先随便给一个占位，以便 imprime 不报错（反正 PRINT=0 不会画）
vxp = 0; vyp = 0; vzp = 0;

x = decode(Ab(:,1:22),varMin,varMax);
y = decode(Ab(:,23:end),varMin,varMax);

fit = f([x(:),y(:)]);

figure
imprime(0,vxp,vyp,vzp,x,y,fit)

% Hypermutation controlling parameters
pma  = pm;
itpm = gen;
pmr  = 0.8;

% General defintions
vfx   = zeros(gen,1);
PRINT = 0;          % 不画 3D 曲面，避免卡 UI

% Generations
for it = 1:gen
    
    % Decode (Step 2)
    x = decode(Ab(:,1:22),varMin,varMax);
    y = decode(Ab(:,23:end),varMin,varMax);
    
    fit = f([x(:),y(:)]);
    
    [a,ind] = sort(fit);
    
    % Select (Step 3)
    valx = x(ind);
    valy = y(ind);
    
    imprime(PRINT,vxp,vyp,vzp,x,y,fit);

    % Clone (Step 4)
    [T,pcs] = reprod(N,beta,ind,Ab);

    % Hypermutation (Step 5)
    % - 逐位生成与 T 同尺寸的掩码 M，概率 pm 置 1，相当于标记“需要翻转的基因位”。
    % - T = T - 2 .* (T.*M) + M 利用 0/1 的代数关系执行批量比特翻转：
    %   * 当 M(i,j)=0 时，T 不变；
    %   * 当 M(i,j)=1 时，T(i,j) = 1 - T(i,j)，即 1→0、0→1。
    % - 超变异后，将每个克隆批次的末尾个体（由 pcs 给出其全局下标）替换为对应原型
    %   Ab(fliplr(ind),:)，确保每批都保留“原型对照”，便于后续在批次内部做择优比较。
    M = rand(size(T,1),44) <= pm;
    T = T - 2 .* (T.*M) + M;

    T(pcs,:) = Ab(fliplr(ind),:);
    
    % Decode (Step 6)
    x = decode(T(:,1:22),varMin,varMax);
    y = decode(T(:,23:end),varMin,varMax);
    
    fit = f([x(:),y(:)]);
    
    pcs = [0 pcs];

    for i = 1:N
        % Re-Selection (Step 7)
        [~,bcs(i)] = min(fit(pcs(i)+1:pcs(i+1)));
        bcs(i) = bcs(i) + pcs(i);
    end

    % Insert
    Ab(fliplr(ind),:) = T(bcs,:);

    % Editing (Repertoire shift)
    nedit = round(d*N);

    % Replace (Step 8)
    Ab(ind(end-(nedit-1):end),:) = cadeia(nedit,44);

    pm = pmcont(pm,pma,pmr,it,itpm);
    
    vfx(it) = a(1);
        
    % fprintf('%2d  f(%6.2f,%6.2f): %7.2f\n',it,valx(1),valy(1),vfx(it))

end

% Minimization problem
x  = valx(1);
y  = valy(1);
fx = vfx(end);

% Plot
figure
semilogy(vfx)
title('Minimization')
xlabel('Iterations')
ylabel('Best f(x,y)')
grid on

txt2 = ['F Best: ', num2str(fbest)];
text(0,1,txt2,'Units','normalized',...
     'HorizontalAlignment','left','VerticalAlignment','bottom');

txt3 = ['F Found: ', num2str(fx)];
text(1,1,txt3,'Units','normalized',...
     'HorizontalAlignment','right','VerticalAlignment','bottom');



% INTERNAL FUNCTIONS

function imprime(PRINT,vx,vy,vz,x,y,fx)

    if PRINT == 1
        meshc(vx,vy,vz)
        hold on
        title('Minimization')
        xlabel('x')
        ylabel('y')
        zlabel('f(x,y)')
        plot3(x,y,fx,'k*')
        colormap jet
        drawnow
        hold off
        pause(0.1)
    end

end

function [T,pcs] = reprod(N,beta,ind,Ab)

    % N    -> number of clones
    % beta -> multiplying factor
    % ind  -> best individuals
    % Ab   -> antibody population

    % T    -> temporary population
    % pcs  -> final position of each clone

    % 说明：
    % - 按当前个体优劣顺序（ind 从差到优）逐个取出，克隆数约为 beta*N。
    % - cs(i) 记录第 i 个原型的克隆数量；pcs(i) 为累积和，标记该批克隆在
    %   临时族群 T 中的结尾位置，便于后续在每个批次内部选择最优突变体。
    % - T 通过纵向拼接各批克隆得到：第 i 批等于第 i 个原型的重复行向量。

    T = [];

   for i = 1:N
      % round -> MATLAB 的“四舍五入”函数（0.5 远离 0 取整），用来把 beta*N 变成整数
      % 克隆份数；若 beta*N=3.2 得到 3，若=3.5 得到 4。
      cs(i) = round(beta*N);
      pcs(i) = sum(cs);
      T = [T; ones(cs(i),1) * Ab(ind(end-i+1),:)];
   end

end

function pm = pmcont(pm,pma,pmr,it,itpm)

    % pma  -> initial value
    % pmr  -> control rate
    % itpm -> iterations for restoring
    
    if rem(it,itpm) == 0
       pm = pm * pmr;
       if rem(it,10*itpm) == 0
          pm = pma;
       end
    end

end

function z = decode(Ab,varMin,varMax)

    % x	-> real value (precision: 6)
    % v	-> binary string (length: 22)
    
    Ab = fliplr(Ab);
    s = size(Ab);
    aux = 0:1:21;
    aux = ones(s(1),1)*aux;
    x1 = sum((Ab.*2.^aux),2);
    
    % Keeping values between bounds
    z = varMin + x1' .* (varMax - varMin)/(2^22 - 1);

end

function Ab = cadeia(n1,s1)

    % Antibody (Ab) chains
    Ab = 2 .* rand(n1,s1) - 1;
    Ab = hardlim(Ab);

end
