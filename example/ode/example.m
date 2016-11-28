opt=[];
opt.RelTol=1e-8;opt.AbsTol=1e-8;
figure(1);
tic;
[tNodes,xNodes,stats]=dopri5Mex(@f,[0,4, 8],1,opt);
toc;
% also possible (auch möglich):
% [tNodes,xNodes,stats]=dopri5Mex('f',[0,3],1,opt);
% [tNodes,xNodes,stats]=dopri5Mex(inline('x','t','x'),[0,3],1,opt);
t=linspace(0,8,200);
plot(t,exp(t),'b-',tNodes,xNodes,'rx');
figure(2);
plot(tNodes,exp(tNodes) - xNodes,'b-');
figure(3);
tic;
[tNodes,xNodes,stats]=dop853Mex(@f,[0,4, 8],1,opt);
toc;
% also possible (auch möglich):
% [tNodes,xNodes,stats]=dopri5Mex('f',[0,3],1,opt);
% [tNodes,xNodes,stats]=dopri5Mex(inline('x','t','x'),[0,3],1,opt);
t=linspace(0,8,200);
plot(t,exp(t),'b-',tNodes,xNodes,'rx');
figure(4);
plot(tNodes,exp(tNodes) - xNodes,'b-');
figure(5)
options = odeset('RelTol',opt.RelTol,'AbsTol',opt.AbsTol);
tic;
[t,y]=ode45(@f,[0, 4, 8],1,options);
toc;
plot(t,exp(t),'b-',t,y,'rx');
figure(6);
plot(t,exp(t) - y,'b-');

