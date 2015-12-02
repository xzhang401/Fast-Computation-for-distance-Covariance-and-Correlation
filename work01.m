% work07_comp3methods

n = 20000; 
x = rand(n,1);
y = rand(n,1);

tic
cor3 = FaDCor(x,y);
toc
