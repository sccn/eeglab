function [IDX,C] = neural_net(clustdata,clus_num)

nmin = min(clustdata);
nmax = max(clustdata);
net = newc([nmin ;nmax].',clus_num);
net = train(net,(clustdata).');
Y = sim(net,(clustdata).');
IDX = vec2ind(Y);
C = zeros(clus_num,size(clustdata,2));
for k = 1:clus_num
    C(k,:) = sum(clustdata(find(IDX == k),:));
end