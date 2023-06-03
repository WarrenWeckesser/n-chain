
figure(1)
clf
hold on
for k = 1:6
curve = ReadNChainSolutionBranch('autodata/s.n05',k);
plot(curve(:,1),curve(:,4))
pause
end

