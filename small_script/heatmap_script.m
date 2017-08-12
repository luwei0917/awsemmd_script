function heatmap_script_rmsd
a=MatrixConvert;
figure(1);

Z=linkage(a,'centroid')
c=dendrogram(Z);

figure(2);
CGobj = clustergram(a);
set(CGobj,'Colormap',redbluecmap,'DisplayRange',6.0:6.0);

function a=MatrixConvert
f=load('matrix')
size(f)
a=zeros(20, 20);
time=0;
for i=1:20
    for j=1:20
        time=time+1;
        a(i,j)=f(time);
    end
end

