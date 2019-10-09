figure(brainfig)
axis vis3d
for ii = 1:720
   %axis([-1 1 -1 1])
   camorbit(0.5,0.5,'data',[0 1 0])
   drawnow
   M(ii)=getframe(gcf,[74 47 450 350]);
end
%%
v=VideoWriter('C:/temp/test_probes.mp4');
v.Quality=100;
v.open()
v.writeVideo(M)
v.close()