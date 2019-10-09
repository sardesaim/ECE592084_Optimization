[X,Y] = meshgrid([-2:.25:2]);
% Plotting the Z-values of the function whose level sets have
% to be determined
Z = X.*exp(-X.^2-Y.^2);
% Plotting the Z-values of the function on which the level
% sets have to be projected
Z1 = X.^2+Y.^2;
% Plot your contour
[cm,c]=contour(X,Y,Z,30);
% Plot the surface on which the level sets have to be projected
s=surface(X,Y,Z1,'EdgeColor',[.8 .8 .8],'FaceColor','none')
% Get the handle to the children i.e the contour lines of the contour
cv=get(c,'children');
% Extract the (X,Y) for each of the contours and recalculate the
% Z-coordinates on the surface on which to be projected.
for i=1:length(cv)
      cc = cv(i);
      xd=get(cc,'XData');
      yd=get(cc,'Ydata');
      zd=xd.^2+yd.^2;
      set(cc,'Zdata',zd);
end
grid off
view(-15,25)
colormap cool