function stackPlots(fig, xlab)

bottomMargin=.05;
totalHeight=1-bottomMargin;
if isempty(fig)
    fig=gcf;
end
deselectAxes(fig)


if nargin<2 || isempty(xlab)
    
    xlab='';
end

 a=getAllAxes(fig);
 
 %%divide the height evenly 

 axesHeight= totalHeight/numel(a);
 
  positions =linspace(1-axesHeight, bottomMargin, numel(a));
  axorder= fliplr(1:numel(positions))
 for j= 1:numel(positions)
     
     set(a(axorder(j)), 'Position', [0, positions(j), 1, axesHeight])
    
      if j~=numel(positions)
         
          z=get(a(axorder(j)), 'XAxis');
          set(z, 'Visible', 'off');


      end
    
    
     if j==numel(positions)
         xlabel(xlab)
        showTicks(a(j));
        
     end         
 end
end

     
     
 