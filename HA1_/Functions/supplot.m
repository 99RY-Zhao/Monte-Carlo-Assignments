function [] = supplot(x,y,min,max,color)
    low = x(find(y ~= 0, 1, "first"));
    upp = x(find(y ~= 0, 1, "last"));

    plot([low low],[min max],":","Color",color,"HandleVisibility","off");
    plot([upp upp],[min max],":","Color",color,"HandleVisibility","off");
end

