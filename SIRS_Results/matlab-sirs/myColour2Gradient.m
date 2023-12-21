function c_map = myColour2Gradient(N,colour_1,colour_2)
    length_p = N;

    length_12 = floor(length_p);
    
    length_1 = floor((1/2)*length_12);
    length_2 = ceil((1/2)*length_12);
    temp_colours_12 = colour_1 + 0.5*(colour_2-colour_1);
    
    colors_1t = [linspace(colour_1(1),temp_colours_12(1),length_1)', linspace(colour_1(2),temp_colours_12(2),length_1)', linspace(colour_1(3),temp_colours_12(3),length_1)'];
    colors_t2 = [linspace(temp_colours_12(1),colour_2(1),length_2)', linspace(temp_colours_12(2),colour_2(2),length_2)', linspace(temp_colours_12(3),colour_2(3),length_2)'];
    
    
    colors_p = [colors_1t; colors_t2];
    c_map = colors_p;
end

