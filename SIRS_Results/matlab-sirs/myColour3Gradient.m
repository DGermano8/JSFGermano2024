function c_map = myColour3Gradient(N,colour_1,colour_2,colour_3)
    length_p = N+1;

    length_12 = floor(length_p/2);
    length_23 = ceil(length_p/2)-1;
    
    length_1 = floor((2/3)*length_12);
    length_2 = ceil((1/3)*length_12);
    temp_colours_12 = colour_1 + 0.5*(colour_2-colour_1);
    
    colors_1t = [linspace(colour_1(1),temp_colours_12(1),length_1)', linspace(colour_1(2),temp_colours_12(2),length_1)', linspace(colour_1(3),temp_colours_12(3),length_1)'];
    colors_t2 = [linspace(temp_colours_12(1),colour_2(1),length_2)', linspace(temp_colours_12(2),colour_2(2),length_2)', linspace(temp_colours_12(3),colour_2(3),length_2)'];
    
    length_3 = floor((2/3)*length_23);
    length_4 = ceil((1/3)*length_23);
    temp_colours_23 = colour_2 + 0.5*(colour_3-colour_2);
    
    colors_2t = [linspace(colour_2(1),temp_colours_23(1),length_3)', linspace(colour_2(2),temp_colours_23(2),length_3)', linspace(colour_2(3),temp_colours_23(3),length_3)'];
    colors_t3 = [linspace(temp_colours_23(1),colour_3(1),length_4)', linspace(temp_colours_23(2),colour_3(2),length_4)', linspace(temp_colours_23(3),colour_3(3),length_4)'];
    
    colors_p = [colors_1t; colors_t2; colors_2t; colors_t3];
    c_map = colors_p;
end

