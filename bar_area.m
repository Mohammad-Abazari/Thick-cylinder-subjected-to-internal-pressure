function x = bar_area(barnumber)
% gives crossectional area of a given bar
% usage: Ab = bar_area(barnumber)
    switch barnumber
        case 3; x = 0.11;
        case 4;  x = 0.20;
        case 5;  x = 0.31;
        case 6;  x = 0.44;
        case 7;  x = 0.60;
        case 8;  x = 0.79;
        case 9;  x = 1.00;
        case 10; x = 1.27;
        case 11; x = 1.56;
        case 14; x = 2.25;
        case 18; x = 4.00;
    end
end