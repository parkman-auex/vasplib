function sym_out=double2sym(double_in,accuracy)
if nargin < 2
    accuracy =5;
end
if isnan(double_in )
    sym_out = nan;
    return;
end



double_single = roundn(double_in^1,-accuracy);
double_square = roundn(double_in^2,-accuracy);
double_cubic =  roundn(double_in^3,-accuracy);


double_single_out  = double2sym_one(double_single,accuracy);
double_square_out  = double2sym_one(double_square,accuracy);
double_cubic_out   = double2sym_one(double_cubic,accuracy);

double_out = [double_single_out;double_square_out;double_cubic_out ];
double_out_length = [length(char(string(double_single_out)));...
    length(char(string(double_square_out)));...
    length(char(string(double_cubic_out)))];

[~,min_label ]= min(double_out_length );

sym_out = double_out(min_label)^(1/min_label);


end

function sym_out = double2sym_one(double_in,accuracy)
%disp(double_in)
integer = fix(double_in);
digit = double_in-fix(double_in);

sym_out = sym(integer )+1/sym(roundn(1/digit,-accuracy+2));
%disp(sym_out);
end
