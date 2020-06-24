function [ string ] = sym2str( symbolic_matrix )
%SYM2STR transforms a symbolic matrix into a string 

sizler=size(symbolic_matrix);
string='[';
for i=1:sizler(1)
    new_line=char(symbolic_matrix(i,:));
    new_line=strrep(new_line,',','');
    string=[string,new_line(10:length(new_line)-3),';'];
end
string(length(string))=']';
end

