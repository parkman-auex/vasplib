function OpenMX_bandplot(EIGENVAL)
arguments
    EIGENVAL double = []
end
if isempty(EIGENVAL)
    EIGENVAL = OpenMX_bandread();
end
bandplot(EIGENVAL)
end