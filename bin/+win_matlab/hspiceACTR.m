function hspiceACTR(filename)
copyfile(filename,strcat(filename,".bk"));
Selectline = win_matlab.FindStrLines(filename,"$&%#");
%
PowerShellCmd = "powershell -command ""type " + filename + "|Select -First "+...
    string(Selectline) + " > "+ filename + ".info"" ";
system(PowerShellCmd);
%
PowerShellCmd = "powershell -command ""Get-Content " + filename + "|Select-Object -Skip "+...
    string(Selectline) + " |Set-Content "+ filename + ".data "" ";
system(PowerShellCmd);
end