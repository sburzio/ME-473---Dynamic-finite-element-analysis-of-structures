@echo off
Echo -------------------------------------
Echo MAFE START, Problem = %1
date /t
time /t
Echo -------------------------------------
if exist %1.lis del %1.lis > NULL >> NULL
if exist %1.plt del %1.plt > NULL >> NULL
if exist %1.hpl del %1.hpl > NULL >> NULL
rem if exist %1.stv del %1.stv > NULL >> NULL
if exist %1.usr del %1.usr > NULL >> NULL
if exist %1.tmp del %1.tmp > NULL >> NULL
if exist %1.mat del %1.mat > NULL >> NULL
echo %1.dat > %1.ini
echo %1.lis >> %1.ini
echo %1.plt >> %1.ini
echo %1.hpl >> %1.ini
echo %1.stv >> %1.ini
echo %1.usr >> %1.ini
echo %1.tmp >> %1.ini
echo %1.mat >> %1.ini
echo %1.mod >> %1.ini
echo %1.mct >> %1.ini
echo lmahp6 (inter) >> %1.ini
rem -------------------------------------
rem Modify path to match executable location
C:\Emplacement\mafe.exe < %1.ini
Echo -------------------------------------
Echo MAFE END, Problem = %1
date /t
time /t
Echo -------------------------------------

