%1 -DWRITE_ARITH_H -DNO_FPINIT %2 %3 %4 %5 %6 %7 %8 %9
if errorlevel 1 goto nolonglong
exit 0
:nolonglong
%1 -DNO_LONG_LONG -DWRITE_ARITH_H -DNO_FPINIT %2 %3 %4 %5 %6 %7 %8 %9
