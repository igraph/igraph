REM http://webcache.googleusercontent.com/search?q=cache:SjoPPpuQxuoJ:www.tcm.phy.cam.ac.uk/~mr349/cygwin_install.html+install+cygwin+ssh+commandline&cd=2&hl=nl&ct=clnk&gl=be&source=www.google.be

REM create the cygwin directory
cmd /c mkdir %SystemDrive%\cygwin

if "%PROCESSOR_ARCHITECTURE%" == "AMD64" (set ARCH=x86_64) else (set ARCH=x86)
set URL=http://cygwin.com/setup-%ARCH%.exe

cmd /c bitsadmin /transfer CygwinSetupExe /download /priority normal %URL% %SystemDrive%\cygwin\cygwin-setup.exe

REM goto a temp directory
cd /D %SystemDrive%\windows\temp

set PACKAGES= alternatives
set PACKAGES=%PACKAGES%,csih
set PACKAGES=%PACKAGES%,cygrunsrv
set PACKAGES=%PACKAGES%,crypt
set PACKAGES=%PACKAGES%,diffutils
set PACKAGES=%PACKAGES%,libasn1_8
set PACKAGES=%PACKAGES%,libattr1
set PACKAGES=%PACKAGES%,libcom_err2
set PACKAGES=%PACKAGES%,libcrypt0
set PACKAGES=%PACKAGES%,libffi6
set PACKAGES=%PACKAGES%,libgcc1
set PACKAGES=%PACKAGES%,libgcrypt11
set PACKAGES=%PACKAGES%,libgmp10
set PACKAGES=%PACKAGES%,libgmp3
set PACKAGES=%PACKAGES%,libgnutls26
set PACKAGES=%PACKAGES%,libgpg-error0
set PACKAGES=%PACKAGES%,libgssapi3
set PACKAGES=%PACKAGES%,libheimbase1
set PACKAGES=%PACKAGES%,libheimntlm0
set PACKAGES=%PACKAGES%,libhx509_5
set PACKAGES=%PACKAGES%,libiconv2
set PACKAGES=%PACKAGES%,libidn11
set PACKAGES=%PACKAGES%,libintl8
set PACKAGES=%PACKAGES%,libkafs0
set PACKAGES=%PACKAGES%,libkrb5_26
set PACKAGES=%PACKAGES%,libmpfr4
set PACKAGES=%PACKAGES%,libncursesw10
set PACKAGES=%PACKAGES%,libopenssl100
set PACKAGES=%PACKAGES%,libp11-kit0
set PACKAGES=%PACKAGES%,libpcre0
set PACKAGES=%PACKAGES%,libpcre1
set PACKAGES=%PACKAGES%,libreadline7
set PACKAGES=%PACKAGES%,libroken18
set PACKAGES=%PACKAGES%,libsqlite3_0
set PACKAGES=%PACKAGES%,libssp0
set PACKAGES=%PACKAGES%,libtasn1_3
set PACKAGES=%PACKAGES%,libwind0
set PACKAGES=%PACKAGES%,libwrap0
set PACKAGES=%PACKAGES%,openssh
set PACKAGES=%PACKAGES%,openssl
set PACKAGES=%PACKAGES%,rebase
set PACKAGES=%PACKAGES%,termcap
set PACKAGES=%PACKAGES%,terminfo
set PACKAGES=%PACKAGES%,wget
set PACKAGES=%PACKAGES%,zlib0

REM run the installation
%SystemDrive%\cygwin\cygwin-setup.exe -a %ARCH% -q -R %SystemDrive%\cygwin -P %PACKAGES% -s http://cygwin.mirrors.pair.com

REM stop the service, instead of attempting to remove it
%SystemDrive%\cygwin\bin\bash -c 'PATH=/usr/local/bin:/usr/bin:/bin:/usr/X11R6/bin cygrunsrv -E sshd'

REM /bin/ash is the right shell for this command
cmd /c %SystemDrive%\cygwin\bin\ash -c /bin/rebaseall

cmd /c %SystemDrive%\cygwin\bin\bash -c 'PATH=/usr/local/bin:/usr/bin:/bin:/usr/X11R6/bin mkgroup -l'>%SystemDrive%\cygwin\etc\group

cmd /c %SystemDrive%\cygwin\bin\bash -c 'PATH=/usr/local/bin:/usr/bin:/bin:/usr/X11R6/bin mkpasswd -l'>%SystemDrive%\cygwin\etc\passwd

%SystemDrive%\cygwin\bin\bash -c 'PATH=/usr/local/bin:/usr/bin:/bin:/usr/X11R6/bin /usr/bin/ssh-host-config -y -c "ntsecbinmode mintty" -w "abc&&123!!" '

cmd /c if exist %Systemroot%\system32\netsh.exe netsh advfirewall firewall add rule name="SSHD" dir=in action=allow program="%SystemDrive%\cygwin\usr\sbin\sshd.exe" enable=yes

cmd /c if exist %Systemroot%\system32\netsh.exe netsh advfirewall firewall add rule name="ssh" dir=in action=allow protocol=TCP localport=22

%SystemDrive%\cygwin\bin\bash -c 'PATH=/usr/local/bin:/usr/bin:/bin:/usr/X11R6/bin ln -s "$(/bin/dirname $(/bin/cygpath -D))" /home/$USERNAME'

net start sshd

REM Put local users home directories in the Windows Profiles directory
%SystemDrive%\cygwin\bin\bash -c 'PATH=/usr/local/bin:/usr/bin:/bin:/usr/X11R6/bin mkpasswd -l -p "$(/bin/cygpath -H)"'>%SystemDrive%\cygwin\etc\passwd

REM Fix corrupt recycle bin
REM http://www.winhelponline.com/blog/fix-corrupted-recycle-bin-windows-7-vista/
cmd /c rd /s /q %SystemDrive%\$Recycle.bin
