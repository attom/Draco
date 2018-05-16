@echo off
rem ---------------------------------------------------------------------------
rem File  : regression/update_regrssion_dir.bat
rem Date  : Tuesday, May 31, 2016, 14:48 pm
rem Author: Kelly Thompson
rem Note  : Copyright (C) 2016, Los Alamos National Security, LLC.
rem         All rights are reserved.
rem ---------------------------------------------------------------------------

set LOGDIR=e:\regress\logs
set GIT="c:\\Program Files\\Git\\bin\\git.exe"
set PATH="c:\\Program Files\\Git\\bin";%PATH%
echo cd e:\regress\draco
cd /d e:\regress\draco
echo %GIT% pull
%GIT% pull
echo cd e:\regress\jayenne
cd /d e:\regress\jayenne
echo %GIT% pull
%GIT% pull

rem set SVN_SSH="c:\\Program Files\\TortoiseSVN\\bin\\TortoisePlink.exe"
rem set SVN="c:\\Program Files\\TortoiseSVN\\bin\\svn.exe"
rem %SVN% update > %LOGDIR%\update_regression_dir.log

rem 2&>1
