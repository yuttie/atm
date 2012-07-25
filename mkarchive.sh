#!/bin/sh

version=$1
format=$2
appname=`sed -n -e "/^APPNAME *= */ { s/^APPNAME *= *\(['\"]\)\([^\1]*\)\1$/\2/; p }" wscript`
tag=${appname}-${version}
cmdline="git archive --prefix=${appname}-${version}/ -o ${appname}-${version}.${format} ${tag}"
echo $cmdline
$cmdline
