#!/bin/bash

# PREFIX definition work with an autotool generated Makefile
PREFIX=`cat Makefile | grep '^prefix =' | awk '{print $3}'`

case $1 in
    '')
        echo '-- building...'
        make -j3
        echo '-- installing...'
        make uninstall 1>/dev/null
        make install 1>/dev/null;;
    'clean')
        echo '-- cleaning...'
        make -j3 clean;;
    *)
        echo 'error: unknown action' 1>&2
        exit 1;;
esac
