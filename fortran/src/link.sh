#!/bin/bash

#Script to create a link for the main executable on the main directory
version=` date +%y.%m.%d `
binname="bin/shw$version"
echo $binname
cp bin/shw $binname
rm -rf shw
if [[ -f shw ]] ; then
    ln -n -s -f $binname shw
    echo "A link for the main executable (named 'shw') was updated " $binname  
else
    if [[ -f bin/shw ]] ; then
	ln -s $binname shw
	echo "A link for the main executable (named 'shw') was created " $binname 
    else
	echo "Could not create link, bin/shw does not exist" $binname 
    fi
fi
echo
