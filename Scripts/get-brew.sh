#!/bin/bash

if ! type brew; then
    echo Do you want me to install Homebrew?
    echo Press RETURN to continue or any other key to abort
    read ans
    if test "$ans"; then
	echo Aborting
	exit 1
    else
	/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
    fi
fi
