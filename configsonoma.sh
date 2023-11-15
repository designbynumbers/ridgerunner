#!/bin/zsh

./configure LDFLAGS="-L/usr/local/lib -L/opt/homebrew/lib -L/opt/homebrew/opt/openblas/lib" CPPFLAGS="-I/usr/local/include -I/opt/homebrew/include -I/opt/homebrew/opt/openblas/include"
