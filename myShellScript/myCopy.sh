#!/bin/bash
if [ $# -ne 1 ]; then
    echo $0: usage: myCopy name
    exit 1
fi

NAME=$1
LOC=$NAME"_ha"
ADDON="_pl"
cp -r ../July-22/$LOC/$NAME/ .
mkdir "$NAME$ADDON"
mv $NAME "$NAME$ADDON"
