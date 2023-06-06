#!/bin/bash

if [[ $# -eq 0 ]]; then
    echo "No basis name provided."; exit 1
fi

if [[ $# -eq 1 ]]; then
    echo "No save folder provided."; exit 1
fi

[[ $(bse get-basis "$1" gaussian94) ]] && bse get-basis "$1" gaussian94 > "$2/${1,,}.g94"
