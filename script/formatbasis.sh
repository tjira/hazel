#!/bin/bash

if [[ $# -eq 0 ]]; then
    echo "No basis file provided."; exit 1
fi

cat "$1" | sed 's/D-/E-/g ; s/D+/E+/g' | while read LINE; do
    if [[ "$LINE" == !* ]]; then echo "$LINE"
    elif [[ -z "$LINE" ]]; then continue
    elif [[ "$LINE" == [0-9]* ]]; then
        if [[ $(echo "$LINE" | awk '{print NF}') == 6 ]]; then
            echo "$LINE" | awk '{printf "% 14.10e % 14.10e % 14.10e % 14.10e % 14.10e % 14.10e\n", $1, $2, $3, $4, $5, $6}'
        elif [[ $(echo "$LINE" | awk '{print NF}') == 5 ]]; then
            echo "$LINE" | awk '{printf "% 14.10e % 14.10e % 14.10e % 14.10e % 14.10e\n", $1, $2, $3, $4, $5}'
        elif [[ $(echo "$LINE" | awk '{print NF}') == 4 ]]; then
            echo "$LINE" | awk '{printf "% 14.10e % 14.10e % 14.10e % 14.10e\n", $1, $2, $3, $4}'
        elif [[ $(echo "$LINE" | awk '{print NF}') == 3 ]]; then
            echo "$LINE" | awk '{printf "% 14.10e % 14.10e % 14.10e\n", $1, $2, $3}'
        elif [[ $(echo "$LINE" | awk '{print NF}') == 2 ]]; then
            echo "$LINE" | awk '{printf "% 14.10e % 14.10e\n", $1, $2}'
        else
            echo "$LINE"
        fi
    else
        if [[ $(echo "$LINE" | awk '{print NF}') == 3 ]]; then
            echo "$LINE" | awk '{printf "%s %d %3.2f\n", $1, $2, $3}'
        elif [[ $(echo "$LINE" | awk '{print NF}') == 2 ]]; then
            echo "$LINE" | awk '{printf "%s %d\n", $1, $2}'
        else
            echo "$LINE"
        fi
    fi
done > "$1.tmp" && mv "$1.tmp" "$1"
