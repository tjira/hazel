#!/bin/bash

cat include/* src/* tool/hview/include/* tool/hview/src/* tool/hplot/include/* tool/hplot/src/* | tr -d "[:blank:]" | sed "/^$/d ; /^\/\//d" | wc -l
