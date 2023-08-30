#!/bin/bash

cat include/* src/* | tr -d "[:blank:]" | sed "/^$/d ; /^\/\//d" | wc -l
