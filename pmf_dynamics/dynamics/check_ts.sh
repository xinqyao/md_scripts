#!/bin/bash

if test $# -gt 0; then
   thres=$1
else
   thres=1
fi

grep -c " 1$" ts_time_life$thres/*
