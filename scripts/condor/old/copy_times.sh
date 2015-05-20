#!/bin/bash

head -1 out/*.out | grep seconds | cut -d' ' -f3

