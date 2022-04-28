#!/bin/bash

set -e

printf "\033[0;32mGenerating MNIST images figure\033[0m\n"

Rscript --vanilla mnist.R
