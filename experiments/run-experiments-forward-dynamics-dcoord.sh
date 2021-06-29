#!/usr/bin/env bash

set -e
set -x

endtime=5.0
BIN_DIR=../build/bin
export PATH=$BIN_DIR:$PATH

function part1()
{
# Part 1 ============================
# Effect of lag:
lags="0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009 0.010"
smootheriterations=15

for lag in $lags; do
  echo "Running: lag=$lag smoother-iterations=$smootheriterations..."
  time mbse-fg-smoother-forward-dynamics \
    --mechanism ../config/mechanisms/fourbars1.yaml \
      --dt 0.001 \
      --lag-time $lag \
      --end-time $endtime \
      --smoother-iterations $smootheriterations \
      --dont-show-error-progress \
      --output-prefix bars4_dc_dt=0.001_lag=${lag}_iters=${smootheriterations}_ #\
	  #> /dev/null
done
}

function part2()
{
# Part 2 ============================
# Effect of iterations
smootheriterationss="3 4 5 6 7 8 9 10"
lag=0.005

for smootheriterations in $smootheriterationss; do
  echo "Running: lag=$lag smoother-iterations=$smootheriterations..."
  mbse-fg-smoother-forward-dynamics \
      --dt 0.001 \
      --lag-time $lag \
      --mechanism 4bars \
      --end-time $endtime \
      --smoother-iterations $smootheriterations \
      --output-prefix bars4_dc_dt=0.001_lag=${lag}_iters=${smootheriterations}_ \
	  > /dev/null
done
}

part1
#part2
