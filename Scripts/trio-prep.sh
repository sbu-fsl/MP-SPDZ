#!/usr/bin/env bash

HERE=$(cd `dirname $0`; pwd)
SPDZROOT=$HERE/..

rm Player-Data/3-trio-*/*

export PLAYERS=3

. $HERE/run-common.sh

run_player trio-prep-party.x $* || exit 1
