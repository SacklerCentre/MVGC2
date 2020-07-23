#!/bin/sh

# $1 is version (eq 1.0), $2 is directory (eg /tmp)

# zip -r -v $2/mvgc_v$1.zip mvgc_v$1 -x mvgc_v$1/extra\* mvgc_v$1/testing\* mvgc_v$1/maintainer\*
zip -r -v $2/mvgc_v$1.zip mvgc_v$1 mvgc_v$1/testing\* mvgc_v$1/maintainer\*
