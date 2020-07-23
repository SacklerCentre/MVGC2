#!/bin/sh

# $1 is version (eq 1.0), $2 is directory (eg /tmp)

# tar czvf $2/mvgc_v$1.tar.gz --exclude=extra --exclude=testing --exclude=maintainer mvgc_v$1/*
tar czvf $2/mvgc_v$1.tar.gz --exclude=testing --exclude=maintainer mvgc_v$1/*
