#!/bin/bash

$(dirname $0)/build.sh

git subtree push --prefix dist origin master
