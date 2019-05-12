#!/bin/bash

$(dirname $0)/build.sh

git add -all
git push origin dev

git subtree push --prefix dist origin master
