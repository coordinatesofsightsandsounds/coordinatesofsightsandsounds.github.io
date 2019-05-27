#!/bin/bash
set -e

$(dirname $0)/build.sh

git add --all
npm version patch --force
git push origin dev --follow-tags

git subtree push --prefix dist origin master
