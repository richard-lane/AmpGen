#!/bin/bash

# Must be run from within the repo
git remote add upstream https://github.com/GooFit/AmpGen.git
git fetch upstream
git checkout master
git merge upstream/master
git remote rm upstream

