#!/bin/sh
set -eu
# Apply once per machine, never during shell startup.
defaults write -g InitialKeyRepeat -int 15
defaults write -g KeyRepeat -int 2
defaults write NSGlobalDomain ApplePressAndHoldEnabled -bool false
