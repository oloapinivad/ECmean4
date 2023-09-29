# Changelog

All notable changes to ECmean4 will be documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)

## [Unreleased]

Unreleased is the current development version, which is currently lying in `main` branch.

- Remove support for python 3.8
- Refactor tests so that comparison is element-wise and not file-wise (Allow xesmf>0.8 and refactor the tests #93)
- Simple support for allow running the code on MacOS and pinning of xesmf<0.8 to avoid code slowdown (support for macos #86)
- Matplotlib pinning for solve temporary seaborn bug (Various updates and fixing the issue of matplotlib #89)

## [v0.1.5]

- Specifying output directory as an argument (Support more general output directory #87)

## [v0.1.4]

- Added the `changelog` file
- Better representation of mask file search and loading (#80)
- Fix of global mean temperature computaton (#84)
- Minor EC-Earth4 related fixes

## [v0.1.3]

- New logging facility
- Support for python 3.11

## [v0.1.0]

First version of ECmean4 published on pypi
