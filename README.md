# Jigsaw Puzzle Solver

Solving Jigsaw Puzzles by using the basic Computer Vision techniques and linear programming method.

## Overview

A jigsaw puzzle is a tiling puzzle that requires the assembly of often oddly shaped interlocking and mosaiced pieces. Typically, each individual piece has a portion of a picture; when assembled, they produce a complete picture.

In this project, I have implemented an algorithm to solve the jigsaw puzzles.

Also, two papers are present in this repository that I have used to get some ideas on how to solve the jigsaw puzzles.

Stages of the algorithm include:
* Input patches
* Generate initial pairwise matches "U"
* Compute matrix "A" based on matrix "U"
* Run minimization on matrix "A"
* Sort the patches based on minimization
* Construct the resulting image

## Usage

First, run cvx_setup by using Matlab.

Then, just run main.m file by using Matlab.
