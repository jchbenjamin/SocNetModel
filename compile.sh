#!/bin/bash

g++ game.cpp -c -g -std=c++0x -rdynamic
g++ gene.cpp -c -g -std=c++0x -rdynamic
g++ generation.cpp -c -g -std=c++0x -rdynamic
g++ node.cpp -c -g -std=c++0x -rdynamic
g++ util.cpp -c -g -std=c++0x -rdynamic

g++ -o game game.o gene.o generation.o node.o util.o -rdynamic -std=c++0x

g++ cycletest.cpp -c -g -std=c++0x -rdynamic

g++ -o cycletest cycletest.o gene.o generation.o node.o util.o -rdynamic -std=c++0x
