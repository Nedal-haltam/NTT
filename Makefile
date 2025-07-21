


c: build runc

python: runp

build: cmain.cpp
	g++ -Wall -Wextra -Wpedantic -o cmain cmain.cpp

runc:
	./cmain

runp:
	python3 nthroot-generator.py

all: c python