


c: build runc

python: runp

build: cmain.cpp
	g++ -Wall -Wextra -Wpedantic -o cmain.exe cmain.cpp

runc:
	./cmain.exe

runp:
	python nthroot-generator.py

all: c python