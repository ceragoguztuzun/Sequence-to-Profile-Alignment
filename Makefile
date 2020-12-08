all: bhw5

bhw5: bhw5.c
	gcc -Wall -o bhw5 bhw5.c -lm -O3 -mcmodel=medium

clean: 
	rm -fr *~ bhw5