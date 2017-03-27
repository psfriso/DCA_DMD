V = 0.2.6
EXEDIR =exe
#
.PHONY: lib install discrete
#
all:  discrete
	
lib: 
	cd lib; make

discrete: lib
	cd discrete; make

install: all
	cp discrete/discrete $(EXEDIR)/dcaDMD

clean:
	cd lib;make clean
	cd discrete; make clean
	rm $(EXEDIR)/dcaDMD

        	
