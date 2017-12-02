RM=rm -f
FORT=gfortran
FORT_FLAGS=

test: marspack.o test_mars.f03
	$(FORT) $(FORT_FLAGS) -o test_mars marspack.o test_mars.f03
	./test_mars

marspack.o: marspack.f
	$(FORT) -c $(FORT_FLAGS) marspack.f

clean:
	$(RM) *.o test_mars

