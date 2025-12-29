all:
	odin build . -out:planet_orbits.exe -o:speed

clean:
	rm -f ./planet_orbits.exe

run:
	./planet_orbits.exe
