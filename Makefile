CFLAGS = -Wall -Werror -g -O3
LDFLAGS = -lm
CONSTS = src/constants.c
FILES = src/random.c src/linalg.c src/bmr.c src/flow.c src/field.c src/init.c src/methods.c src/io.c main.c

main: $(CONSTS) $(FILES)
	gcc $(CFLAGS) $^ $(LDFLAGS) -o $@

clean:
	rm main
