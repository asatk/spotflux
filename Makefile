CFLAGS = -Wall -Werror -Wextra -O3
LDFLAGS = -lm
CONSTS = src/constants.c
FILES = src/random.c src/linalg.c src/bmr.c src/flow.c src/field.c src/init.c src/methods.c src/io.c
DEPS = $(FILES) src/constants.h $(patsubst %.c,%.h, $(FILES))

main: $(CONSTS) $(DEPS)
	gcc $(CFLAGS) $(CONSTS) $(FILES) main.c $(LDFLAGS) -o $@

clean:
	rm -f main
