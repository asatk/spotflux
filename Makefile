CFLAGS = -Wall -Werror -g -O3
LDFLAGS = -lm

main: constants.c random.c linalg.c bmr.c flow.c field.c init.c methods.c io.c main.c
	gcc $(CFLAGS) $^ $(LDFLAGS) -o $@

clean:
	rm main
