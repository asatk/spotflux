CFLAGS = -Wall -Werror -g -O3
LDFLAGS = -lm

main: flow.c field.c bmr.c init.c methods.c io.c main.c
	gcc $(CFLAGS) $^ $(LDFLAGS) -o $@

clean:
	rm main
