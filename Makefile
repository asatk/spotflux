CFLAGS = -Wall -Werror -g -O3
LDFLAGS = -lm

main: random.c bmr.c flow.c field.c init.c methods.c io.c main.c
	gcc $(CFLAGS) $^ $(LDFLAGS) -o $@

clean:
	rm main
