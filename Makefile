CFLAGS = -Wall -Werror -g -O2
LDFLAGS = -lm

main: flow.c field.c init.c methods.c io.c main.c
	gcc $(CFLAGS) $^ $(LDFLAGS) -o $@

clean:
	rm main *.o *~ core.[1-9]*
