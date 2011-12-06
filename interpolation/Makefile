all: nninterp linterp tcinterp nninterp-ct linterp-ct tcinterp-ct lekien-tcinterp

nninterp: nninterp.c
	vcc nninterp.c -o nninterp
linterp: linterp.c
	vcc linterp.c -o linterp
tcinterp: tcinterp.c
	vcc tcinterp.c -o tcinterp
nninterp-ct: nninterp-ct.c
	vcc nninterp-ct.c -o nninterp-ct
linterp-ct: linterp-ct.c
	vcc linterp-ct.c -o linterp-ct
tcinterp-ct: tcinterp-ct.c
	vcc tcinterp-ct.c -o tcinterp-ct
lekien-tcinterp: lekien-tcinterp.c
	vcc lekien-tcinterp.c -o lekien-tcinterp

clean:
	rm nninterp linterp tcinterp nninterp-ct linterp-ct tcinterp-ct lekientcinterp -f
