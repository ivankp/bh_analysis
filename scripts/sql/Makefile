DB := ntuples.db

.PHONY: all H AA clean

all: $(DB) H AA

$(DB):
	./init $@

H:
	./add_bh_Hj $(DB)

AA:
	./add_bh_AAj $(DB)

clean:
	@rm -rv $(DB)

