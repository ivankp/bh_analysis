DB := ntuples.db

.PHONY: all clean H AA remove_bad wt optimize

all: $(DB) H AA remove_bad optimize

$(DB):
	./init $@

H:
	./add_bh_Hj $(DB)

AA:
	./add_bh_AAj $(DB)

remove_bad:
	./remove_bad.sh

wt:
	./add_wt.py $(DB)

hist:
	./add_hist.py $(DB)

optimize:
	@sqlite3 $(DB) '.dump' > ~tmp.sql
	@sqlite3 ~tmp.db < ~tmp.sql
	@mv ~tmp.db $(DB)
	@rm -f ~tmp.sql ~tmp.db

clean:
	@rm -rv $(DB)

