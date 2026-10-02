.PHONY: install test figures figures-rapide sql interactif
install:        ; pip install -r requirements.txt
test:           ; python -m pytest -q
figures:        ; python -m magnetisme figures
figures-rapide: ; python -m magnetisme figures --rapide
sql:            ; python -m magnetisme sql
interactif:     ; python -m magnetisme interactif
