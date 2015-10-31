
BIBFILES =
FIGS = $(wildcard src/figs/*)

SOURCE_DIR = src
BUILD_DIR = build
PRODUCT_DIR = ..

.PHONY: all html pdf clean

all: index.html archipelago.pdf

archipelago.pdf: src/archipelago.tex $(BIBFILES) $(FIGS)
	cd $(SOURCE_DIR) \
		&& latexmk -gg -halt-on-error -interaction=nonstopmode -file-line-error -pdf -use-make -output-directory=../$(BUILD_DIR) $(notdir $<) \
		&& cd .. \
		&& cp $(BUILD_DIR)/$(notdir $@) $(PRODUCT_DIR)/$(notdir $@)

index.html: archipelago.pdf
	pdf2htmlEX $(PRODUCT_DIR)/archipelago.pdf $(PRODUCT_DIR)/index.html

pdf: archipelago.pdf

html: index.html

clean:
	rm -rf build
	rm -f archipelago.pdf
	rm -f index.html

