# vim: ft=make
.PHONY: RichardKirchofer._graphics
RichardKirchofer.aux RichardKirchofer.aux.make RichardKirchofer.d RichardKirchofer.pdf: $(call path-norm,/usr/share/texlive/texmf-dist/tex/latex/base/article.cls)
RichardKirchofer.aux RichardKirchofer.aux.make RichardKirchofer.d RichardKirchofer.pdf: $(call path-norm,RichardKirchofer.tex)
.SECONDEXPANSION:
