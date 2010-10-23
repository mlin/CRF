.PHONY : all clean install uninstall reinstall

all:
	ocamlc -c -ccopt -O3 -ccopt -msse2 -ccopt -mfpmath=sse Elog_impl.c
	ocamlc -c Elog.mli
	ocamlmklib -o Elog Elog_impl.o Elog.ml
	
clean:
	rm -f *.cm* *.o *.a *.so
	
install: all
	ocamlfind install Elog META Elog.cmi Elog.cma Elog.cmxa dllElog.so Elog.a libElog.a
	
uninstall:
	ocamlfind remove Elog
	
reinstall:
	make uninstall
	make install
