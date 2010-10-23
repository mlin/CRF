all: Elog CRF

.PHONY: Elog CRF

Elog:
	cd lib/Elog; $(MAKE) $(MFLAGS) reinstall
	
CRF:
	cd lib/CRF; $(MAKE) $(MFLAGS) reinstall

clean:
	cd lib/Elog; $(MAKE) clean
	cd lib/CRF; $(MAKE) clean
