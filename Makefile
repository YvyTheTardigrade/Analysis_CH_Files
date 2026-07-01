.PHONY: run test build

run:
	PYTHONPATH=. python src/Main.py 


test:
	PYTHONPATH=. pytest	


build: clean
	pyinstaller --onefile --windowed --name "PepSynth_Analyzer" \
	--paths=src                      \
    --distpath "./build/executable/" \
    --workpath "./build/tmp/"        \
    --specpath "./build/spec/"       \
	src/Main.py

run_build:
	./build/executable/PepSynth_Analyzer 


clean :
	rm -rf build/tmp/PepSynth_Analyzer/
	#rm -f build/executable/PepSynth_Analyzer 
	#rm -f build/spec/PepSynth_Analyzer.spec



