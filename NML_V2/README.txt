# General Framework

synonyms.dat --> Reserved Words Creator (createReservedWords.jl) --> reservedWords.jl
preprocessorTest.dat --> Preprocessor (= tokenizer + lemmatizer)

grammar.dat --> Grammar Generator (does not use grammar in the prototype)
senRecognizerTest.dat --> Sentence Recognizer

BioSym Recognizer + Key Verb Extractor
  1) general BioSym recognizer? --> deal unique features in IR generation?


Semantic checking (for each type of sentence)
  -> symbol conversion dictionary
  -> decomposing into general bio-symbols & syntactic checking for each symbol
     (type is in defined types)
  -> semantic checking for each sentence type
      1) biologically meaningful
      2) have all components for each type, if not, add default component
      3) Ready for IR generation
  ? decomposing for checking is necessary? for convenience?


Intermediate Representation (IR) Generator (based on previous version,
try to make this part reusable)
  -> generate RNX-list, RNX-species, txtl-list, txtl-species
  -> sort all species by type, set up Dictionary of species and position
  -> calling Julia strategy to generate Stoichiometric Table, Kinetics,
     Data Dictionary, Simulation
  -> DONE!


  1. check julia strategy:
    1) add biomass-growth
  2. extend to other language: python3 python2 and MATLAB
  3. clear up pre-IR part:
    1) improving error checking system and collecting results -> error report
    2) find out potential bugs <-- needs to learn more about compiler, not an easy task
-->
  4. what are default settings?

1. change file names
2. run washout ??? 

NOTE:
1) keep a good annotation: full description of every functions --> FIOBA
2) encourage the use of functions
3) keep a good code structure
all is about logic flow and error checking!!!
