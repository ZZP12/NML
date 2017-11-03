# General Framework

synonyms.dat --> Reserved Words Creator --> reservedWords.jl
preprocessorTest.dat --> Preprocessor (= tokenizer + lemmatizer)

grammar.dat --> Grammar Generator
senRecognizerTest.dat --> Sentence Recognizer

BioSym Recognizer + Key Verb Extractor
  1) general BioSym recognizer? --> deal unique features in IR generation?


>>>Semantic checking (for each type of sentence)
  -> symbol conversion dictionary
  -> decomposing into general bio-symbols & syntactic checking for each symbol (type is in defined types)
---> semantic checking for each sentence type
      1) biologically meaningful
      2) have all components for each type, if not, add default component
      3) Ready for IR generation
  ? decomposing for checking is necessary?

>>>Intermediate Representation (IR) Generator



NOTE:
1) keep a good annotation: full description of every functions --> FIOBA
2) encourage the use of functions
3) keep a good code structure
all is about logic flow and error checking!!!
