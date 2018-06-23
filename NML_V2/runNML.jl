# run NML
#= tags for comments
FIXME BUG TODO HACK NOTE REVIEW QUESTION
XXX: Person's name - ************
@FIXME: the code should be modified/refactored to achieve some goal (higher
maintainability, better performance, and so on)
@BUG - Indicates something where a situation on the input or flow of the program
leads to outright breaking of the program (freeze, crash or other severe problem)
or of the output (failing to output, or outputting wrong data).
@TODO - Indicates something that needs to be done, but that isn't a fix (may or
may not be important).
@HACK - Indicates a quick-fix or workaround that was done, either by necessity
or convenience, and that should be replaced or improved if possible.
@NOTE - Indicates any kind of note about something relevant enough to be
noteworthy, but that doesn't involve an action.
@REVIEW ???
@QUESTION
NOTE: good annotation --> the Function, what are the Input and Output (what
features the input should have, what features the output have), potential Bugs,
and the processing Approach. --> F I O B A
reduced the use of For-Loop, instead use comprehension (vector computation?)
=#

inputFilePath = "DATfiles/testInput.jl"
reservedWordsPath = "reservedWords.jl"
grammarFilePath = "DATfiles/grammar3.jl"

include("preprocessor3.jl")
println("\n------------loading sentences------------")
inputSentencesArray = getRidOfNewlineHashcommentBlankline(inputFilePath)
tokenizedInputSentencesArray = sentenceTokenization2(inputSentencesArray)
println("\n------------normalizing and tagging------------")
taggedSentencesArray = tokenClassification2(tokenizedInputSentencesArray, reservedWordsPath)
# reshape for output observation
printTagSen = [["$(y[1])/$(y[2])" for y in x] for x in taggedSentencesArray]
println(typeof(printTagSen))
foreach(println, [join(ts, "  ") for ts in printTagSen])


include("grammarGenerator3.jl")
GrammarTree = generateGrammarTree(grammarFilePath)
include("senRecognizer3.jl")
println("\n------------searching results-----------------------")
for (id, taggedSen) in enumerate(taggedSentencesArray)
  taggedSenTokens = [x[1] for x in taggedSen]
  taggedSenTags = [x[2] for x in taggedSen]
  searchSentenceInGrammarTree(GrammarTree, taggedSenTags, taggedSenTokens,
                                                  inputSentencesArray[id])
end

include("biosymDecoder.jl")
println("\n------------information extraction------------")
include(reservedWordsPath)
BioSymVerbInfo = extractBioSymVerbInformation(taggedSentencesArray, reservedWords["SentenceType"])
println("print by foreach():")
foreach(println, [s for s in BioSymVerbInfo])

println("\n------------decoding bio symbols------------")
decodingBioSymGroups(BioSymVerbInfo)
println("\n------------decoding bio symbols results------------")
printArrayOfTripleToken(BioSymVerbInfo, AbstractString)


# generate Intermediate Representation (IR)

# model generation

# code generation for different languages
