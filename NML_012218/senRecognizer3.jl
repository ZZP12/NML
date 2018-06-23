
#= NOTE
F: throw an error message.
I: error type, position of error token in target sentence, tokens of target sentence.
O: message 1 if error found at the 1st token, otherwise, message 2.
B: the purpose of arguments "tokens" and "sentence"? just for error message,
   hence, repeated inputs here.
A:
=#
function generateUnkownSentenceError(errorType::AbstractString, errorPoint::Int,
  tokens::Array, sentence::AbstractString)
  if errorPoint == 0
    println(errorType * " found from the beginning of sentence \""
        * sentence * "\"")
  else
    println(errorType * " found after \"" * tokens[errorPoint] * "\" in sentence \""
      * sentence * "\"")
  end
end


#= NOTE
F: recursively search a token thru a grammar tree.
I: grammar tree in nested Dict, target position (level) in the tree, tokens to
   be matched, tokens for error message.
O: "SUCCESS" if success, otherwise throw out corresponding error message.
B: names of function arguments are inappropriate, should be more general, e.g.
   searchToken, tagForError. Repeated inputs (arguments) here
A: searching step-by-step w/ alternative breadth and depth first search (BF-WF
   alternative search), call generateUnkownSentenceError if encountering errors,
   call handleBioSymToken if not the end, print "SUCCESS" if success.
=#
function searchForNextToken(gTree::Dict, nextId::Int, tokens::Array,
  tagForError::Array, sentence::AbstractString)
  if nextId <= length(tokens)  # search
    if in(tokens[nextId], gTree[nextId-1][tokens[nextId-1]])  # in the children list
      if nextId == length(tokens)  # reach the end of sentence -> check the end of the branch?
        if in("END", gTree[nextId][tokens[nextId]])  # confirmed
          println("SUCCESS")
        else  # error
          generateUnkownSentenceError("IMCOMPLETE SENTENCE", nextId, tagForError, sentence)
        end
      else  # continue searching: check BioSym
        # TODO (reorganize) BioSym? handleBioSymToken : searchForNextToken
        handleBioSymToken(gTree, nextId+1, tokens, tagForError, sentence)
      end
    else   #  error
      generateUnkownSentenceError("GRMMAR ERROR", nextId-1, tagForError, sentence)
    end
  else  # reach the end of sentence in last token, need confirmation
    if in("END", gTree[nextId-1][tokens[nextId-1]])  # confirmed
      println("SUCCESS")
    else  # error
      generateUnkownSentenceError("IMCOMPLETE SENTENCE", nextId-1, tagForError, sentence)
    end
  end
end


#= NOTE
F: recursively processing BioSym token(s)
I: grammar tree in nested Dict, target position (level) in the tree, tokens to
   be matched, tokens for error message, original sentence.
O: calling searchForNextToken if not the end, otherwise, print "SUCCESS" or throw
   error message.
B: .
A: if the previous token is BioSym, enter "recursively BioSym matching" process,
   remove extra BioSyms in a sequence away one-by-one until reach next non-biosym
   or the end. E.g. "BioSym BioSym BioSym" -> "BioSym"; "BioSym BioSym induce
   BioSym BioSym" -> "BioSym induce BioSym"
=#
function handleBioSymToken(gTree::Dict, nextId::Int, tokens::Array,
  tagForError::Array, sentence::AbstractString)
  if nextId <= length(tokens)
    if tokens[nextId-1] == "BioSym"  # try recursively matching for BioSym
      while (nextId <= length(tokens) && tokens[nextId] == "BioSym")  # remove BioSym from tokens
        splice!(tokens, nextId)
        splice!(tagForError, nextId)
      end
    end
    # continue searching
    searchForNextToken(gTree, nextId, tokens, tagForError, sentence)
  else  # reach the end of sentence in last token, need confirmation
    if in("END", gTree[nextId-1][tokens[nextId-1]])  # confirmed
      println("SUCCESS")
    else  # error
      generateUnkownSentenceError("IMCOMPLETE SENTENCE", nextId, tagForError, sentence)
    end
  end
end


#= NOTE
F: recursively search a token thru a grammar tree.
I: grammar tree in nested Dict, tags of tokens to be matched, sentence in tokens
   (for error message).
O: "SUCCESS" if success, otherwise throw out corresponding error message.
B: the purpose of arguments "senTokens" and "sentence"? just for error message,
   hence, repeated inputs here.
A: searching step-by-step w/ alternative breadth and depth first search (BF-WF
   alternative search), calling searchForNextToken if no error occurs, calling
   generateUnkownSentenceError if encountering errors.
=#
function searchSentenceInGrammarTree(gTree::Dict, senTags::Array,
                        senTokens::Array, sentence::AbstractString)
  if senTags[1] in keys(gTree[1])
    handleBioSymToken(gTree, 2, senTags, senTokens, sentence)
  else
    generateUnkownSentenceError("GRMMAR ERROR", 0, senTokens, sentence)
  end
end


# # NOTE: test searching method
# include("grammarGenerator3.jl")
# GrammarTree = generateGrammarTree("DATfiles/grammar3.dat")
# open("DATfiles/senRecognizerTest.dat", "r") do test3
#   println("\n----------test results-----------")
#   for eachline in eachline(test3)
#     if !contains(eachline, "#") && eachline[1] != '\n'
#       eachTest = chomp(eachline)
#       tokens = split(eachTest)
#       if length(tokens) == 0  # get rid of lines with only spaces
#         continue
#       end
#       # start searching
#       searchSentenceInGrammarTree(GrammarTree, tokens, tokens, eachTest)
#     end
#   end
# end
