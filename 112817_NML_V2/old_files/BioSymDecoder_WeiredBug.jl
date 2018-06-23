
#= NOTE
F: extracting BioSym & Verb of each sentence;
I: each sentence in Array of tuples, sentenceType dictionary;
O: sentence type, & array of biosymbol groups;
BA
=#
function extractBioSymGroups(tokenTagPairs::Array, keyVerbs::Set)
  BioSymGroupsArray = []
  typeMarker = ""
  i = 1
  while i <= length(tokenTagPairs)
    if tokenTagPairs[i][2] == "BioSym"  # find a BioSym group
      tmpBioSymGroup = [tokenTagPairs[i][1]]
      i += 1  # go for next one
      # repeat until find the whole group
      while (i <= length(tokenTagPairs)) && (tokenTagPairs[i][2] == "BioSym")
        push!(tmpBioSymGroup, tokenTagPairs[i][1])
        i += 1  # go for next one
      end
      push!(BioSymGroupsArray, tmpBioSymGroup)
    elseif in(tokenTagPairs[i][2], keyVerbs)
      typeMarker = tokenTagPairs[i][2]
      i += 1  # go for next one
    else  # go for next one
      i += 1
    end
  end
  if typeMarker == ""
    println("ERROR: SENTENCE TYPE UNDEFINED for: $(tokenTagPairs)")
  end
  return typeMarker, BioSymGroupsArray
end

#= NOTE
F: extracting BioSym & Key Verb information from tagged sentences, call
   extractBioSymGroups() to process each sentence;
I: Array of sentence, with each sentence in form of array of tuples (token, tag);
O: Array of Tuples, with each tuple in form of (sentence type, array of biosymbol groups);
BA
=#
function extractBioSymVerbInformation(taggedSenteces::Array, sentenceTypeKeyVerbs::Set)
  BioSymGroupsKeyVerbArray = []
  for i in 1:length(taggedSenteces)
    senType, BioSym = extractBioSymGroups(taggedSenteces[i], sentenceTypeKeyVerbs)
    push!(BioSymGroupsKeyVerbArray, (senType, BioSym))
  end
  return BioSymGroupsKeyVerbArray
end



#= NOTE
F: decoding biosymbol groups by enumerating all possible combinations
I: Array of tuples (Verb, biosymbol groups), "biosymbol groups" is array of
   biosymbol arrays.
O: biosymbol groups are replaced with corresponding array of biosymbol array of
   biosymbol arrays. (2-layer array for simple case, 3-layer array for symbol
   groups with parenthesis)
BA
=#
function decodingBioSymGroups(VerbBioSymArray::Array)
  println(VerbBioSymArray)
  for vbsa in VerbBioSymArray
    for i in 1:length(vbsa[2])
      vbsa[2][i] = decodingABioSymGroupByLogicalRules(vbsa[2][i])
    end
  end
end



#= NOTE ambiguity problem
F: just handling simple cases (one-layer of parenthesis),
   e.g. (A or B) & C & D & (E or F), (A & B) or C or D.
   --> currently, try to do the right thing, but does not consider all possible cases.
I: array of biological tokens (including logcial tokens & bio-tokens)
O: Array of Strings for simple cases; or array of arrays of strings for cases
   with parenthesis.
BA
=#
function decodingABioSymGroupByLogicalRules(tokens::Array)
  println("\n\ninside decodingABioSymGroupByLogicalRules")
  println(tokens)
  andSet = Set(["and", "&", ","])
  orSet = Set(["or", "|"])

  finalBioSymGroups = []  # elements in parallel relation
  BioSymGroup = []  # for insiders, each element is in the form of [logicalRelation, token1, token2, ...]

  tokensOutParenthesis = []  # for the outsider
  i = 1
  while i <= length(tokens)
    if tokens[i] == "("  # find a ()
      tokensInsideParenthesis = []  # for tokens inside a pair of ()
      i += 1
      while (i<=length(tokens)) && (tokens[i] != ")") # assume one layer of ()
        push!(tokensInsideParenthesis, tokens[i])
        i += 1
      end
      # out of while loop, i is at ")"
      if length(tokensInsideParenthesis) != 0  #  has tokens insider () -> handle it
        logicInParenthesis, tmpBioSymGroup = decodingASequenceOfTokensWithOneLogicalRelation(
        tokensInsideParenthesis, andSet, orSet)
        if logicInParenthesis == ""
          logicInParenthesis = "and"
        end
        prepend!(tmpBioSymGroup, [logicInParenthesis])
        push!(BioSymGroup, tmpBioSymGroup)
      else  # what kind of errors?
        println("ERRORS: incomplete or empty \'()\' group")
      end
    else  # the outsiders
      push!(tokensOutParenthesis, tokens[i])
    end
    i += 1  # move forward
  end
  # handling outsiders
  logicOutParen, outsiderGroup = decodingASequenceOfTokensWithOneLogicalRelation(
  tokensOutParenthesis, andSet, orSet)
  println("outsiderGroup: $outsiderGroup")
  println("logicOutParen: $logicOutParen")
  if logicOutParen == ""  # no relation for outsiders
    if length(BioSymGroup) == 0   # only one token
      finalBioSymGroups = deepcopy(outsiderGroup)
      println("ONE token biosymbol")
    else
      println("ERROR: unspecified relationship")
    end
  else  # check logics in/outside () and reshape BioSymGroup for return
    goOnAllow =  true
    for (id, biosym) in enumerate(BioSymGroup)  # check logics in/outside ()
      if biosym[1] == logicOutParen  # logic check
        println("ERROR: logic inside and outside \'()\' should be different")
        goOnAllow =  false
      end
    end
    if goOnAllow  # reshape
      BioSymGroup = [subg[2:end] for subg in BioSymGroup]  # reshape
      reshapeBioSymbolGroups(finalBioSymGroups, BioSymGroup, outsiderGroup, logicOutParen)
      println("finalBioSymGroups : $finalBioSymGroups")
    end
  end
  println("\nfinalBioSymGroups: ")
  println(s for s in finalBioSymGroups)

  return finalBioSymGroups
end


#=
FIOBA
=#
function reshapeBioSymbolGroups(finalGroups, inParenGroups, outParenGroup, logicOutParen)
  if logicOutParen == "and"
    println("inParenGroups: $inParenGroups")
    println("outParenGroup: $outParenGroup")
    if length(inParenGroups) == 0  # simple case: Array of strings
      finalGroups = deepcopy(outParenGroup)
      println("finalGroups: $finalGroups")
      println("Logic AND simple bio-symbol group")
    else  # complicated case: array of arrays of strings
      finalGroups = [outParenGroup]
      # println("finalGroups: $finalGroups")
      for i = 1:length(inParenGroups)
        tmp2 = []
        for biosym in inParenGroups[i]  # append each token to the previous finalGroups
          tmp1 = deepcopy(finalGroups)
          # println("tmp1: $tmp1")
          for j = 1:length(tmp1)
            push!(tmp1[j], biosym)
          end
          # println("tmp1: $tmp1")
          append!(tmp2, tmp1)
        end
        # println("tmp2: $tmp2")
        finalGroups = deepcopy(tmp2)  # update finalGroups after processing each outParenGroup element
        # println("finalGroups: $finalGroups")
      end
      println("finalGroups: $finalGroups")
      println("Logic AND complicated bio-symbol group")
    end
  else  # "or" --> split outsiders into arrays, add insiders
    if length(inParenGroups) == 0 # simple case
      finalGroups = outsiderGroup
      println("Logic OR simple bio-symbol group")
    else
      for biosym in outParenGroup
        push!(finalGroups, [biosym])
      end
      for biosymgroup in inParenGroups
        push!(finalGroups, biosymgroup)
      end
      println("Logic OR complicated bio-symbol group")
    end
  end
  # return finalGroups
end
#= NOTE: You can change elements of the array, but you can't change the variable
so that it points to a different array. In other words, your function isn't
allowed to change the binding of the argument. --> passing by reference, means
inside the funtion, create a new pointer to the memory, if modification is allowed,
then change is saved, otherwise, take up a new memory space which is deleted after
the return, hence, modification is not saved. 
=#

#=
F: take in a sequence of tokens, seperate bio tokens and logical tokens, return
   logical relation and bio symbol array. Assume only one relation in this sequence.
IOBA
=#
function decodingASequenceOfTokensWithOneLogicalRelation(tokens::Array, andSet::Set, orSet::Set)
  # println("inside decodingASequenceOfTokensWithOneLogicalRelation")
  biosymGroup = []
  logicalRelation = [""]
  for i = 1:length(tokens)
    if in(tokens[i], andSet)  # find "and"
      if logicalRelation[1] == ""
        logicalRelation[1] = "and"
      elseif logicalRelation[1] == "or"
        println("LOGICAL INCONSISTENCY outside ()")
      else
        continue
      end
    elseif in(tokens[i], orSet)  # find "or"
      if logicalRelation[1] == ""
        logicalRelation[1] = "or"
      elseif logicalRelation[1] == "and"
        println("LOGICAL INCONSISTENCY outside ()")
      else
        continue
      end
    else
      push!(biosymGroup, tokens[i])
    end
  end
  # println(logicalRelation[1])
  # println("jump out decodingASequenceOfTokensWithOneLogicalRelation")
  return logicalRelation[1], biosymGroup
end





# NOTE test
include("preprocessor3.jl")
println("\n---------loading sentences------------")
inputSentencesArray = getRidOfNewlineHashcommentBlankline("DATfiles/testInput.jl")
tokenizedInputSentencesArray = sentenceTokenization2(inputSentencesArray)

println("\n------------normalizing and tagging------------")
reservedWordsPath = "reservedWords.jl"
taggedSentencesArray = tokenClassification2(tokenizedInputSentencesArray, reservedWordsPath)
# reshape for output observation
printTagSen = [["$(y[1])/$(y[2])" for y in x] for x in taggedSentencesArray]
println(typeof(printTagSen))
for ts in printTagSen
  println(join(ts, "  "))
end

println("\n------------information extraction------------")
include(reservedWordsPath)
BioSymVerbInfo = extractBioSymVerbInformation(taggedSentencesArray, reservedWords["SentenceType"])
for bsgkv in BioSymVerbInfo
  println(bsgkv)
end

println("\n------------decoding bio symbols------------")
decodingBioSymGroups(BioSymVerbInfo)
println("\n------------decoding bio symbols results------------")
# println(BioSymVerbInfo)
for tuple in BioSymVerbInfo
  println(tuple[1])
  println()
  if length(tuple[2]) == 0
    println("THIS sentence is wrong")
  else
    for symArray in tuple[2]
      if length(symArray) == 0
        println("THIS sentence is wrong")
      else
        if typeof(symArray[1]) <: AbstractString
          # println("print 1")
          println(symArray)
        else
          for symArray2 in symArray
            if typeof(symArray2[1]) <: AbstractString
              println("print 2")
              println(symArray2)
            else
              for symArray3 in symArray2
                if typeof(symArray3[1]) <: AbstractString
                  # println("print 3")
                  println(symArray3)
                else
                  for symArray4 in symArray3
                    if typeof(symArray4[1]) <: AbstractString
                      # println("NEED PRINT 4")
                      println(symArray4)
                    else
                      println(symArray4)
                      println("NEED PRINT 5")
                    end
                  end
                end
              end
            end
          end
        end
      end
      println("")
    end
  end
  println("-----------------")
end
