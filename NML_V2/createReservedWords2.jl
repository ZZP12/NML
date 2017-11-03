
#= NOTE preprocessor3.jl also has this function
F: get rid of the trailing newlines, comments, blank lines in a file;
I: file path, read line-by-line;
O: array of strings, each element maps to a line in the file;
B: load the whole file at the same time, may require large memory for large file;
A: read line-by-line, ingnore blank lines and pure comments, throw away in-line comments;
=#
function getRidOfNewlineHashcommentBlankline(inputFilePath::String)
  f = open(inputFilePath, "r")
  lines = readlines(f)
  newlines = Array{String,1}()  # container for results
  for (id, val) in enumerate(lines)
    val = chomp(val)  # trailing newline
    hashId = searchindex(val, '#')
    if hashId != 0  # get rid of comments
      val = val[1:hashId-1]
    end
    if length(split(val)) != 0  # non-empty
      # throw away ugly trailing spaces, although can be done when splitting
      while(val[end] == ' ')
        val = val[1:end-1]
      end
      push!(newlines, val)
      println(val)
    end
  end
  println(typeof(newlines))
  return newlines
end


#= NOTE
F: create dictionary for synonyms and negligible words from a file, e.g. synonyms.dat;
I: file path, read line-by-line. In each line, the first token is the value,
   while all following tokens are keys. Delimiter between tokens can be equal sign
   or comma.
O: a Julia Dict structure mapping synonyms to the designated word, mapping
   negligible words to "negligible";
B: if there are repeated keys for a designated word, the output will have
   repeated lines; (while when generated jl file was imported, Julia compiler
   will automatically get rid of those repeated ones, so this is just an
   appearance issue)
A: read line-by-line, split by "=|,", and then go thru each token;
=#
function createReservedWordsDict(inputFilePath::String, outputPath::String)
  synSens = getRidOfNewlineHashcommentBlankline(inputFilePath)
  buffer = "\nreservedWords = Dict{String, Any}()"
  reservedWordsForSentenceType =  false
  for sen in synSens
    sen = replace(sen, r"=|,", " ")
    tokens = split(sen)
    println(tokens)
    if tokens[1] == "SentenceType"
      if !reservedWordsForSentenceType
        buffer *= "\nreservedWords[\"$(tokens[1])\"] = Set{AbstractString}()"
        reservedWordsForSentenceType = true
      end
      for i in 2:length(tokens)
        buffer *= "\npush!(reservedWords[\"$(tokens[1])\"], \"$(tokens[i])\")"
      end
    else
      for i in 2:length(tokens)
        buffer *= "\nreservedWords[\"$(tokens[i])\"] = \"$(tokens[1])\""
      end
    end
  end
  write(outputPath, buffer)
end


# NOTE test
createReservedWordsDict("DATfiles/synonyms3.dat", "reservedWords.jl")
