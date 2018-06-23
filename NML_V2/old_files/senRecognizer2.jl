# grammar generator & sentence recognizer: self-defined data structure approach  

# Q? not only position & name, also related to its parent
# also data structure should be based on dictionary
keyReservedWords = Set{AbstractString}(["activate", "activates",
                   "activating", "activated","and", "or"])

unReservedWords = Set{AbstractString}(["the", "expression", "of"])

type GrammarToken
  parent::AbstractString
  name::AbstractString
  position::Int
  children::Set{AbstractString}

  function GrammarToken(name::AbstractString, position::Int, parent="START", children=Set{AbstractString}())
    new(parent, name, position, children)
  end
end

function getGT(a::Array{GrammarToken}, b::AbstractString, c::Int)
  for gt in a
    if b == gt.name && c == gt.position
      return gt
      break
    end
  end
  return false
end

grammarTree = GrammarToken[]
# load grammar
open("grammar2.dat", "r") do grammar
  for eachline in eachline(grammar)
    if !contains(eachline, "#") && eachline[1] != '\n'
      eachSentence = chomp(eachline)
      tokens = split(eachSentence)
      for (id, token) in enumerate(tokens)
        gt = getGT(grammarTree, token, id)
        if gt != false # already exit
          (id == length(tokens)) ? push!(gt.children, "END") : push!(gt.children, tokens[id+1])
        else  # create new one
          if id == 1
            gt = GrammarToken(token, id)
          else
            gt = GrammarToken(token, id, tokens[id-1])
          end
          push!(grammarTree, gt)
        end
      end
    end
  end
end
println("-------------------------")
for gt in grammarTree
  println(gt)
end
