# grammar generator: nested dictionary


#= NOTE
F: grows a new branch in a grammar tree, may not be brand-new (see B);
I: target position, target sentence in tokens form, grammar tree;
O: grammar tree in nested dictionary --> position.name.children;
B: some branches may merge together, but seems serve our purpose very well;
   @HACK what happened previously was not "may merge together" but if "can merge
   together" occurs, new branch rewrite old one.
A: grow a new branch directly;
=#
function growABranch(id::Int, tokens::Array, gTree::Dict)
 if id < length(tokens) # not the end, add a branch
   if id in keys(gTree)
     if tokens[id] in keys(gTree[id])  # merge branches
       push!(gTree[id][tokens[id]], tokens[id+1])
     else  # new branches
       gTree[id][tokens[id]] = Set{AbstractString}([tokens[id+1]])
     end
     growABranch(id+1, tokens, gTree)
   else  # grow a branch from first key
     gTree[id] = Dict()
     gTree[id][tokens[id]] = Set{AbstractString}([tokens[id+1]])
     growABranch(id+1, tokens, gTree)
   end
 elseif id == length(tokens)  # the end
   if id in keys(gTree)
     if tokens[id] in keys(gTree[id])  # merge branches
       push!(gTree[id][tokens[id]], "END")
     else  # new branches
       gTree[id][tokens[id]] = Set{AbstractString}(["END"])
     end
   else  # grow a branch from first key
     gTree[id] = Dict()
     gTree[id][tokens[id]] = Set{AbstractString}(["END"])
   end
 else
   return
 end
end


#= NOTE
F: attach a sentence (in form of Array of tokens) to a grammar tree;
I: target position, target token, target sentence in tokens form, grammar tree;
O: grammar tree in nested dictionary --> position.name.children;
B: some branches may merge together, but seems serve our purpose very well;
A: try to match tokens recursively with the tree from the beginning, upon unmatched point,
   grow a new branch by calling growABranch;
=#
function attachToken2Tree(idx::Int, token::AbstractString, tokens::Array, gTree::Dict)
 if idx in keys(gTree)  # has postion
   if token in keys(gTree[idx]) # has token, add child
     if idx == length(tokens)  # reach the end
       push!(gTree[idx][token], "END")
       return    # end
     else  # add a child, go on attaching
       push!(gTree[idx][token], tokens[idx+1])
       attachToken2Tree(idx+1, tokens[idx+1], tokens, gTree)
     end
   else  # grow a branch
     gTree[idx][token] = Set{AbstractString}()
     if idx == length(tokens)  # reach the end
       push!(gTree[idx][token], "END")
     else   # add a child, then continue grow the branch
       push!(gTree[idx][token], tokens[idx+1])
       growABranch(idx+1, tokens, gTree)
     end
   end
 else  # add a position, then grow a branch
   gTree[idx] = Dict()
   gTree[idx][token] = Set{AbstractString}()
   if idx == length(tokens)  # reach the end
     push!(gTree[idx][token], "END")
   else   # add a child, then continue grow the branch
     push!(gTree[idx][token], tokens[idx+1])
     growABranch(idx+1, tokens, gTree)
   end
 end
end


#= NOTE preprocessing of grammar.dat can be merged with others
F: generate grammar tree;
I: file path;
O: grammar tree in nested dictionary --> position.name.children;
B: some branches may merge together, but seems serve our purpose very well;
A: read line-by-line, ingnore blank lines and pure comments, call attachToken2Tree;
=#
function generateGrammarTree(grammarFilePath)
  grammarTree = Dict{Int, Dict}()  # position.name.children
  # load grammar
  open(grammarFilePath, "r") do grammar
    for eachline in eachline(grammar)
      if !contains(eachline, "#") && eachline[1] != '\n'
        eachSentence = chomp(eachline)
        tokens = split(eachSentence)
        attachToken2Tree(1, tokens[1], tokens, grammarTree)
      end
    end
  end

  println("\n----------Grammar Tree-----------")
  for (key, val) in grammarTree
    for (key2, val2) in grammarTree[key]
      println(key, "==>", key2, "==>", val2)
    end
  end

  return grammarTree
end


# # XXX test
# generateGrammarTree("DATfiles/grammar3.jl")
