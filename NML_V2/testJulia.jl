type generalBioSym
  bioName::AbstractString
  bioType::AbstractString
  coeff::Float64  # Float64 has default value 0.0?
  function generalBioSym()
    new()
  end
end

# a = generalBioSym()
# println(a)
# println(isdefined(a, :bioName))
#
# println(isdefined(a, :bioType))
#
# println(isdefined(a, :coeff))
# println(a.coeff)

# b = ("a", "bc", "def")
# c = ("a", "def")
# println(isdefined(b, 3))
# println(isdefined(c, 3))

# testComArr = [["P_A", "P_B", "P_C"], ["P_E", "P_F"], ["P_H"]]
# testSimArr = ["P_A", "P_B", "P_C"]
#
# function processBioSymIntoGeneralBioSym(token::AbstractString)
#   bioSym = generalBioSym()
#   bioSym.bioName = token
#   bioSym.bioType = split(token, "_")[1]*"_"
#   return bioSym
# end
#
# function replaceAllSymbolIntoGeneralBioSym(testComArr::Array)
#   postTestComArr = []
#   if typeof(testComArr[1]) <: AbstractString
#     postTestComArr = [processBioSymIntoGeneralBioSym(s) for s in testComArr]
#   else
#     for arr in testComArr
#       if typeof(arr[1]) <: AbstractString
#         push!(postTestComArr, [processBioSymIntoGeneralBioSym(s) for s in arr])
#       else
#         println("Currently in land with no man")
#       end
#     end
#   end
#   return postTestComArr
# end
#
# postTestComArr = replaceAllSymbolIntoGeneralBioSym(testComArr)
# postTestSimArr = replaceAllSymbolIntoGeneralBioSym(testSimArr)

# function printArrayOfTupleOfTokenAndMutipleLayersArray(generalBioSymArray::Array, targetDataType::Type)
#   if typeof(generalBioSymArray[1]) <: targetDataType
#     println("print 1")
#     foreach(println, generalBioSymArray)
#   else
#     for arr1 in generalBioSymArray
#       if typeof(arr1[1]) <: targetDataType
#         println("print 2")
#         foreach(println, arr1)
#       else
#         println("world of no man")
#       end
#     end
#   end
#   println("---------------------")
# end

# printArrayOfTupleOfTokenAndMutipleLayersArray(postTestComArr, generalBioSym)
# printArrayOfTupleOfTokenAndMutipleLayersArray(postTestSimArr, generalBioSym)
