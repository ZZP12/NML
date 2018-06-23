# reserved words in Dict[normalized word] = Set{Synonyms}

# normalize synonyms
import JSON
# why need JSON?

# create dict in Dict
# JSON.json
# write
# JSON.parsefile <-- load the file

reservedWords = Dict{String, Any}()

reservedWords["induce"] = Set(["activate", "activates", "activating", "activated",
                            "promote", "promotes", "promoting", "promoted"])
reservedWords["repress"] = Set(["inhibit", "inhibits", "inhibiting", "inhibited",
                            "repress", "represses", "repressing", "repressed"])

reservedWords["BioSym"] = Set(["sBioSym", "mBioSym", "BioSym", "BioSymType1", "BioSymType2",
                            "BioSymType3"])

reservedWords["negligible"] = Set(["the", "expression", "of", "transcription", "translation",
                            ])

# save as JSON file
# Q: the writing format is very poor.
rWJSON = JSON.json(reservedWords)
fileName = "ReservedWords.json"
open(fileName, "w") do f
  write(f, rWJSON)
end
ReservedWords = JSON.parsefile(fileName)
