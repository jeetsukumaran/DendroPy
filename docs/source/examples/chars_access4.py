import dendropy

dna = dendropy.DnaCharacterMatrix.get(
        path="primates.chars.nexus",
        schema="nexus")

s1 = dna[0]

print(type(s1))
# <class 'dendropy.datamodel.charmatrixmodel.DnaCharacterDataSequence'>

print(len(s1))
# 898

for v in s1:
    print("{}, {}".format(repr(v), str(v)))
# <<StateIdentity at 0x10134a290: 'A'>, A
# <<StateIdentity at 0x10134a290: 'A'>, A
# <<StateIdentity at 0x10134a350: 'G'>, G
# ...
# ...

print(s1.values())
# [<StateIdentity at 0x101b4a290: 'A'>, <StateIdentity at 0x101b4a290: 'A'>, <StateIdentity at 0x101b4a350: 'G'>, ...

print(s1.symbols_as_list())
# ['A', 'A', 'G', 'C', 'T', 'T', 'C', 'A', 'T', ...

print(s1.symbols_as_string())
# AAGCTTCATAGGAGCAACCATTCT ...


