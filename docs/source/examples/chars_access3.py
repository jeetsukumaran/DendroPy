import dendropy

cc = dendropy.ContinuousCharacterMatrix.get(
        path="pythonidae_continuous.chars.nexml",
        schema="nexml")

s1 = cc[0]

print(type(s1))
# <class 'dendropy.datamodel.charmatrixmodel.ContinuousCharacterDataSequence'>

print(len(s1))
# 100

for v in s1:
    print("{}, {}".format(type(v), str(v)))
# <type 'float'>, -0.0230088801573
# <type 'float'>, -0.327376261257
# <type 'float'>, -0.483676644025
# ...
# ...

print(s1.values())
# [-0.0230088801573, -0.327376261257, -0.483676644025, ...

print(s1.symbols_as_list())
# ['-0.0230088801573', '-0.327376261257', '-0.483676644025', ...

print(s1.symbols_as_string())
# -0.0230088801573 -0.327376261257 -0.483676644025 0.0868649474847 ...


