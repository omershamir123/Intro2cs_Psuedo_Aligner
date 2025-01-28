dict1 = {"a": 1, "b": 2}

print(next(iter(dict1))[0])
for item in dict1.items():
    print(item)