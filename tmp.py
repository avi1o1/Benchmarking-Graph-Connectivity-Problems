final = []
with open("./datasets/slashdot.txt", "r") as f:
    lines = f.readlines()
    
    for line in lines:
        final.append(" ".join(line.split()[:2]))

# print(final)

with open("./datasets/slashdot.txt", "w") as f:
    for line in final:
        f.write(line + "\n")

# 82167