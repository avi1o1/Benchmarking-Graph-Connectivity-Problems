import csv

with open("large_twitch_edges.csv", "r") as f:
    reader = csv.reader(f)
    next(reader)
    for row in reader:
        with open("datasets/large.txt", "a") as f:
            f.write(f"{row[0]} {row[1]}\n")